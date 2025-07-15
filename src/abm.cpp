#include "abm.h"
#pragma omp declare reduction(merge_int_pair_vecs : std::vector<std::pair<int, int>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(merge_str_int_pair_vecs : std::vector<std::pair<std::string, int>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(merge_int_vecs : std::vector<int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

int ABM::WriteToLogFile(std::string message, Log message_type) {
    if(this->log_level >= message_type) {
        std::chrono::time_point<std::chrono::steady_clock> now = std::chrono::steady_clock::now();
        std::string log_message_prefix;
        if(message_type == Log::info) {
            log_message_prefix = "[INFO]";
        } else if(message_type == Log::debug) {
            log_message_prefix = "[DEBUG]";
        } else if(message_type == Log::error) {
            log_message_prefix = "[ERROR]";
        }
        auto days_elapsed = std::chrono::duration_cast<std::chrono::days>(now - this->start_time);
        auto hours_elapsed = std::chrono::duration_cast<std::chrono::hours>(now - this->start_time - days_elapsed);
        auto minutes_elapsed = std::chrono::duration_cast<std::chrono::minutes>(now - this->start_time - days_elapsed - hours_elapsed);
        auto seconds_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - this->start_time - days_elapsed - hours_elapsed - minutes_elapsed);
        auto total_seconds_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - this->start_time);
        log_message_prefix += "[";
        log_message_prefix += std::to_string(days_elapsed.count());
        log_message_prefix += "-";
        log_message_prefix += std::to_string(hours_elapsed.count());
        log_message_prefix += ":";
        log_message_prefix += std::to_string(minutes_elapsed.count());
        log_message_prefix += ":";
        log_message_prefix += std::to_string(seconds_elapsed.count());
        log_message_prefix += "]";

        log_message_prefix += "(t=";
        log_message_prefix += std::to_string(total_seconds_elapsed.count());
        log_message_prefix += "s)";
        this->log_file_handle << log_message_prefix << " " << message << '\n';

        if(this->num_calls_to_log_write % 1 == 0) {
            std::flush(this->log_file_handle);
        }
        this->num_calls_to_log_write ++;
    }
    return 0;
}

void ABM::ReadPlantedNodes() {
    char delimiter = ',';
    std::ifstream planted_nodes_stream(this->planted_nodes);
    std::string line;
    int line_no = 1;
    while(std::getline(planted_nodes_stream, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        std::string year = current_line[0];
        std::string fitness_lag_duration = current_line[1];
        std::string fitness_peak_value = current_line[2];
        std::string fitness_peak_duration = current_line[3];
        std::string count = current_line[4];
        this->planted_nodes_map[std::stoi(year)][line_no]["fitness_lag_duration"] = std::stoi(fitness_lag_duration);
        this->planted_nodes_map[std::stoi(year)][line_no]["fitness_peak_value"] = std::stoi(fitness_peak_value);
        this->planted_nodes_map[std::stoi(year)][line_no]["fitness_peak_duration"] = std::stoi(fitness_peak_duration);
        this->planted_nodes_map[std::stoi(year)][line_no]["count"] = std::stoi(count);
        line_no += 1;
        /* std::cout << "adding " << current_line[1] << " to  bag " << std::endl; */
    }
}

void ABM::ReadOutDegreeBag() {
    char delimiter = ',';
    std::ifstream out_degree_bag_stream(this->out_degree_bag);
    std::string line;
    int line_no = 0;
    while(std::getline(out_degree_bag_stream, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        if(line_no != 0) {
            this->out_degree_bag_vec.push_back(std::stoi(current_line[1]));
        }
        line_no ++;
        /* std::cout << "adding " << current_line[1] << " to  bag " << std::endl; */
    }
}


void ABM::ReadRecencyProbabilities() {
    char delimiter = ',';
    std::ifstream recency_probabilities_stream(this->recency_probabilities);
    std::string line;
    int line_no = 0;
    while(std::getline(recency_probabilities_stream, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        if(line_no != 0) {
            int integer_year_diff = std::stoi(current_line[0]);
            double probability = std::stod(current_line[1]);
            this->recency_probabilities_map[integer_year_diff] = probability;
        }
        line_no ++;
    }
}

std::unordered_map<int, int> ABM::BuildContinuousNodeMapping(Graph* graph) {
    int next_node_id = 0;
    std::unordered_map<int, int> continuous_node_mapping;
    for(auto const& node : graph->GetNodeSet()) {
        continuous_node_mapping[node] = next_node_id;
        next_node_id ++;
    }
    return continuous_node_mapping;
}

std::unordered_map<int, int> ABM::ReverseMapping(std::unordered_map<int, int> mapping) {
    std::unordered_map<int, int> reverse_mapping;
    for(auto const& [key,val] : mapping) {
        reverse_mapping[val] = key;
    }
    return reverse_mapping;
}


void ABM::FillInDegreeArr(Graph* graph, const std::unordered_map<int, int>& continuous_node_mapping, int* in_degree_arr) {
    for(auto const& node: graph->GetNodeSet()) {
        int continuous_id = continuous_node_mapping.at(node);
        in_degree_arr[continuous_id] = graph->GetInDegree(node);
    }
}

void ABM::InitializeFitness(Graph* graph) {
    this->AssignPeakFitnessValues(graph, graph->GetNodeSet());
    this->AssignFitnessLagDuration(graph, graph->GetNodeSet());
    this->AssignFitnessPeakDuration(graph, graph->GetNodeSet());
}

void ABM::FillFitnessArr(Graph* graph, const std::unordered_map<int, int>& continuous_node_mapping, int current_year, int* fitness_arr) {
    for(auto const& node : graph->GetNodeSet()) {
        int fitness_peak_value = graph->GetIntAttribute("fitness_peak_value", node);
        int fitness_lag_duration = graph->GetIntAttribute("fitness_lag_duration", node);
        int fitness_peak_duration = graph->GetIntAttribute("fitness_peak_duration", node);
        int published_year = graph->GetIntAttribute("year", node);
        int continuous_index = continuous_node_mapping.at(node);
        if (published_year + fitness_lag_duration > current_year) {
            fitness_arr[continuous_index] = 1;
        } else if (published_year + fitness_lag_duration + fitness_peak_duration >= current_year) {
            fitness_arr[continuous_index] = fitness_peak_value;
        } else {
            double decayed_fitness_value = fitness_peak_value / pow(current_year - published_year - fitness_lag_duration - fitness_peak_duration + 1, this->fitness_decay_alpha);
            fitness_arr[continuous_index] = decayed_fitness_value;
        }
    }
}


void ABM::FillRecencyArr(Graph* graph, const std::unordered_map<int, int>& reverse_continuous_node_mapping, int current_year, double* recency_arr) {
    std::unordered_map<int, int> year_count;
    double unique_year_sum = 0.0;
    for(auto const& node : graph->GetNodeSet()) {
        int current_published_year = graph->GetIntAttribute("year", node);
        int year_diff = current_year - current_published_year;
        if(!year_count.contains(year_diff)) {
            unique_year_sum += this->recency_probabilities_map[year_diff];
        }
        year_count[year_diff] ++;
    }
    // Mark: removed for node-level
    #pragma omp parallel for
    for(size_t i = 0; i < graph->GetNodeSet().size(); i ++) {
        int node_id = reverse_continuous_node_mapping.at(i);
        int current_published_year = graph->GetIntAttribute("year", node_id);
        int year_diff = current_year - current_published_year;
        recency_arr[i] = (float)this->recency_probabilities_map[year_diff] / year_count[year_diff];
    }
    // Mark: removed for node-level
    #pragma omp parallel for
    for(size_t i = 0; i < graph->GetNodeSet().size(); i ++) {
        recency_arr[i] /= unique_year_sum;
    }
}


int ABM::GetMaxYear(Graph* graph) {
    int max_year = -1;
    bool is_first = true;
    for(auto const& node : graph->GetNodeSet()) {
        int current_node_year = graph->GetIntAttribute("year", node);
        if (is_first) {
            max_year = current_node_year;
            is_first = false;
        }
        if (current_node_year > max_year) {
            max_year = current_node_year;
        }
    }
    return max_year;
}

int ABM::GetMaxNode(Graph* graph) {
    int max_node = -1;
    bool is_first = true;
    for(auto const& node : graph->GetNodeSet()) {
        if (is_first) {
            max_node = node;
            is_first = false;
        }
        if (node > max_node) {
            max_node = node;
        }
    }
    return max_node;
}

int ABM::GetFinalGraphSize(Graph* graph) {
    int current_graph_size = graph->GetNodeSet().size();
    for(int i = 0; i < this->num_cycles; i ++) {
        int num_new_nodes = std::ceil(current_graph_size * this->growth_rate);
        current_graph_size += num_new_nodes;
    }
    return current_graph_size;
}
void ABM::PopulateAlphaArr(double* alpha_arr, int len) {
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);

    if(this->alpha < 0) {
        for(int i = 0; i < len; i ++) {
            double alpha_uniform = this->alpha_uniform_distribution(generator);
            alpha_arr[i] = alpha_uniform;
        }
    } else {
        for(int i = 0; i < len; i ++) {
            alpha_arr[i] = this->alpha;
        }
    }
}

void ABM::PopulateWeightArrs(double* pa_weight_arr, double* rec_weight_arr, double* fit_weight_arr, int len) {
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    if(this->preferential_weight != -1 && this->recency_weight != -1 && this->fitness_weight != -1) {
        for(int i = 0; i < len; i ++) {
            double pa_uniform = this->preferential_weight;
            double rec_uniform = this->recency_weight;
            double fit_uniform = this->fitness_weight;
            double sum = pa_uniform + rec_uniform + fit_uniform;
            pa_weight_arr[i] = (double)pa_uniform / sum;
            rec_weight_arr[i] = (double)rec_uniform / sum;
            fit_weight_arr[i] = (double)fit_uniform / sum;
        }
    } else {
        for(int i = 0; i < len; i ++) {
            double pa_uniform = this->weights_uniform_distribution(generator);
            double rec_uniform = this->weights_uniform_distribution(generator);
            double fit_uniform = this->weights_uniform_distribution(generator);
            double sum = pa_uniform + rec_uniform + fit_uniform;
            pa_weight_arr[i] = (double)pa_uniform / sum;
            rec_weight_arr[i] = (double)rec_uniform / sum;
            fit_weight_arr[i] = (double)fit_uniform / sum;
        }
    }
}

void ABM::PopulateOutDegreeArr(int* out_degree_arr, int len) {
    std::uniform_int_distribution<int> outdegree_index_uniform_distribution{0, (int)(this->out_degree_bag_vec.size() - 1)};
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    for(int i = 0; i < len; i ++) {
        int index_uniform = outdegree_index_uniform_distribution(generator);
        out_degree_arr[i] = this->out_degree_bag_vec[index_uniform];
    }
}

void ABM::UpdateGraphAttributesWeights(Graph* graph, int next_node_id, double* pa_weight_arr, double* rec_weight_arr, double* fit_weight_arr, int len) {
    for(int i = 0; i < len; i ++) {
        int current_node_id = next_node_id + i;
        graph->SetDoubleAttribute("preferential_attachment_weight", current_node_id, pa_weight_arr[i]);
        graph->SetDoubleAttribute("recency_weight", current_node_id, rec_weight_arr[i]);
        graph->SetDoubleAttribute("fitness_weight", current_node_id, fit_weight_arr[i]);
    }
}

void ABM::UpdateGraphAttributesAlphas(Graph* graph, int next_node_id, double* alpha_arr, int len) {
    for(int i = 0; i < len; i ++) {
        int current_node_id = next_node_id + i;
        graph->SetDoubleAttribute("alpha", current_node_id, alpha_arr[i]);
    }
}

void ABM::UpdateGraphAttributesOutDegrees(Graph* graph, int next_node_id, int* out_degree_arr, int len) {
    for(int i = 0; i < len; i ++) {
        int current_node_id = next_node_id + i;
        graph->SetIntAttribute("assigned_out_degree", current_node_id, out_degree_arr[i]);
    }
}

std::vector<int> ABM::GetGraphAttributesGeneratorNodes(Graph* graph, int new_node) const {
    std::vector<int> generator_nodes;
    std::string generator_node_string = graph->GetStringAttribute("generator_node_string", new_node);
    std::stringstream ss(generator_node_string);
    std::string current_value;
    while(std::getline(ss, current_value, ';')) {
        generator_nodes.push_back(std::stoi(current_value));
    }
    return generator_nodes;
}

void ABM::UpdateGraphAttributesGeneratorNodes(Graph* graph, int new_node, const std::vector<int>& generator_nodes) {
    std::string generator_node_string;
    generator_node_string += std::to_string(generator_nodes.at(0));
    for(size_t i = 1; i < generator_nodes.size(); i ++) {
        generator_node_string += ";";
        generator_node_string += std::to_string(generator_nodes.at(i));
    }
    graph->SetStringAttribute("generator_node_string", new_node, generator_node_string);
}

void ABM::CalculateScores(std::unordered_map<int, double>& cached_results, int* src_arr, double* dst_arr, int len) {
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < len; i ++) {
        double current_dst = -1;
        if (src_arr[i] < 10000) {
            current_dst = cached_results[src_arr[i]];
        } else {
            current_dst = pow(src_arr[i], this->gamma) + 1;
            /* current_dst = (src_arr[i] * this->gamma) + 1; */
        }
        dst_arr[i] = current_dst;
        sum += current_dst;
    }
    #pragma omp parallel for
    for(int i = 0; i < len; i ++) {
        dst_arr[i] /= sum;
    }
}

void ABM::FillSameYearSourceNodes(std::set<int>& same_year_source_nodes, int current_year_new_nodes) {
    size_t num_same_year_source_nodes = (size_t)std::floor(current_year_new_nodes * this->same_year_proportion);
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    std::uniform_int_distribution<int> int_uniform_distribution(0, current_year_new_nodes - 1);
    while(same_year_source_nodes.size() != num_same_year_source_nodes) {
        int current_source = int_uniform_distribution(generator);
        if (same_year_source_nodes.count(current_source) == 0) {
            same_year_source_nodes.insert(current_source);
        }
    }
}

int ABM::MakeSameYearCitations(int num_new_nodes, const std::unordered_map<int, int>& reverse_continuous_node_mapping, int* citations, int current_graph_size) {
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    std::uniform_int_distribution<int> int_uniform_distribution(0, num_new_nodes - 1);
    int current_citation = int_uniform_distribution(generator);
    citations[0] = reverse_continuous_node_mapping.at(current_graph_size + current_citation);
    return 1;
}

int ABM::MakeUniformRandomCitations(Graph* graph, const std::unordered_map<int, int>& reverse_continuous_node_mapping, std::vector<int>& generator_nodes, int* citations, int num_cited_so_far, int num_citations) {
    if (num_citations <= 0) {
        return 0;
    }
    int actual_num_cited = num_citations;
    if ((int)graph->GetNodeSet().size() - num_cited_so_far - (int)generator_nodes.size() < num_citations) {
        actual_num_cited = (int)graph->GetNodeSet().size() - num_cited_so_far - (int)generator_nodes.size();
    }
    if (actual_num_cited == 0) {
        return actual_num_cited;
    }
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    std::uniform_int_distribution<int> int_uniform_distribution(0, (int)(graph->GetNodeSet().size() - 1));
    std::set<int> selected;
    for(int i = 0; i < num_cited_so_far; i ++) {
        selected.insert(citations[i]);
    }
    for(size_t i = 0; i < generator_nodes.size(); i ++) {
        selected.insert(generator_nodes.at(i));
    }
    int current_citation_index = 0;
    while(selected.size() != num_cited_so_far + generator_nodes.size() + actual_num_cited) {
        int current_citation = int_uniform_distribution(generator);
        if (!selected.contains(reverse_continuous_node_mapping.at(current_citation))) {
            citations[num_cited_so_far + current_citation_index] = reverse_continuous_node_mapping.at(current_citation);
            selected.insert(reverse_continuous_node_mapping.at(current_citation));
            current_citation_index ++;
        }
    }
    return actual_num_cited;
}

int ABM::MakeCitations(Graph* graph, const std::unordered_map<int, int>& continuous_node_mapping, int current_year, const std::vector<int>& candidate_nodes, int* citations, double* pa_arr, double* recency_arr, double* fit_arr, double pa_weight, double rec_weight, double fit_weight, int current_graph_size, int num_citations) {
    if (num_citations <= 0) {
        return 0;
    }
    if (candidate_nodes.size() <= 0) {
        return 0;
    }

    int actual_num_cited = num_citations;
    if (candidate_nodes.size() < (size_t)num_citations) {
        actual_num_cited = candidate_nodes.size();
    }
    std::vector<std::pair<double, int>> element_index_vec;
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    std::uniform_real_distribution<double> wrs_uniform_distribution{std::numeric_limits<double>::min(), 1};

    Eigen::MatrixXd current_scores(candidate_nodes.size(), 3);
    Eigen::Vector3d current_weights(pa_weight, rec_weight, fit_weight);
    for(size_t i = 0; i < candidate_nodes.size(); i ++) {
        int continuous_node_id = continuous_node_mapping.at(candidate_nodes.at(i));
        double current_pa = pa_arr[continuous_node_id];
        double current_rec = recency_arr[continuous_node_id];
        double current_fit = fit_arr[continuous_node_id];
        current_scores(i, 0) = current_pa;
        current_scores(i, 1) = current_rec;
        current_scores(i, 2) = current_fit;
    }
    Eigen::MatrixXd current_weighted_scores = current_scores * current_weights;
    auto current_wrs_uniform = [&] () {return wrs_uniform_distribution(generator);};
    Eigen::ArrayXd current_bases = Eigen::ArrayXd::NullaryExpr(candidate_nodes.size(), current_wrs_uniform);
    Eigen::ArrayXd weighted_random_sampling_results = current_bases.pow(1 / current_weighted_scores.array());
    for(size_t i = 0; i < candidate_nodes.size(); i ++) {
        element_index_vec.push_back({weighted_random_sampling_results(i), candidate_nodes.at(i)});
    }

    std::partial_sort(element_index_vec.begin(), element_index_vec.begin() + actual_num_cited, element_index_vec.end(), [](auto& left, auto& right){
        return left.first > right.first; // read
    });
    for (int i = 0; i < actual_num_cited; i ++) {
        citations[i] = element_index_vec[i].second;
    }
    return actual_num_cited;
}

std::vector<int> ABM::GetGeneratorNodes(Graph* graph, const std::unordered_map<int, int>& reverse_continuous_node_mapping) {
    std::vector<int> generator_nodes;
    std::uniform_int_distribution<int> generator_uniform_distribution{0, (int)(graph->GetNodeSet().size() - 1)};
    int num_generator_nodes = 1;
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    for(int i = 0; i < num_generator_nodes; i ++) {
        int continuous_generator_node = generator_uniform_distribution(generator);
        int generator_node = reverse_continuous_node_mapping.at(continuous_generator_node);
        generator_nodes.push_back(generator_node);
    }
    return generator_nodes;
}

std::unordered_map<int, std::vector<int>> ABM::GetOneAndTwoHopNeighborhood(Graph* graph, const std::vector<int>& generator_nodes, const std::unordered_map<int, int>& reverse_continuous_node_mapping) {
    std::unordered_map<int, std::vector<int>> one_and_two_hop_neighborhood_map;
    one_and_two_hop_neighborhood_map[1] = std::vector<int>();
    one_and_two_hop_neighborhood_map[2] = std::vector<int>();
    one_and_two_hop_neighborhood_map[1].reserve(this->neighborhood_sample);
    one_and_two_hop_neighborhood_map[2].reserve(this->neighborhood_sample);
    if (this->neighborhood_sample == -1) {
        std::set<int> visited;
        int num_hops = 2;
        for(size_t i = 0; i < generator_nodes.size(); i ++) {
            int generator_node = generator_nodes.at(i);
            std::queue<std::pair<int, int>> to_visit;
            to_visit.push({generator_node, 0});
            visited.insert(generator_node);
            while(!to_visit.empty()) {
                std::pair<int, int> current_pair = to_visit.front();
                to_visit.pop();
                int current_node = current_pair.first;
                int current_distance = current_pair.second;
                one_and_two_hop_neighborhood_map[current_distance].push_back(current_node);
                if (current_distance < num_hops) {
                    if(graph->GetOutDegree(current_node) > 0) {
                        for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(current_node)) {
                            if(!visited.contains(outgoing_neighbor)) {
                                visited.insert(outgoing_neighbor);
                                to_visit.push({outgoing_neighbor, current_distance + 1});
                            }
                        }
                    }
                    if(graph->GetInDegree(current_node) > 0) {
                        for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(current_node)) {
                            if(!visited.contains(incoming_neighbor)) {
                                visited.insert(incoming_neighbor);
                                to_visit.push({incoming_neighbor, current_distance + 1});
                            }
                        }
                    }
                }
            }
        }
    } else {
        size_t max_neighborhood_size = this->neighborhood_sample;
        std::set<int> visited;
        pcg_extras::seed_seq_from<std::random_device> rand_dev;
        pcg32 generator(rand_dev);
        for(size_t i = 0; i < generator_nodes.size(); i ++) {
            // get the 1-hop first
            int generator_node = generator_nodes.at(i);
            visited.insert(generator_node);
            std::vector<int> current_one_hop_neighborhood;
            if(graph->GetOutDegree(generator_node) > 0) {
                for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(generator_node)) {
                    if (!visited.contains(outgoing_neighbor)) {
                        current_one_hop_neighborhood.push_back(outgoing_neighbor);
                        visited.insert(outgoing_neighbor);
                    }
                }
            }
            if(graph->GetInDegree(generator_node) > 0) {
                for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(generator_node)) {
                    if (!visited.contains(incoming_neighbor)) {
                        current_one_hop_neighborhood.push_back(incoming_neighbor);
                        visited.insert(incoming_neighbor);
                    }
                }
            }
            // pick random nodes to get to 2-hop
            if (current_one_hop_neighborhood.size() > max_neighborhood_size) {
                std::vector<int> sampled_one_hop_neighborhood;
                std::sample(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), std::back_inserter(sampled_one_hop_neighborhood), max_neighborhood_size, generator);
                current_one_hop_neighborhood = sampled_one_hop_neighborhood;
            }
            // until here should be fast
            std::copy(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), std::back_inserter(one_and_two_hop_neighborhood_map[1]));
            std::shuffle(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), generator);
            for(size_t j = 0; j < current_one_hop_neighborhood.size(); j ++) {
                /* int current_two_hop_size = graph->GetForwardAdjMap().at(current_one_hop_neighborhood[j]).size(); */
                /* current_two_hop_size += graph->GetBackwardAdjMap().at(current_one_hop_neighborhood[j]).size(); */
                // worst case for one node we might only grab the outgoing edges
                std::vector<int> current_two_hop_neighborhood;
                if (graph->GetOutDegree(current_one_hop_neighborhood[j]) > 0) {
                    for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(current_one_hop_neighborhood[j])) {
                        if (!visited.contains(outgoing_neighbor)) {
                            visited.insert(outgoing_neighbor);
                            current_two_hop_neighborhood.push_back(outgoing_neighbor);
                        }
                    }
                }
                if (graph->GetInDegree(current_one_hop_neighborhood[j]) > 0) {
                    for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(current_one_hop_neighborhood[j])) {
                        if (!visited.contains(incoming_neighbor)) {
                            visited.insert(incoming_neighbor);
                            current_two_hop_neighborhood.push_back(incoming_neighbor);
                        }
                    }
                }
                /* if (one_and_two_hop_neighborhood_map[2].size() < max_neighborhood_size) { */

                if (one_and_two_hop_neighborhood_map[2].size() + current_two_hop_neighborhood.size() < max_neighborhood_size) {
                    std::copy(current_two_hop_neighborhood.begin(), current_two_hop_neighborhood.end(), std::back_inserter(one_and_two_hop_neighborhood_map[2]));
                } else {
                    size_t num_to_insert = max_neighborhood_size - one_and_two_hop_neighborhood_map[2].size();
                    std::sample(current_two_hop_neighborhood.begin(), current_two_hop_neighborhood.end(), std::back_inserter(one_and_two_hop_neighborhood_map[2]), num_to_insert, generator);
                    return one_and_two_hop_neighborhood_map;
                }
                /* } else { */
                /*     size_t num_to_insert = max_neighborhood_size - one_and_two_hop_neighborhood_map[2].size(); */
                /*     std::set<int> currently_sampled; */
                /*     std::uniform_int_distribution<int> generator_uniform_distribution{0, (int)(current_two_hop_neighborhood.size() - 1)}; */
                /*     while(currently_sampled.size() != num_to_insert) { */
                /*         int current_random_index = generator_uniform_distribution(generator); */
                /*         currently_sampled.insert(current_two_hop_neighborhood[current_random_index]); */
                /*     } */
                /*     std::copy(currently_sampled.begin(), currently_sampled.end(), std::back_inserter(one_and_two_hop_neighborhood_map[2])); */
                /*     return one_and_two_hop_neighborhood_map; */
                /* } */
                /* if (one_and_two_hop_neighborhood_map[2].size() < max_neighborhood_size) { */
                /*     return one_and_two_hop_neighborhood_map; */
                /* } */
                /* } else { */
                /*     size_t num_to_insert = max_neighborhood_size - one_and_two_hop_neighborhood_map[2].size(); */
                /*     /1* std::shuffle(current_two_hop_neighborhood.begin(), current_two_hop_neighborhood.end(), generator); *1/ */
                /*     /1* std::copy(current_two_hop_neighborhood.begin(), current_two_hop_neighborhood.begin() + num_to_insert, std::back_inserter(one_and_two_hop_neighborhood_map[2])); *1/ */
                /*     std::sample(current_two_hop_neighborhood.begin(), current_two_hop_neighborhood.end(), std::back_inserter(one_and_two_hop_neighborhood_map[2]), num_to_insert, generator); */
                /*     return one_and_two_hop_neighborhood_map; */
                /* } */
            }
        }
    }

    return one_and_two_hop_neighborhood_map;
}

void ABM::LogTime(int current_year, std::string label, int time_elapsed) {
    this->timing_file_handle << std::to_string(current_year) << "," << label << "," << std::to_string(time_elapsed) << "\n";
    std::flush(this->timing_file_handle);
}

void ABM::LogTime(int current_year, std::string label) {
    std::chrono::time_point<std::chrono::steady_clock> current_time = std::chrono::steady_clock::now();
    auto milliseconds_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - this->prev_time);
    this->timing_map[current_year].push_back({label, milliseconds_elapsed.count()});
    this->prev_time = current_time;
    this->timing_file_handle << std::to_string(current_year) << "," << label << "," << std::to_string(milliseconds_elapsed.count()) << "\n";
    std::flush(this->timing_file_handle);
}


std::chrono::time_point<std::chrono::steady_clock> ABM::LocalLogTime(std::vector<std::pair<std::string, int>>& local_parallel_stage_time_vec, std::chrono::time_point<std::chrono::steady_clock> local_prev_time, std::string label) {
    std::chrono::time_point<std::chrono::steady_clock> local_current_time = std::chrono::steady_clock::now();
    auto milliseconds_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(local_current_time - local_prev_time);
    local_parallel_stage_time_vec.push_back({label, milliseconds_elapsed.count()});
    return local_current_time;
}

void ABM::WriteTimingFile(int start_year, int end_year) {
    for(int i = start_year; i < end_year; i ++) {
        for(size_t j = 0; j < this->timing_map[i].size(); j ++) {
            this->timing_file_handle << std::to_string(i) << "," << (this->timing_map[i][j]).first << "," << std::to_string((this->timing_map[i][j]).second) << "\n";
        }
    }
}

std::unordered_map<int, std::vector<int>> ABM::GetOneAndTwoHopNeighborhoodFromMatrix(Graph* graph, int current_graph_size, const Eigen::SparseMatrix<int, 0, int>& two_hop_matrix, std::vector<int> generator_nodes, const std::unordered_map<int, int>& continuous_node_mapping, const std::unordered_map<int, int>& reverse_continuous_node_mapping) {
    std::unordered_map<int, std::vector<int>> one_and_two_hop_neighborhood_map;
    for(size_t i = 0; i < generator_nodes.size(); i ++) {
        // find 1-hop first
        if (graph->GetOutDegree(generator_nodes.at(i)) > 0) {
            for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(generator_nodes.at(i))) {
                one_and_two_hop_neighborhood_map[1].push_back(outgoing_neighbor);
            }
        }
        if (graph->GetInDegree(generator_nodes.at(i)) > 0) {
            for(auto const& incoming_neighbor : graph->GetForwardAdjMap().at(generator_nodes.at(i))) {
                one_and_two_hop_neighborhood_map[1].push_back(incoming_neighbor);
            }
        }
        // find 2-hop next
        int row_id = continuous_node_mapping.at(generator_nodes.at(i));
        /* for(int col = 0; col < current_graph_size; col ++) { */
        /*     if (two_hop_matrix(row_id, col) > 0) { */
        /*         int col_node_id = reverse_continuous_node_mapping.at(col); */
        /*         one_and_two_hop_neighborhood_map[2].push_back(col_node_id); */
        /*     } */
        /* } */
        for(Eigen::SparseMatrix<int, 0, int>::InnerIterator it(two_hop_matrix, row_id); it; ++it) {
            one_and_two_hop_neighborhood_map[2].push_back(it.value());
        }
    }
    return one_and_two_hop_neighborhood_map;
}

int ABM::main() {
    Graph* graph = new Graph(this->edgelist, this->nodelist);
    this->WriteToLogFile("loaded graph", Log::info);
    /* this->InitializeSeedFitness(graph); */
    this->InitializeFitness(graph);

    /* node ids to continous integer from 0 */
    std::unordered_map<int, int> continuous_node_mapping = this->BuildContinuousNodeMapping(graph);

    /* continous integer from 0 to node ids*/
    std::unordered_map<int, int> reverse_continuous_node_mapping = this->ReverseMapping(continuous_node_mapping);

    int start_year = this->GetMaxYear(graph) + 1;
    int next_node_id = this->GetMaxNode(graph) + 1;
    int initial_next_node_id = next_node_id;

    /* get input to score arrays based on continuous_node_mapping */
    int initial_graph_size = graph->GetNodeSet().size();
    int final_graph_size = this->GetFinalGraphSize(graph);
    this->WriteToLogFile("final graph size is " + std::to_string(final_graph_size), Log::info);
    int* in_degree_arr = new int[final_graph_size];
    int* fitness_arr = new int[final_graph_size];
    double* pa_arr = new double[final_graph_size];
    double* fit_arr = new double[final_graph_size];
    double* recency_arr = new double[final_graph_size];
    double* random_weight_arr = new double[final_graph_size];
    double* current_score_arr = new double[final_graph_size];

    // the first new agent node has index 0 but is actually index initial_graph_size in the continuous mapping
    double* pa_weight_arr = new double[final_graph_size - initial_graph_size];
    double* rec_weight_arr = new double[final_graph_size - initial_graph_size];
    double* fit_weight_arr = new double[final_graph_size - initial_graph_size];
    double* alpha_arr = new double[final_graph_size - initial_graph_size];
    int* out_degree_arr = new int[final_graph_size - initial_graph_size];

    this->PopulateWeightArrs(pa_weight_arr, rec_weight_arr, fit_weight_arr, final_graph_size - initial_graph_size);
    this->PopulateAlphaArr(alpha_arr, final_graph_size - initial_graph_size);
    this->PopulateOutDegreeArr(out_degree_arr, final_graph_size - initial_graph_size);

    std::vector<int> new_nodes_vec;
    std::set<int> same_year_source_nodes;
    std::vector<std::pair<int, int>> new_edges_vec;
    std::unordered_map<int, double> cached_results;
    for(int i = 0; i < 10000; i ++) {
        cached_results[i] = pow(i, this->gamma) + 1;
        /* cached_results[i] = (i * this->gamma) + 1; */
    }
    Eigen::setNbThreads(1);
    for (int current_year = start_year; current_year < start_year + this->num_cycles; current_year ++) {
        this->prev_time = std::chrono::steady_clock::now();
        int current_graph_size = graph->GetNodeSet().size();
        this->WriteToLogFile("current year is: " + std::to_string(current_year) + " and the graph is " + std::to_string(current_graph_size) + " nodes large", Log::info);
        this->FillInDegreeArr(graph, continuous_node_mapping, in_degree_arr);
        this->LogTime(current_year, "Fill in-degree array");
        this->FillFitnessArr(graph, continuous_node_mapping, current_year, fitness_arr);
        this->LogTime(current_year, "Fill fitness array");
        this->FillRecencyArr(graph, reverse_continuous_node_mapping, current_year, recency_arr);
        this->LogTime(current_year, "Fill recency array");
        this->CalculateScores(cached_results, in_degree_arr, pa_arr, current_graph_size);
        this->LogTime(current_year, "Process in-degree array");
        this->CalculateScores(cached_results, fitness_arr, fit_arr, current_graph_size);
        this->LogTime(current_year, "Process fitness array");

        /* initialize new nodes */
        int num_new_nodes = std::ceil(current_graph_size * this->growth_rate);
        this->WriteToLogFile("making " + std::to_string(num_new_nodes) + " nodes this year", Log::info);
        for(int i = 0; i < num_new_nodes; i ++) {
            continuous_node_mapping[next_node_id] = current_graph_size + i;
            reverse_continuous_node_mapping[current_graph_size + i] = next_node_id;
            new_nodes_vec.push_back(next_node_id);
            graph->SetIntAttribute("year", next_node_id, current_year);
            graph->SetStringAttribute("type", next_node_id, "agent");
            next_node_id ++;
        }
        this->LogTime(current_year, "Create new node ids");
        this->FillSameYearSourceNodes(same_year_source_nodes, new_nodes_vec.size());
        this->LogTime(current_year, "Pick same year nodes");

        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            int new_node = new_nodes_vec[i];
            std::vector<int> generator_nodes = this->GetGeneratorNodes(graph, reverse_continuous_node_mapping);
            this->UpdateGraphAttributesGeneratorNodes(graph, new_node, generator_nodes);
        }
        this->LogTime(current_year, "Pick generator nodes");
        std::vector<std::pair<int, int>> sampled_neighborhood_sizes_map(new_nodes_vec.size());

        std::vector<std::pair<std::string, int>> parallel_stage_time_vec;
        #pragma omp parallel for reduction(merge_int_pair_vecs: new_edges_vec) reduction(merge_str_int_pair_vecs: parallel_stage_time_vec)
        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            std::chrono::time_point<std::chrono::steady_clock> local_prev_time = std::chrono::steady_clock::now();
            std::vector<std::pair<int, int>> local_new_edges_vec;
            std::vector<std::pair<std::string, int>> local_parallel_stage_time_vec;


            int citations[250]; // out-degree assumed to be max 249
            int new_node = new_nodes_vec[i];
            int weight_arr_index = continuous_node_mapping[new_node] - initial_graph_size;
            double pa_weight = pa_weight_arr[weight_arr_index];
            double rec_weight = rec_weight_arr[weight_arr_index];
            double fit_weight = fit_weight_arr[weight_arr_index];
            double alpha = alpha_arr[weight_arr_index];
            std::vector<int> generator_nodes = this->GetGraphAttributesGeneratorNodes(graph, new_node);
            std::unordered_map<int, std::vector<int>> one_and_two_hop_neighborhood_map = this->GetOneAndTwoHopNeighborhood(graph, generator_nodes, reverse_continuous_node_mapping);
            sampled_neighborhood_sizes_map[i] = {one_and_two_hop_neighborhood_map[1].size(), one_and_two_hop_neighborhood_map[2].size()};
            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "retrieve one and two hop neighborhoods");
            int num_generator_node_citation = generator_nodes.size(); // should be 1 for now
            int same_year_citation = same_year_source_nodes.count(i); // could be 0 or 1
            int num_fully_random_cited_reserved = std::floor(this->fully_random_citations * out_degree_arr[weight_arr_index]); // e.g., 5% of out-degree. some small number

            // now the number of things to cite from distance 1 is the remaining citations * alpha. Call this remaining citations R for later.
            // unless distance 1 neighborhood is too small
            int num_citations_inside = std::ceil((out_degree_arr[weight_arr_index] - num_generator_node_citation - same_year_citation - num_fully_random_cited_reserved) * alpha);
            num_citations_inside = std::min(num_citations_inside, (int)one_and_two_hop_neighborhood_map[1].size());

            // now the number of things to cite from distance 2 is the remaining citations
            // unless distance 2 neighborhood is too small
            int num_citations_outside = out_degree_arr[weight_arr_index] - num_generator_node_citation - same_year_citation - num_fully_random_cited_reserved - num_citations_inside;
            num_citations_outside = std::min(num_citations_outside, (int)one_and_two_hop_neighborhood_map[2].size());

            // if it turns out that the 2-hop neighborhood (including 1 and 2) is small than R from earlier, then the leftover citations get cited randomly from the graph
            int num_fully_random_cited = out_degree_arr[weight_arr_index] - num_generator_node_citation - same_year_citation - num_citations_inside - num_citations_outside;

            int num_actually_cited = 0;
            if (same_year_citation) {
                num_actually_cited += this->MakeSameYearCitations(new_nodes_vec.size(), reverse_continuous_node_mapping, citations, current_graph_size);
            }

            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make same year citations");
            num_actually_cited += this->MakeCitations(graph, continuous_node_mapping, current_year, one_and_two_hop_neighborhood_map[1], citations + num_actually_cited, pa_arr, recency_arr, fit_arr, pa_weight, rec_weight, fit_weight, current_graph_size, num_citations_inside);
            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make 1-hop citations");
            num_actually_cited += this->MakeCitations(graph, continuous_node_mapping, current_year, one_and_two_hop_neighborhood_map[2], citations + num_actually_cited, pa_arr, recency_arr, fit_arr, pa_weight, rec_weight, fit_weight, current_graph_size, num_citations_outside);
            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make 2-hop citations");
            num_actually_cited += this->MakeUniformRandomCitations(graph, reverse_continuous_node_mapping, generator_nodes, citations, num_actually_cited, num_fully_random_cited);
            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make random citations");

            for(size_t j = 0; j < generator_nodes.size(); j ++) {
                local_new_edges_vec.push_back({new_node, generator_nodes[j]});
            }
            for(int j = 0; j < num_actually_cited; j ++) {
                local_new_edges_vec.push_back({new_node, citations[j]});
            }
            new_edges_vec.insert(new_edges_vec.end(), local_new_edges_vec.begin(), local_new_edges_vec.end());
            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "record edges");
            parallel_stage_time_vec.insert(parallel_stage_time_vec.end(), local_parallel_stage_time_vec.begin(), local_parallel_stage_time_vec.end());
        } // end of omp parallel loop
        std::map<std::string, int> per_stage_time_map;
        for(size_t i = 0; i < parallel_stage_time_vec.size(); i ++) {
            per_stage_time_map[parallel_stage_time_vec[i].first] +=  parallel_stage_time_vec[i].second;
        }
        this->LogTime(current_year, "find one and two hop neighborhoods", per_stage_time_map["retrieve one and two hop neighborhoods"]);
        this->LogTime(current_year, "make same year citations", per_stage_time_map["make same year citations"]);
        this->LogTime(current_year, "make 1-hop citations", per_stage_time_map["make 1-hop citations"]);
        this->LogTime(current_year, "make 2-hop citations", per_stage_time_map["make 2-hop citations"]);
        this->LogTime(current_year, "make random citations", per_stage_time_map["make random citations"]);
        this->LogTime(current_year, "record edges", per_stage_time_map["record edges"]);
        for(size_t i = 0; i < new_edges_vec.size(); i ++) {
            int new_node = new_edges_vec[i].first;
            int destination_id = new_edges_vec[i].second;
            graph->AddEdge({new_node, destination_id});
        }
        this->LogTime(current_year, "Add edges to graph");

        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            int new_node = new_nodes_vec[i];
            graph->SetIntAttribute("sampled_one_hop_neighborhood_size", new_node, sampled_neighborhood_sizes_map[i].first);
            graph->SetIntAttribute("sampled_two_hop_neighborhood_size", new_node, sampled_neighborhood_sizes_map[i].second);
        }

        this->LogTime(current_year, "Update graph attributes (neighborhood sizes)");
        this->AssignPeakFitnessValues(graph, new_nodes_vec);
        this->AssignFitnessLagDuration(graph, new_nodes_vec);
        this->AssignFitnessPeakDuration(graph, new_nodes_vec);
        this->PlantNodes(graph, new_nodes_vec, current_year - start_year + 1);
        this->LogTime(current_year, "Assign fitness values to new nodes");
        new_nodes_vec.clear();
        new_edges_vec.clear();
    }

    this->WriteToLogFile("finished sim", Log::info);
    graph->WriteGraph(this->output_file);

    this->UpdateGraphAttributesWeights(graph, initial_next_node_id, pa_weight_arr, rec_weight_arr, fit_weight_arr, final_graph_size - initial_graph_size);
    this->UpdateGraphAttributesAlphas(graph, initial_next_node_id, alpha_arr, final_graph_size - initial_graph_size);
    this->UpdateGraphAttributesOutDegrees(graph, initial_next_node_id, out_degree_arr, final_graph_size - initial_graph_size);

    for(auto const& node_id : graph->GetNodeSet()) {
        graph->SetIntAttribute("in_degree", node_id, graph->GetInDegree(node_id));
        graph->SetIntAttribute("out_degree", node_id, graph->GetOutDegree(node_id));
    }

    graph->WriteAttributes(this->auxiliary_information_file);
    this->WriteToLogFile("wrote graph", Log::info);
    delete[] in_degree_arr;
    delete[] fitness_arr;
    delete[] pa_arr;
    delete[] fit_arr;
    delete[] recency_arr;
    delete[] pa_weight_arr;
    delete[] rec_weight_arr;
    delete[] fit_weight_arr;
    delete[] alpha_arr;
    delete[] out_degree_arr;
    /* delete[] citations; */
    delete[] random_weight_arr;
    delete[] current_score_arr;

    delete graph;
    return 0;
}
