#include "abm.h"
#pragma omp declare reduction(merge_int_pair_vecs : std::vector<std::pair<int, int>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(merge_int_vecs : std::vector<int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

int ABM::WriteToLogFile(std::string message, Log message_type) {
    if(this->log_level >= message_type) {
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
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
    while(std::getline(out_degree_bag_stream, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        std::string index = current_line[0];
        if(index[0] == '#') {
            continue;
        }
        this->out_degree_bag_vec.push_back(std::stoi(current_line[1]));
        /* std::cout << "adding " << current_line[1] << " to  bag " << std::endl; */
    }
}


void ABM::ReadRecencyProbabilities() {
    char delimiter = ',';
    std::ifstream recency_probabilities_stream(this->recency_probabilities);
    std::string line;
    while(std::getline(recency_probabilities_stream, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        std::string year = current_line[0];
        if(year[0] == '#') {
            continue;
        }
        int integer_year_diff = std::stoi(current_line[0]);
        double probability = std::stod(current_line[1]);
        this->recency_probabilities_map[integer_year_diff] = probability;
    }
}

std::map<int, int> ABM::BuildContinuousNodeMapping(Graph* graph) {
    int next_node_id = 0;
    std::map<int, int> continuous_node_mapping;
    for(auto const& node : graph->GetNodeSet()) {
        continuous_node_mapping[node] = next_node_id;
        next_node_id ++;
    }
    return continuous_node_mapping;
}

std::map<int, int> ABM::ReverseMapping(std::map<int, int> mapping) {
    std::map<int, int> reverse_mapping;
    for(auto const& [key,val] : mapping) {
        reverse_mapping[val] = key;
    }
    return reverse_mapping;
}


void ABM::FillInDegreeArr(Graph* graph, const std::map<int, int>& continuous_node_mapping, int* in_degree_arr) {
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

void ABM::FillFitnessArr(Graph* graph, const std::map<int, int>& continuous_node_mapping, int current_year, int* fitness_arr) {
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


void ABM::FillRecencyArr(Graph* graph, const std::map<int, int>& reverse_continuous_node_mapping, int current_year, double* recency_arr) {
    std::map<int, int> year_count;
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
    /* #pragma omp parallel for simd */
    for(size_t i = 0; i < graph->GetNodeSet().size(); i ++) {
        int node_id = reverse_continuous_node_mapping.at(i);
        int current_published_year = graph->GetIntAttribute("year", node_id);
        int year_diff = current_year - current_published_year;
        recency_arr[i] = (float)this->recency_probabilities_map[year_diff] / year_count[year_diff];
    }
    // Mark: removed for node-level
    /* #pragma omp parallel for simd */
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
    std::random_device rand_dev;
    std::minstd_rand generator{rand_dev()};
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
    std::random_device rand_dev;
    std::minstd_rand generator{rand_dev()};
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
    std::uniform_int_distribution<int> outdegree_index_uniform_distribution{0, this->out_degree_bag_vec.size() - 1};
    std::random_device rand_dev;
    std::minstd_rand generator{rand_dev()};
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

void ABM::UpdateGraphAttributesGeneratorNodes(Graph* graph, int new_node, const std::vector<int>& generator_nodes) {
    std::string generator_node_string;
    generator_node_string += std::to_string(generator_nodes.at(0));
    for(size_t i = 1; i < generator_nodes.size(); i ++) {
        generator_node_string += ";";
        generator_node_string += std::to_string(generator_nodes.at(i));
    }
    graph->SetStringAttribute("generator_node_string", new_node, generator_node_string);
}

void ABM::CalculateScores(int* src_arr, double* dst_arr, int len) {
    double sum = 0;
    // Mark: removed for node-level
    /* #pragma omp parallel for simd */
    for(int i = 0; i < len; i ++) {
        dst_arr[i] = pow(src_arr[i], this->gamma) + 1;
    }
    // Mark: removed for node-level
    /* #pragma omp parallel for reduction(+:sum) */
    for(int i = 0; i < len; i ++) {
        sum += dst_arr[i];
    }
    // Mark: removed for node-level
    /* #pragma omp parallel for simd */
    for(int i = 0; i < len; i ++) {
        dst_arr[i] /= sum;
    }
}

void ABM::FillSameYearSourceNodes(std::set<int>& same_year_source_nodes, int current_year_new_nodes) {
    int num_same_year_source_nodes = int(current_year_new_nodes * this->same_year_proportion);
    std::random_device rand_dev;
    std::minstd_rand generator{rand_dev()};
    std::uniform_int_distribution<int> int_uniform_distribution(0, current_year_new_nodes - 1);
    while(same_year_source_nodes.size() != num_same_year_source_nodes) {
        int current_source = int_uniform_distribution(generator);
        if (same_year_source_nodes.count(current_source) == 0) {
            same_year_source_nodes.insert(current_source);
        }
    }
}

int ABM::MakeSameYearCitations(int num_new_nodes, const std::map<int, int>& reverse_continuous_node_mapping, int* citations, int current_graph_size) {
    std::random_device rand_dev;
    std::minstd_rand generator{rand_dev()};
    std::uniform_int_distribution<int> int_uniform_distribution(0, num_new_nodes - 1);
    int current_citation = int_uniform_distribution(generator);
    citations[0] = reverse_continuous_node_mapping.at(current_graph_size + current_citation);
    return 1;
}

int ABM::MakeUniformRandomCitations(Graph* graph, const std::map<int, int>& reverse_continuous_node_mapping, std::vector<int>& generator_nodes, int* citations, int num_cited_so_far, int num_citations) {
    if (num_citations == 0) {
        return 0;
    }
    int actual_num_cited = num_citations;
    this->WriteToLogFile("trying to uniformly cite " + std::to_string(actual_num_cited), Log::debug);
    if (graph->GetNodeSet().size() - num_cited_so_far - generator_nodes.size() < num_citations) {
        actual_num_cited = graph->GetNodeSet().size();
    }
    /* this->WriteToLogFile("can only cite " + std::to_string(actual_num_cited), Log::debug); */
    std::random_device rand_dev;
    std::minstd_rand generator{rand_dev()};
    std::uniform_int_distribution<int> int_uniform_distribution(0, graph->GetNodeSet().size() - 1);
    std::set<int> selected;
    int current_citation_index = 0;
    /* this->WriteToLogFile("currently chose " + std::to_string(selected.size()) + " things", Log::debug); */
    for(int i = 0; i < num_cited_so_far; i ++) {
        selected.insert(citations[i]);
    }
    for(int i = 0; i < generator_nodes.size(); i ++) {
        selected.insert(generator_nodes.at(i));
    }
    while(selected.size() != num_cited_so_far + actual_num_cited + generator_nodes.size()) {
        int current_citation = int_uniform_distribution(generator);
        /* this->WriteToLogFile("picking a node", Log::debug); */
        if (selected.count(reverse_continuous_node_mapping.at(current_citation)) == 0) {
            citations[num_cited_so_far + current_citation_index] = reverse_continuous_node_mapping.at(current_citation);
            selected.insert(current_citation);
            current_citation_index ++;
            /* this->WriteToLogFile("currently chose " + std::to_string(selected.size()) + " things", Log::debug); */
        }
    }
    /* this->WriteToLogFile("finished citing " + std::to_string(selected.size()), Log::debug); */
    /* std::cout << "cited " << actual_num_cited << " things from outside" << std::endl; */
    return actual_num_cited;
}

int ABM::MakeCitations(Graph* graph, const std::map<int, int>& continuous_node_mapping, int current_year, const std::vector<int>& candidate_nodes, int* citations, double* pa_arr, double* recency_arr, double* fit_arr, double pa_weight, double rec_weight, double fit_weight, int current_graph_size, int num_citations) {
    if (num_citations <= 0) {
        return 0;
    }
    if (candidate_nodes.size() <= 0) {
        return 0;
    }
    double* current_scores = new double[candidate_nodes.size() + 1];
    current_scores[candidate_nodes.size()] = 0.0;
    int local_continuous_node_id = 0;
    // Mark: removed for node-level
    /* #pragma omp parallel for simd */
    for(size_t i = 0; i < candidate_nodes.size(); i ++) {
        int continuous_node_id = continuous_node_mapping.at(candidate_nodes.at(i));
        double current_pa = pa_arr[continuous_node_id];
        double current_rec = recency_arr[continuous_node_id];
        double current_fit = fit_arr[continuous_node_id];
        current_scores[i] = (current_pa * pa_weight) + (current_rec * rec_weight) + (current_fit * fit_weight);
    }

    int actual_num_cited = num_citations;
    if (candidate_nodes.size() < num_citations) {
        actual_num_cited = candidate_nodes.size();
    }

    /* double sum = 0.0; */
    /* #pragma omp parallel for reduction(+:sum) */
    /* for(int i = 0; i < candidate_nodes.size(); i ++) { */
    /*     sum += current_scores[i]; */
    /* } */
    /* #pragma omp parallel for */
    /* for(int i = 0; i < candidate_nodes.size(); i ++) { */
    /*     current_scores[i] /= sum; */
    /* } */
    // begin
    /* for(int i = 0; i < actual_num_cited; i ++) { */
    /*     std::discrete_distribution<int> int_discrete_distribution(current_scores, current_scores + candidate_nodes.size()); */
    /*     int current_citation = int_discrete_distribution(generator); */
    /*     citations[i] = candidate_nodes.at(current_citation); */
    /*     current_scores[current_citation] = 0.0; */
    /* } */
    // BEGIN code that scales with size of array
    /* int highest_score_node_index = 0; */
    /* for(int i = 0; i < candidate_nodes.size(); i ++) { */
    /*     if(current_scores[i] > current_scores[highest_score_node_index]) { */
    /*         highest_score_node_index = i; */
    /*     } */
    /* } */
    /* std::cout << "highest score node id: " << std::to_string(candidate_nodes[highest_score_node_index]) << std::endl; */
    /* std::cout << "highest score out of 1: " << std::to_string(current_scores[highest_score_node_index]) << std::endl; */
    /* int highest_continuous_node_id = continuous_node_mapping.at(candidate_nodes.at(highest_score_node_index)); */
    /* std::cout << "highest score original pa: " << std::to_string(pa_arr[highest_continuous_node_id]) << std::endl; */
    /* std::cout << "highest score original rec: " << std::to_string(recency_arr[highest_continuous_node_id]) << std::endl; */
    /* std::cout << "highest score original fit: " << std::to_string(fit_arr[highest_continuous_node_id]) << std::endl; */
    /* std::cout << "highest weight original pa: " << std::to_string(pa_weight) << std::endl; */
    /* std::cout << "highest weight original rec: " << std::to_string(rec_weight) << std::endl; */
    /* std::cout << "highest weight original fit: " << std::to_string(fit_weight) << std::endl; */

    /* std::random_device rand_dev; */
    /* std::minstd_rand generator{rand_dev()}; */
    double* random_weight_arr = new double[candidate_nodes.size()];
    /* random_weight_arr[candidate_nodes.size()] = 0.0; */
    /* std::uniform_real_distribution<double> wrs_uniform_distribution{std::numeric_limits<double>::min(), 1}; */
    /* #pragma omp parallel for */
    std::random_device rand_dev;
    std::vector<std::minstd_rand> generator_vec;
    // num processors should be the same as omp max num
    for (int i = 0; i < this->num_processors; i ++) {
        generator_vec.push_back(std::minstd_rand(rand_dev()));
    }

    // Mark: removed for node-level
    /* #pragma omp parallel for simd */
    /* std::random_device rand_dev; */
    std::minstd_rand generator{rand_dev()};
    for(int i = 0; i < candidate_nodes.size(); i ++) {
        std::uniform_real_distribution<double> wrs_uniform_distribution{std::numeric_limits<double>::min(), 1};
        /* double wrs_uniform = wrs_uniform_distribution(generator_vec.at(omp_get_thread_num())); */
        double wrs_uniform = wrs_uniform_distribution(generator);
        if (current_scores[i] != 0) {
            random_weight_arr[i] = pow(wrs_uniform, 1.0/current_scores[i]);
        } else {
            random_weight_arr[i] = 0;
        }
    }

    // TODO: probably omp parallel it
    std::vector<std::pair<double, int>> element_index_vec(candidate_nodes.size());
    // Mark: removed for node-level
    /* #pragma omp parallel for simd */
    for (size_t i = 0; i < candidate_nodes.size(); i ++) {
        /* element_index_vec.push_back({random_weight_arr[i], candidate_nodes.at(i)}); */
        element_index_vec[i] = {random_weight_arr[i], candidate_nodes.at(i)};
    }
    /* std::cout << "before sorted citing" << std::endl; */
    /* std::cout << "num cited should be " << actual_num_cited << std::endl; */
    /* for (size_t i = 0; i < candidate_nodes.size(); i ++) { */
    /*     std::cout << element_index_vec.at(i).first << ":" << element_index_vec.at(i).second << std::endl; */
    /* } */
    /* std::partial_sort(random_weight_arr, random_weight_arr + actual_num_cited, random_weight_arr + candidate_nodes.size(), std::greater<std::pair<int, int>>()); */
    std::partial_sort(element_index_vec.begin(), element_index_vec.begin() + actual_num_cited, element_index_vec.end(), [](auto& left, auto& right){
        return left.first > right.first;
    });

    /* std::cout << "after sorted citing" << std::endl; */
    /* for (size_t i = 0; i < candidate_nodes.size(); i ++) { */
    /*     std::cout << element_index_vec.at(i).first << ":" << element_index_vec.at(i).second << std::endl; */
    /* } */
    /* std::partial_sort(random_weight_arr, element_index_vec.end(), [](auto& left, auto& right) { */
    /*     return left.first < right.first; */
    /* }); */
    for (int i = 0; i < actual_num_cited; i ++) {
        citations[i] = element_index_vec[i].second;
    }
    /* return actual_num_cited; */

    // end
    delete[] current_scores;
    delete[] random_weight_arr;
    return actual_num_cited;
}

std::vector<int> ABM::GetComplement(Graph* graph, const std::vector<int>& base_vec, const std::map<int, int>& reverse_continuous_node_mapping) {
    std::vector<int> complement;
    std::uniform_int_distribution<int> generator_uniform_distribution{0, graph->GetNodeSet().size() - 1};
    std::set<int> base_set(base_vec.begin(), base_vec.end());
    // Mark: removed for node-level
    /* #pragma omp parallel for reduction(merge_int_vecs: complement) */
    for(size_t i = 0; i < graph->GetNodeSet().size(); i ++) {
        int node_id = reverse_continuous_node_mapping.at(i);
        if (!base_set.contains(node_id)) {
            complement.push_back(node_id);
        }
    }
    return complement;
}

std::vector<int> ABM::GetGeneratorNodes(Graph* graph, const std::map<int, int>& reverse_continuous_node_mapping) {
    std::vector<int> generator_nodes;
    std::uniform_int_distribution<int> generator_uniform_distribution{0, graph->GetNodeSet().size() - 1};
    int num_generator_nodes = 1;
    std::random_device rand_dev;
    std::minstd_rand generator{rand_dev()};
    for(int i = 0; i < num_generator_nodes; i ++) {
        int continuous_generator_node = generator_uniform_distribution(generator);
        int generator_node = reverse_continuous_node_mapping.at(continuous_generator_node);
        generator_nodes.push_back(generator_node);
    }
    return generator_nodes;
}

std::map<int, std::vector<int>> ABM::GetOneAndTwoHopNeighborhood(Graph* graph, const std::vector<int>& generator_nodes, const std::map<int, int>& reverse_continuous_node_mapping) {
    std::map<int, std::vector<int>> one_and_two_hop_neighborhood_map;
    std::set<int> visited;
    int num_hops = 2;
    for(int i = 0; i < generator_nodes.size(); i ++) {
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
    return one_and_two_hop_neighborhood_map;
}

std::vector<int> ABM::GetNeighborhood(Graph* graph, const std::vector<int>& generator_nodes, const std::map<int, int>& reverse_continuous_node_mapping) {
    std::vector<int> neighborhood;
    std::set<int> visited;
    std::uniform_int_distribution<int> generator_uniform_distribution{0, graph->GetNodeSet().size() - 1};
    int num_hops = 1;
    /* std::cout << "getting 2-hop neighborhood from: " << generator_nodes.at(0) << std::endl; */
    for(int i = 0; i < generator_nodes.size(); i ++) {
        int generator_node = generator_nodes.at(i);
        std::queue<std::pair<int, int>> to_visit;
        to_visit.push({generator_node, 0});
        visited.insert(generator_node);
        while(!to_visit.empty()) {
            std::pair<int, int> current_pair = to_visit.front();
            to_visit.pop();
            int current_node = current_pair.first;
            int current_distance = current_pair.second;
            if (current_distance < num_hops) {
                if(graph->GetOutDegree(current_node) > 0) {
                    for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(current_node)) {
                        if(!visited.contains(outgoing_neighbor)) {
                            neighborhood.push_back(outgoing_neighbor);
                            visited.insert(outgoing_neighbor);
                            to_visit.push({outgoing_neighbor, current_distance + 1});
                        }
                    }
                }
                if(graph->GetInDegree(current_node) > 0) {
                    for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(current_node)) {
                        if(!visited.contains(incoming_neighbor)) {
                            neighborhood.push_back(incoming_neighbor);
                            visited.insert(incoming_neighbor);
                            to_visit.push({incoming_neighbor, current_distance + 1});
                        }
                    }
                }
            }
        }
    }
    return neighborhood;
}

/*
void ABM::InitializeAuthors(Graph* graph, const std::map<int, int>& continuous_node_mapping, std::map<int, std::vector<int>>& author_to_publication_map, std::map<int, std::vector<int>>& number_published_to_author_map, std::map<int, int>& author_remaining_years_map) {
    for(int i = 1; i <= this->max_author_lifetime; i ++) {
        number_published_to_author_map[i] = std::set<int>();
    }
    std::uniform_int_distribution<int> generator_uniform_lifetime_distribution{1, this->max_author_lifetime};
    for(size_t i = 0; i < graph->GetNodeSet().size(); i ++) {
        int current_graph_size = i + 1;
        // int current_pub_count is that the current node should be authored by someone with this number of publications
        for(int current_pub_count = 1; current_pub_count <= this->max_author_lifetime; current_pub_count ++) {
            int minimum_num_authors_required = (float)(current_graph_size) / pow(current_pub_count, this->k);
            int actual_num_authors = number_published_to_author_map[current_pub_count].size();
            if (actual_num_authors < minimum_num_authors_required) {
                break;
            }
        }
        number_published_to_author_map[current_pub_count] += 1
        if(current_pub_count > 1 && number_published_to_author_map[current_pub_count - 1] > 0) {
            // an author wrote another article
            number_published_to_author_map[current_pub_count - 1] -= 1
            std::uniform_int_distribution<int> generator_uniform_distribution{0, number_published_to_author_map[current_pub_count - 1] - 1};
            int random_author_index = generator_uniform_distribution(this->generator);
            int random_author_id = number_published_to_author_map[current_pub_count - 1][random_author_index];
            author_to_publication_map[random_author_id] = continuous_node_mapping.at(i);
        } else {
            // should be a new author
            author_to_publication_map[this->next_author_id] = continuous_node_mapping.at(i);
            int random_year = generator_uniform_lifetime_distribution(this->generator);
            author_remaining_years_map[this->next_author_id] = random_year;
            this->next_author_id += 1
        }
    }
}

void ABM::AssignAuthor(Graph* graph, int new_node, std::map<int, std::vector<int>>& author_to_publication_map, std::map<int, std::vector<int>>& number_published_to_author_map, std::map<int, int>& author_remaining_years_map) {
    int current_graph_size = graph->GetNodeSet().size();
    // int current_pub_count is that the current node should be authored by someone with this number of publications
    std::uniform_int_distribution<int> generator_uniform_lifetime_distribution{1, this->max_author_lifetime};
    for(int current_pub_count = 1; current_pub_count <= this->max_author_lifetime; current_pub_count ++) {
        int minimum_num_authors_required = (float)(current_graph_size) / pow(current_pub_count, this->k);
        int actual_num_authors = number_published_to_author_map[current_pub_count].size();
        if (actual_num_authors < minimum_num_authors_required) {
            break;
        }
    }
    number_published_to_author_map[current_pub_count] += 1
    if(current_pub_count > 1 && number_published_to_author_map[current_pub_count - 1] > 0) {
        // an author wrote another article
        number_published_to_author_map[current_pub_count - 1] -= 1
        std::uniform_int_distribution<int> generator_uniform_distribution{0, number_published_to_author_map[current_pub_count - 1] - 1};
        int random_author_index = generator_uniform_distribution(this->generator);
        int random_author_id = number_published_to_author_map[current_pub_count - 1][random_author_index];
        author_to_publication_map[random_author_id] = new_node;
    } else {
        // should be a new author
        author_to_publication_map[this->next_author_id] = new_node;
        int random_year = generator_uniform_lifetime_distribution(this->generator);
        author_remaining_years_map[this->next_author_id] = random_year;
        this->next_author_id += 1
    }
}

void ABM::AgeAuthors(std::map<int, std::vector<int>>& author_to_publication_map, std::map<int, std::vector<int>>& number_published_to_author_map, std::map<int, int>& author_remaining_years_map) {
    for(const auto& [author_id, remaining_lifetime] : author_remaining_years_map) {
        author_remaining_years_map[author_id] -= 1;
        if (remaining_lifetime == 1) {
            // remove author
            //int num_pubs_by_author = author_to_publication_map[author_id].size();
            //number_published_to_author_map[num_pubs_by_author] -= 1;
            //author_to_publication_map.erase(author_id);
        }
    }
}
*/

int ABM::main() {
    /* reading input edgelist, nodelist, outdegree bag, recency probabilities */
    Graph* graph = new Graph(this->edgelist, this->nodelist);
    /* graph->PrintGraph(); */
    this->WriteToLogFile("loaded graph", Log::info);
    this->InitializeFitness(graph);
    this->WriteToLogFile("initialized fitness for the seed graph", Log::debug);

    /* node ids to continous integer from 0 */
    std::map<int, int> continuous_node_mapping = this->BuildContinuousNodeMapping(graph);

    this->WriteToLogFile("forward built", Log::debug);
    /* continous integer from 0 to node ids*/
    std::map<int, int> reverse_continuous_node_mapping = this->ReverseMapping(continuous_node_mapping);
    this->WriteToLogFile("reverse mapping built", Log::debug);

    int start_year = this->GetMaxYear(graph) + 1;
    int next_node_id = this->GetMaxNode(graph) + 1;
    int initial_next_node_id = next_node_id;

    /* get input to score arrays based on continuous_node_mapping */
    int initial_graph_size = graph->GetNodeSet().size();
    int final_graph_size = this->GetFinalGraphSize(graph);
    this->WriteToLogFile("final graph size is " + std::to_string(final_graph_size), Log::info);
    int* in_degree_arr = new int[final_graph_size];
    this->WriteToLogFile("alloc 1", Log::debug);
    int* fitness_arr = new int[final_graph_size];
    this->WriteToLogFile("alloc 2", Log::debug);
    double* pa_arr = new double[final_graph_size];
    this->WriteToLogFile("alloc 3", Log::debug);
    double* fit_arr = new double[final_graph_size];
    this->WriteToLogFile("alloc 4", Log::debug);
    double* recency_arr = new double[final_graph_size];
    this->WriteToLogFile("alloc 5", Log::debug);
    double* random_weight_arr = new double[final_graph_size];
    this->WriteToLogFile("alloc 6", Log::debug);
    double* current_score_arr = new double[final_graph_size];
    this->WriteToLogFile("alloc 7", Log::debug);
    int* citations = new int[250]; // maximum size container with 250 assumption for now
    this->WriteToLogFile("alloc 8", Log::debug);

    // the first new agent node has index 0 but is actually index initial_graph_size in the continuous mapping
    double* pa_weight_arr = new double[final_graph_size - initial_graph_size];
    this->WriteToLogFile("alloc 9", Log::debug);
    double* rec_weight_arr = new double[final_graph_size - initial_graph_size];
    this->WriteToLogFile("alloc 10", Log::debug);
    double* fit_weight_arr = new double[final_graph_size - initial_graph_size];
    this->WriteToLogFile("alloc 11", Log::debug);
    double* alpha_arr = new double[final_graph_size - initial_graph_size];
    this->WriteToLogFile("alloc 12", Log::debug);
    int* out_degree_arr = new int[final_graph_size - initial_graph_size];
    this->WriteToLogFile("alloc 13", Log::debug);
    this->WriteToLogFile("allocated arrays", Log::debug);

    this->PopulateWeightArrs(pa_weight_arr, rec_weight_arr, fit_weight_arr, final_graph_size - initial_graph_size);
    this->WriteToLogFile("populated weight arrays", Log::debug);
    this->PopulateAlphaArr(alpha_arr, final_graph_size - initial_graph_size);
    this->WriteToLogFile("populated alpha array", Log::debug);
    this->PopulateOutDegreeArr(out_degree_arr, final_graph_size - initial_graph_size);
    this->WriteToLogFile("populated out degree array", Log::debug);

    /* std::map<int, std::set<int>> author_to_publication_map; // map[1] = {98} meaning author 1 published the node 98 */
    /* std::map<int, std::set<int>> number_published_to_author_map; // map[1] = {1, 3, 4} meaning authors 1, 3, and 3 have 1 publications each */
    /* std::map<int, int> author_remaining_years_map; // map[1] = 2 meaning  author 1 has 2 more years left to live */
    /* this->InitializeAuthors(graph, continuous_node_mapping, author_to_publication_map, number_published_to_author_map, author_remaining_years_map); */
    /* graph->SetIntAttribute("fitness_peak_value", 4920089, 1000); */
    /* graph->SetIntAttribute("fitness_lag_duration", 4920089, 0); */
    /* graph->SetIntAttribute("fitness_peak_duration", 4920089, 1000); */

    std::vector<int> new_nodes_vec;
    std::set<int> same_year_source_nodes;
    std::vector<std::pair<int, int>> new_edges_vec;
    for (int current_year = start_year; current_year < start_year + this->num_cycles; current_year ++) {
        /* std::cout << "new year" << std::endl; */
        int current_graph_size = graph->GetNodeSet().size();
        this->WriteToLogFile("current year is: " + std::to_string(current_year) + " and the graph is " + std::to_string(current_graph_size) + " nodes large", Log::info);
        this->FillInDegreeArr(graph, continuous_node_mapping, in_degree_arr);
        this->WriteToLogFile("indegree for current year filled", Log::debug);
        this->FillFitnessArr(graph, continuous_node_mapping, current_year, fitness_arr);
        this->WriteToLogFile("fitness for current year filled", Log::debug);
        this->FillRecencyArr(graph, reverse_continuous_node_mapping, current_year, recency_arr);
        this->WriteToLogFile("recency for current year calculated", Log::debug);
        this->CalculateScores(in_degree_arr, pa_arr, current_graph_size);
        /* double current_max_pa_score = pa_arr[0]; */
        /* int current_max_pa_node_id = reverse_continuous_node_mapping.at(0); */
        /* for(int i = 0; i < current_graph_size; i ++) { */
        /*     if(pa_arr[i] > current_max_pa_score) { */
        /*         current_max_pa_score = pa_arr[i]; */
        /*         current_max_pa_node_id = reverse_continuous_node_mapping.at(i); */
        /*     } */
        /* } */
        /* this->WriteToLogFile("max indegree pa score: " + std::to_string(current_max_pa_score), Log::debug); */
        /* this->WriteToLogFile("max indegree pa node id: " + std::to_string(current_max_pa_node_id), Log::debug); */

        this->WriteToLogFile("indegree gammad", Log::debug);
        this->CalculateScores(fitness_arr, fit_arr, current_graph_size);
        this->WriteToLogFile("fitness gammad", Log::debug);

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
        this->WriteToLogFile("all new nodes initialized with years and mapped", Log::debug);
        this->FillSameYearSourceNodes(same_year_source_nodes, new_nodes_vec.size());


        #pragma omp parallel for reduction(merge_int_pair_vecs: new_edges_vec)
        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            int citations[250];
            /* int* citations = new int[250]; // maximum size container with 250 assumption for now */
            this->WriteToLogFile("starting node " + std::to_string(i) + "/" + std::to_string(new_nodes_vec.size()), Log::debug);
            int new_node = new_nodes_vec[i];
            /* this->AssignAuthor(graph, new_node, author_to_publication_map, number_published_to_author_map, author_remaining_years_map); */
            int weight_arr_index = continuous_node_mapping[new_node] - initial_graph_size;
            double pa_weight = pa_weight_arr[weight_arr_index];
            double rec_weight = rec_weight_arr[weight_arr_index];
            double fit_weight = fit_weight_arr[weight_arr_index];
            double alpha = alpha_arr[weight_arr_index];
            std::vector<int> generator_nodes = this->GetGeneratorNodes(graph, reverse_continuous_node_mapping);
            #pragma omp critical
            {
                this->UpdateGraphAttributesGeneratorNodes(graph, new_node, generator_nodes);
            }
            /* std::vector<int> neighborhood = this->GetNeighborhood(graph, generator_nodes, reverse_continuous_node_mapping); */
            std::map<int, std::vector<int>> one_and_two_hop_neighborhood_map = this->GetOneAndTwoHopNeighborhood(graph, generator_nodes, reverse_continuous_node_mapping);
            /* this->WriteToLogFile("neighborhood size is " + std::to_string(neighborhood.size()), Log::debug); */


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

            num_actually_cited += this->MakeCitations(graph, continuous_node_mapping, current_year, one_and_two_hop_neighborhood_map[1], citations + num_actually_cited, pa_arr, recency_arr, fit_arr, pa_weight, rec_weight, fit_weight, current_graph_size, num_citations_inside);
            num_actually_cited += this->MakeCitations(graph, continuous_node_mapping, current_year, one_and_two_hop_neighborhood_map[2], citations + num_actually_cited, pa_arr, recency_arr, fit_arr, pa_weight, rec_weight, fit_weight, current_graph_size, num_citations_outside);
            /* #pragma omp critical */
            /* { */
            /* this->WriteToLogFile("new node id:" + std::to_string(new_node), Log::debug); */
            /* this->WriteToLogFile("assigned out-degree:" + std::to_string(out_degree_arr[weight_arr_index]), Log::debug); */
            /* this->WriteToLogFile("total num cited not including generator: " + std::to_string(num_actually_cited), Log::debug); */
            /* this->WriteToLogFile("cited " + std::to_string(num_generator_node_citation) + " as generator", Log::debug); */
            /* this->WriteToLogFile("cited " + std::to_string(same_year_citation) + " as same year", Log::debug); */
            /* this->WriteToLogFile("initially reserved " + std::to_string(num_fully_random_cited_reserved) + " as random", Log::debug); */
            /* this->WriteToLogFile("cited " + std::to_string(num_citations_inside) + " as inside", Log::debug); */
            /* this->WriteToLogFile("cited " + std::to_string(num_citations_outside) + " as outside", Log::debug); */
            /* this->WriteToLogFile("cited " + std::to_string(num_fully_random_cited) + " as random", Log::debug); */
            /* this->WriteToLogFile("1-hop size: " + std::to_string(one_and_two_hop_neighborhood_map[1].size()), Log::debug); */
            /* this->WriteToLogFile("2-hop size: " + std::to_string(one_and_two_hop_neighborhood_map[2].size()), Log::debug); */
            /* } */
            num_actually_cited += this->MakeUniformRandomCitations(graph, reverse_continuous_node_mapping, generator_nodes, citations, num_actually_cited, num_fully_random_cited);

            for(size_t i = 0; i < generator_nodes.size(); i ++) {
                new_edges_vec.push_back({new_node, generator_nodes[i]});
            }
            for(int i = 0; i < num_actually_cited; i ++) {
                /* if(citations[i] == current_max_pa_node_id) { */
                /*     this->WriteToLogFile("max indegree pa node id cited", Log::debug); */
                /* } */
                /* if(citations[i] == 4920089) { */
                /*     std::cout << "4920089 cited" << std::endl; */
                /* } */
                new_edges_vec.push_back({new_node, citations[i]});
            }
            /* this->WriteToLogFile("edge vec now has " + std::to_string(new_edges_vec.size()) + " edges", Log::debug); */
        }

        /* this->WriteToLogFile("edges saved to vector", Log::debug); */
        for(size_t i = 0; i < new_edges_vec.size(); i ++) {
            int new_node = new_edges_vec[i].first;
            int destination_id = new_edges_vec[i].second;
            graph->AddEdge({new_node, destination_id});
        }
        this->WriteToLogFile("edges saved to graph", Log::debug);
        this->AssignPeakFitnessValues(graph, new_nodes_vec);
        this->WriteToLogFile("assigned peak fitness for new nodes", Log::debug);
        this->AssignFitnessLagDuration(graph, new_nodes_vec);
        this->WriteToLogFile("assigned fitness lag duration for new nodes", Log::debug);
        this->AssignFitnessPeakDuration(graph, new_nodes_vec);
        this->WriteToLogFile("assigned fitness peak duration for new nodes", Log::debug);
        this->PlantNodes(graph, new_nodes_vec, current_year - start_year + 1);
        /* this->WriteToLogFile("assigned fitness peak duration for new nodes", Log::debug); */
        /* this->AgeAuthors(author_to_publication_map, number_published_to_author_map, author_remaining_years_map); */
        new_nodes_vec.clear();
        new_edges_vec.clear();
        this->WriteToLogFile("writing temp graph", Log::info);
        graph->WriteGraph(this->output_file + "_" + std::to_string(current_year));
    }

    this->WriteToLogFile("finished sim", Log::info);
    graph->WriteGraph(this->output_file);

    this->UpdateGraphAttributesWeights(graph, initial_next_node_id, pa_weight_arr, rec_weight_arr, fit_weight_arr, final_graph_size - initial_graph_size);
    this->WriteToLogFile("update 1", Log::debug);
    this->UpdateGraphAttributesAlphas(graph, initial_next_node_id, alpha_arr, final_graph_size - initial_graph_size);
    this->WriteToLogFile("update 2", Log::debug);
    this->UpdateGraphAttributesOutDegrees(graph, initial_next_node_id, out_degree_arr, final_graph_size - initial_graph_size);
    this->WriteToLogFile("update 3", Log::debug);

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
