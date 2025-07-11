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
    /* std::random_device rand_dev; */
    /* std::minstd_rand generator{rand_dev()}; */
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
    /* std::random_device rand_dev; */
    /* std::minstd_rand generator{rand_dev()}; */
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
    /* std::random_device rand_dev; */
    /* std::minstd_rand generator{rand_dev()}; */
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
    /* std::random_device rand_dev; */
    /* std::minstd_rand generator{rand_dev()}; */
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
    /* std::random_device rand_dev; */
    /* std::minstd_rand generator{rand_dev()}; */
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
    this->WriteToLogFile("trying to uniformly cite " + std::to_string(actual_num_cited), Log::debug);
    /* std::cerr << "trying to uniformly cite " << std::to_string(actual_num_cited) << std::endl; */
    if ((int)graph->GetNodeSet().size() - num_cited_so_far - (int)generator_nodes.size() < num_citations) {
        actual_num_cited = (int)graph->GetNodeSet().size() - num_cited_so_far - (int)generator_nodes.size();
    }
    if (actual_num_cited == 0) {
        return actual_num_cited;
    }
    /* if (actual_num_cited < 0) { */
    /*     std::cerr << "trying to uniformly cite " << std::to_string(num_citations) << std::endl; */
    /*     std::cerr << "but it got capped to " << std::to_string(actual_num_cited) << std::endl; */
    /*     std::cerr << "graph size: " << std::to_string(graph->GetNodeSet().size()) << std::endl; */
    /*     std::cerr << "num cited so far in other places:  " << std::to_string(num_cited_so_far) << std::endl; */
    /*     std::cerr << "num generator nodes:  " << std::to_string(generator_nodes.size()) << std::endl; */
    /*     exit(1); */
    /* } else { */
    /*     std::cerr << "trying to uniformly cite " << std::to_string(num_citations) << std::endl; */
    /*     std::cerr << "but it got capped to " << std::to_string(actual_num_cited) << std::endl; */
    /*     std::cerr << "graph size: " << std::to_string(graph->GetNodeSet().size()) << std::endl; */
    /*     std::cerr << "num cited so far in other places:  " << std::to_string(num_cited_so_far) << std::endl; */
    /*     std::cerr << "num generator nodes:  " << std::to_string(generator_nodes.size()) << std::endl; */
    /* } */
    /* std::cerr << "but it got capped to " << std::to_string(actual_num_cited) << std::endl; */
    /* std::cerr << "graph size: " << std::to_string(graph->GetNodeSet().size()) << std::endl; */
    /* std::cerr << "num cited so far in other places:  " << std::to_string(num_cited_so_far) << std::endl; */
    /* std::cerr << "num generator nodes:  " << std::to_string(generator_nodes.size()) << std::endl; */
    /* this->WriteToLogFile("can only cite " + std::to_string(actual_num_cited), Log::debug); */
    /* std::random_device rand_dev; */
    /* std::minstd_rand generator{rand_dev()}; */
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    std::uniform_int_distribution<int> int_uniform_distribution(0, (int)(graph->GetNodeSet().size() - 1));
    std::set<int> selected;
    int current_citation_index = 0;
    /* this->WriteToLogFile("currently chose " + std::to_string(selected.size()) + " things", Log::debug); */
    for(int i = 0; i < num_cited_so_far; i ++) {
        selected.insert(citations[i]);
    }
    for(size_t i = 0; i < generator_nodes.size(); i ++) {
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
    /* std::cerr << "cited " << actual_num_cited << " things uniformly" << std::endl; */
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
    /* std::random_device rand_dev; */
    /* std::minstd_rand generator{rand_dev()}; */
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    std::uniform_real_distribution<double> wrs_uniform_distribution{std::numeric_limits<double>::min(), 1};

    Eigen::MatrixXd current_scores(candidate_nodes.size(), 3);
    Eigen::Vector3d current_weights(pa_weight, rec_weight, fit_weight);
    /* auto cmp = [](const std::pair<double, int> &left, const std::pair<double, int> &right) { */
    /*     return left.first > right.first; */
    /* }; */
    /* std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, decltype(cmp)> min_heap(cmp); */

    for(size_t i = 0; i < candidate_nodes.size(); i ++) {
        int continuous_node_id = continuous_node_mapping.at(candidate_nodes.at(i));
        double current_pa = pa_arr[continuous_node_id];
        double current_rec = recency_arr[continuous_node_id];
        double current_fit = fit_arr[continuous_node_id];
        current_scores(i, 0) = current_pa;
        current_scores(i, 1) = current_rec;
        current_scores(i, 2) = current_fit;
        /* double weighted_score = pa_weight * current_pa + rec_weight * current_rec + fit_weight * current_fit; */
        /* double base = wrs_uniform_distribution(generator); */
        /* double wrs_score = std::pow(base, 1.0 / weighted_score); */
        /* if ((int)min_heap.size() < actual_num_cited) { */
        /*     min_heap.emplace(wrs_score, candidate_nodes.at(i)); */
        /* } else if (wrs_score > min_heap.top().first) { */
        /*     min_heap.pop(); */
        /*     min_heap.emplace(wrs_score, candidate_nodes.at(i)); */
        /* } */
    }
    /* for (int i = 0; i < actual_num_cited; i++) { */
    /*     citations[i] = min_heap.top().second; */
    /*     min_heap.pop(); */
    /* } */
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
    /* std::nth_element(element_index_vec.begin(), element_index_vec.begin() + actual_num_cited, element_index_vec.end(), [](auto& left, auto& right) { */
    /*     return left.first > right.first; */
    /* }); */
    for (int i = 0; i < actual_num_cited; i ++) {
        citations[i] = element_index_vec[i].second;
    }
    return actual_num_cited;
}

std::vector<int> ABM::GetComplement(Graph* graph, const std::vector<int>& base_vec, const std::unordered_map<int, int>& reverse_continuous_node_mapping) {
    std::vector<int> complement;
    std::uniform_int_distribution<int> generator_uniform_distribution{0, (int)(graph->GetNodeSet().size() - 1)};
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

std::vector<int> ABM::GetGeneratorNodes(Graph* graph, const std::unordered_map<int, int>& reverse_continuous_node_mapping) {
    std::vector<int> generator_nodes;
    std::uniform_int_distribution<int> generator_uniform_distribution{0, (int)(graph->GetNodeSet().size() - 1)};
    int num_generator_nodes = 1;
    /* std::random_device rand_dev; */
    /* std::minstd_rand generator{rand_dev()}; */
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
        // MARK: currently assumes single generator node
        int max_neighborhood_size = this->neighborhood_sample;
        std::set<int> visited;
        pcg_extras::seed_seq_from<std::random_device> rand_dev;
        pcg32 generator(rand_dev);
        for(size_t i = 0; i < generator_nodes.size(); i ++) {
            // get the 1-hop first
            int generator_node = generator_nodes.at(i);
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
            /* int num_one_hop_nodes = one_and_two_hop_neighborhood_map[1].size(); */
            int num_one_hop_nodes = current_one_hop_neighborhood.size();
            if (num_one_hop_nodes > max_neighborhood_size) {
                std::set<int> sampled_one_hop_neighborhood;
                std::sample(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), std::inserter(sampled_one_hop_neighborhood, sampled_one_hop_neighborhood.begin()), max_neighborhood_size, generator);
                current_one_hop_neighborhood = std::vector(sampled_one_hop_neighborhood.begin(), sampled_one_hop_neighborhood.end());
            }
            // until here should be fast even if sampled since 1-hop node list is constant time sort of
            visited.insert(generator_node);
            std::shuffle(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), generator);
            std::copy(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), std::back_inserter(one_and_two_hop_neighborhood_map[1]));
            for(size_t j = 0; j < current_one_hop_neighborhood.size(); j ++) {
                if (graph->GetOutDegree(current_one_hop_neighborhood.at(j)) > 0) {
                    for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(current_one_hop_neighborhood.at(j))) {
                        if (!visited.contains(outgoing_neighbor)) {
                            one_and_two_hop_neighborhood_map[2].push_back(outgoing_neighbor);
                        }
                        if (one_and_two_hop_neighborhood_map[2].size() == max_neighborhood_size) {
                            return one_and_two_hop_neighborhood_map;
                        }
                    }
                }
                if (graph->GetInDegree(current_one_hop_neighborhood.at(j)) > 0) {
                    for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(current_one_hop_neighborhood.at(j))) {
                        if (!visited.contains(incoming_neighbor)) {
                            one_and_two_hop_neighborhood_map[2].push_back(incoming_neighbor);
                        }
                        if (one_and_two_hop_neighborhood_map[2].size() == max_neighborhood_size) {
                            return one_and_two_hop_neighborhood_map;
                        }
                    }
                }
            }
        }
    }
    return one_and_two_hop_neighborhood_map;
}

void ABM::SampleFromNeighborhoods(std::unordered_map<int, std::vector<int>>& one_and_two_hop_neighborhood_map, int neighborhood_size_threshold, int max_neighborhood_size) {
    /* std::cerr << "sample function start" << std::endl; */
    /* std::random_device rand_dev; */
    /* std::minstd_rand generator{rand_dev()}; */
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    for(int i = 1; i < 3; i ++) {
        int current_neighborhood_size = one_and_two_hop_neighborhood_map[i].size();
        if (current_neighborhood_size > neighborhood_size_threshold) {
            if (current_neighborhood_size > max_neighborhood_size) {
                std::unordered_set<int> sampled_neighborhood;
                if (current_neighborhood_size < max_neighborhood_size * 10) {
                    // current neighborhood size is 10001, max_neighborhood size is 10000, and 10x is 100,000
                    /* std::cerr << "std::sample" << std::endl; */
                    std::sample(one_and_two_hop_neighborhood_map.at(i).begin(), one_and_two_hop_neighborhood_map.at(i).end(), std::inserter(sampled_neighborhood, sampled_neighborhood.begin()), max_neighborhood_size, generator);
                } else {
                    /* std::cerr << "random repeated sample" << std::endl; */
                    std::uniform_int_distribution<int> generator_uniform_distribution{0, (int)(one_and_two_hop_neighborhood_map.at(i).size() - 1)};
                    while (sampled_neighborhood.size() != (size_t) max_neighborhood_size) {
                        int new_index = generator_uniform_distribution(generator);
                        sampled_neighborhood.insert(one_and_two_hop_neighborhood_map.at(i).at(new_index));
                    }
                }
                /* assert(sampled_neighborhood.size() == max_neighborhood_size); */
                one_and_two_hop_neighborhood_map[i] = std::vector(sampled_neighborhood.begin(), sampled_neighborhood.end());
                /* assert(one_and_two_hop_neighborhood_map.at(i).size() == max_neighborhood_size); */
            }
        }
    }
    /* std::cerr << "sample function end" << std::endl; */
}

std::vector<int> ABM::GetNeighborhood(Graph* graph, const std::vector<int>& generator_nodes, const std::unordered_map<int, int>& reverse_continuous_node_mapping) {
    std::vector<int> neighborhood;
    std::set<int> visited;
    std::uniform_int_distribution<int> generator_uniform_distribution{0, (int)(graph->GetNodeSet().size() - 1)};
    int num_hops = 1;
    /* std::cout << "getting 2-hop neighborhood from: " << generator_nodes.at(0) << std::endl; */
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
void ABM::InitializeAuthors(Graph* graph, const std::unordered_map<int, int>& continuous_node_mapping, std::unordered_map<int, std::vector<int>>& author_to_publication_map, std::unordered_map<int, std::vector<int>>& number_published_to_author_map, std::unordered_map<int, int>& author_remaining_years_map) {
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

void ABM::AssignAuthor(Graph* graph, int new_node, std::unordered_map<int, std::vector<int>>& author_to_publication_map, std::unordered_map<int, std::vector<int>>& number_published_to_author_map, std::unordered_map<int, int>& author_remaining_years_map) {
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

void ABM::AgeAuthors(std::unordered_map<int, std::vector<int>>& author_to_publication_map, std::unordered_map<int, std::vector<int>>& number_published_to_author_map, std::unordered_map<int, int>& author_remaining_years_map) {
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
    /* this->timing_map[current_year].push_back({"initialize indegree array", seconds_elapsed}); */
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
    /* reading input edgelist, nodelist, outdegree bag, recency probabilities */
    Graph* graph = new Graph(this->edgelist, this->nodelist);
    /* graph->PrintGraph(); */
    this->WriteToLogFile("loaded graph", Log::info);
    this->InitializeFitness(graph);
    this->WriteToLogFile("initialized fitness for the seed graph", Log::debug);

    /* node ids to continous integer from 0 */
    std::unordered_map<int, int> continuous_node_mapping = this->BuildContinuousNodeMapping(graph);

    this->WriteToLogFile("forward built", Log::debug);
    /* continous integer from 0 to node ids*/
    std::unordered_map<int, int> reverse_continuous_node_mapping = this->ReverseMapping(continuous_node_mapping);
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
    /* int* citations = new int[250]; // maximum size container with 250 assumption for now */
    /* this->WriteToLogFile("alloc 8", Log::debug); */

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

    /* std::unordered_map<int, std::set<int>> author_to_publication_map; // map[1] = {98} meaning author 1 published the node 98 */
    /* std::unordered_map<int, std::set<int>> number_published_to_author_map; // map[1] = {1, 3, 4} meaning authors 1, 3, and 3 have 1 publications each */
    /* std::unordered_map<int, int> author_remaining_years_map; // map[1] = 2 meaning  author 1 has 2 more years left to live */
    /* this->InitializeAuthors(graph, continuous_node_mapping, author_to_publication_map, number_published_to_author_map, author_remaining_years_map); */
    /* graph->SetIntAttribute("fitness_peak_value", 4920089, 1000); */
    /* graph->SetIntAttribute("fitness_lag_duration", 4920089, 0); */
    /* graph->SetIntAttribute("fitness_peak_duration", 4920089, 1000); */

    std::vector<int> new_nodes_vec;
    std::set<int> same_year_source_nodes;
    std::vector<std::pair<int, int>> new_edges_vec;
    std::unordered_map<int, double> cached_results;
    for(int i = 0; i < 10000; i ++) {
        cached_results[i] = pow(i, this->gamma) + 1;
    }
    /* std::unordered_map<int, std::unordered_map<std::string, int> timing_map; */
    /* std::chrono::steady_clock::time_point<std::chrono::steady_clock> start_time = std::chrono::steady_clock::now(); */
    /* auto seconds_elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_time); */
    Eigen::setNbThreads(1);
    for (int current_year = start_year; current_year < start_year + this->num_cycles; current_year ++) {
        this->prev_time = std::chrono::steady_clock::now();
        /* std::cout << "new year" << std::endl; */
        int current_graph_size = graph->GetNodeSet().size();
        this->WriteToLogFile("current year is: " + std::to_string(current_year) + " and the graph is " + std::to_string(current_graph_size) + " nodes large", Log::info);
        this->FillInDegreeArr(graph, continuous_node_mapping, in_degree_arr);
        this->LogTime(current_year, "Fill in-degree array");
        this->WriteToLogFile("indegree for current year filled", Log::debug);
        this->FillFitnessArr(graph, continuous_node_mapping, current_year, fitness_arr);
        this->LogTime(current_year, "Fill fitness array");
        this->WriteToLogFile("fitness for current year filled", Log::debug);
        this->FillRecencyArr(graph, reverse_continuous_node_mapping, current_year, recency_arr);
        this->LogTime(current_year, "Fill recency array");
        this->WriteToLogFile("recency for current year calculated", Log::debug);
        this->CalculateScores(cached_results, in_degree_arr, pa_arr, current_graph_size);
        this->LogTime(current_year, "Process in-degree array");
        this->WriteToLogFile("indegree gammad", Log::debug);
        this->CalculateScores(cached_results, fitness_arr, fit_arr, current_graph_size);
        this->LogTime(current_year, "Process fitness array");
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
        this->LogTime(current_year, "Create new node ids");
        this->WriteToLogFile("all new nodes initialized with years and mapped", Log::debug);
        this->FillSameYearSourceNodes(same_year_source_nodes, new_nodes_vec.size());
        this->LogTime(current_year, "Pick same year nodes");

        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            int new_node = new_nodes_vec[i];
            std::vector<int> generator_nodes = this->GetGeneratorNodes(graph, reverse_continuous_node_mapping);
            this->UpdateGraphAttributesGeneratorNodes(graph, new_node, generator_nodes);
        }
        this->LogTime(current_year, "Pick generator nodes");
        std::vector<std::pair<int, int>> original_neighborhood_sizes_map(new_nodes_vec.size());

        std::vector<std::pair<std::string, int>> parallel_stage_time_vec;
        /*
        Eigen::setNbThreads(this->num_processors);
        Eigen::SparseMatrix<int, 0, int> adjacency_matrix = Eigen::SparseMatrix<int, 0, int>(current_graph_size, current_graph_size);
        std::vector<Eigen::Triplet<int>> triplet_list;
        for(int row = 0; row < current_graph_size; row ++) {
            int row_current_node_id = reverse_continuous_node_mapping[row];
            if (graph->GetOutDegree(row_current_node_id) > 0) {
                for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(row_current_node_id)) {
                    int col_integer_node_id = continuous_node_mapping.at(outgoing_neighbor);
                    triplet_list.push_back(Eigen::Triplet<int>(row, col_integer_node_id, 1));
                }
            }
            if (graph->GetInDegree(row_current_node_id) > 0) {
                for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(row_current_node_id)) {
                    int col_integer_node_id = continuous_node_mapping.at(incoming_neighbor);
                    triplet_list.push_back(Eigen::Triplet<int>(row, col_integer_node_id, 1));
                }
            }
        }
        adjacency_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
        this->LogTime(current_year, "Build adjacency matrix");
        Eigen::SparseMatrix<int, 0, int> two_hop_matrix = (adjacency_matrix * adjacency_matrix).pruned();
        this->LogTime(current_year, "Build 2-hop matrix");
        Eigen::setNbThreads(1);
        */
        /* std::cerr<< "process nodes start" <<std::endl; */
        #pragma omp parallel for reduction(merge_int_pair_vecs: new_edges_vec) reduction(merge_str_int_pair_vecs: parallel_stage_time_vec)
        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            std::chrono::time_point<std::chrono::steady_clock> local_prev_time = std::chrono::steady_clock::now();
            std::vector<std::pair<int, int>> local_new_edges_vec;
            std::vector<std::pair<std::string, int>> local_parallel_stage_time_vec;


            int citations[250]; // out-degree assumed to be max 249
            this->WriteToLogFile("starting node " + std::to_string(i) + "/" + std::to_string(new_nodes_vec.size()), Log::debug);
            int new_node = new_nodes_vec[i];
            /* this->AssignAuthor(graph, new_node, author_to_publication_map, number_published_to_author_map, author_remaining_years_map); */
            int weight_arr_index = continuous_node_mapping[new_node] - initial_graph_size;
            double pa_weight = pa_weight_arr[weight_arr_index];
            double rec_weight = rec_weight_arr[weight_arr_index];
            double fit_weight = fit_weight_arr[weight_arr_index];
            double alpha = alpha_arr[weight_arr_index];
            std::vector<int> generator_nodes = this->GetGraphAttributesGeneratorNodes(graph, new_node);
            std::unordered_map<int, std::vector<int>> one_and_two_hop_neighborhood_map = this->GetOneAndTwoHopNeighborhood(graph, generator_nodes, reverse_continuous_node_mapping);
            /* std::unordered_map<int, std::vector<int>> one_and_two_hop_neighborhood_map = this->GetOneAndTwoHopNeighborhood(graph, generator_nodes, reverse_continuous_node_mapping, 10000); */
            /* std::unordered_map<int, std::vector<int>> one_and_two_hop_neighborhood_map = this->GetOneAndTwoHopNeighborhoodFromMatrix(graph, current_graph_size, two_hop_matrix, generator_nodes, continuous_node_mapping, reverse_continuous_node_mapping); */
            original_neighborhood_sizes_map[i] = {one_and_two_hop_neighborhood_map[1].size(), one_and_two_hop_neighborhood_map[2].size()};
            /* std::cerr << "got generator neighborhoods" << std::endl; */

            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "retrieve one and two hop neighborhoods");

            /* this->SampleFromNeighborhoods(one_and_two_hop_neighborhood_map, 10000, 10000); */
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
            /* std::cerr << "finished calculating num nodes to cite" << std::endl; */

            int num_actually_cited = 0;
            if (same_year_citation) {
                num_actually_cited += this->MakeSameYearCitations(new_nodes_vec.size(), reverse_continuous_node_mapping, citations, current_graph_size);
            }

            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make same year citations");

            /* std::cerr << "finished same year" << std::endl; */

            num_actually_cited += this->MakeCitations(graph, continuous_node_mapping, current_year, one_and_two_hop_neighborhood_map[1], citations + num_actually_cited, pa_arr, recency_arr, fit_arr, pa_weight, rec_weight, fit_weight, current_graph_size, num_citations_inside);

            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make 1-hop citations");
            /* std::cerr << "finished 1-hop" << std::endl; */

            num_actually_cited += this->MakeCitations(graph, continuous_node_mapping, current_year, one_and_two_hop_neighborhood_map[2], citations + num_actually_cited, pa_arr, recency_arr, fit_arr, pa_weight, rec_weight, fit_weight, current_graph_size, num_citations_outside);

            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make 2-hop citations");
            /* std::cerr << "finished 2-hop" << std::endl; */

            num_actually_cited += this->MakeUniformRandomCitations(graph, reverse_continuous_node_mapping, generator_nodes, citations, num_actually_cited, num_fully_random_cited);
            /* std::cerr << "finished uniform" << std::endl; */

            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make random citations");

            for(size_t j = 0; j < generator_nodes.size(); j ++) {
                local_new_edges_vec.push_back({new_node, generator_nodes[j]});
            }
            for(int j = 0; j < num_actually_cited; j ++) {
                local_new_edges_vec.push_back({new_node, citations[j]});
            }
            new_edges_vec.insert(new_edges_vec.end(), local_new_edges_vec.begin(), local_new_edges_vec.end());
            /* std::cerr << "finished recroding edges" << std::endl; */

            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "record edges");
            parallel_stage_time_vec.insert(parallel_stage_time_vec.end(), local_parallel_stage_time_vec.begin(), local_parallel_stage_time_vec.end());
            /* std::cerr << "updating time log" << std::endl; */
        } // end of omp parallel loop
        /* std::cerr<< "process nodes end" <<std::endl; */
        std::map<std::string, int> per_stage_time_map;
        for(size_t i = 0; i < parallel_stage_time_vec.size(); i ++) {
            per_stage_time_map[parallel_stage_time_vec[i].first] +=  parallel_stage_time_vec[i].second;
        }
        /* this->timing_map[current_year].push_back({"find one and two hop neighborhoods", per_stage_time_map["find one and two hop neighborhoods"]}); */
        this->LogTime(current_year, "find one and two hop neighborhoods", per_stage_time_map["retrieve one and two hop neighborhoods"]);
        /* this->timing_map[current_year].push_back({"make same year citations", per_stage_time_map["make same year citations"]}); */
        this->LogTime(current_year, "make same year citations", per_stage_time_map["make same year citations"]);
        /* this->timing_map[current_year].push_back({"make 1-hop citations", per_stage_time_map["make 1-hop citations"]}); */
        this->LogTime(current_year, "make 1-hop citations", per_stage_time_map["make 1-hop citations"]);
        /* this->timing_map[current_year].push_back({"make 2-hop citations", per_stage_time_map["make 2-hop citations"]}); */
        this->LogTime(current_year, "make 2-hop citations", per_stage_time_map["make 2-hop citations"]);
        /* this->timing_map[current_year].push_back({"make random citations", per_stage_time_map["make random citations"]}); */
        this->LogTime(current_year, "make random citations", per_stage_time_map["make random citations"]);
        /* this->timing_map[current_year].push_back({"record edges", per_stage_time_map["record edges"]}); */
        this->LogTime(current_year, "record edges", per_stage_time_map["record edges"]);
        /* std::chrono::time_point<std::chrono::steady_clock> current_time = std::chrono::steady_clock::now(); */
        /* this->prev_time = current_time; */


        /* this->WriteToLogFile("edges saved to vector", Log::debug); */
        for(size_t i = 0; i < new_edges_vec.size(); i ++) {
            int new_node = new_edges_vec[i].first;
            int destination_id = new_edges_vec[i].second;
            graph->AddEdge({new_node, destination_id});
        }
        this->LogTime(current_year, "Add edges to graph");

        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            int new_node = new_nodes_vec[i];
            graph->SetIntAttribute("one_hop_neighborhood_size", new_node, original_neighborhood_sizes_map[i].first);
            graph->SetIntAttribute("two_hop_neighborhood_size", new_node, original_neighborhood_sizes_map[i].second);
        }

        this->LogTime(current_year, "Update graph attributes (neighborhood sizes)");

        this->WriteToLogFile("edges saved to graph", Log::debug);
        this->AssignPeakFitnessValues(graph, new_nodes_vec);
        this->WriteToLogFile("assigned peak fitness for new nodes", Log::debug);
        this->AssignFitnessLagDuration(graph, new_nodes_vec);
        this->WriteToLogFile("assigned fitness lag duration for new nodes", Log::debug);
        this->AssignFitnessPeakDuration(graph, new_nodes_vec);
        this->WriteToLogFile("assigned fitness peak duration for new nodes", Log::debug);
        this->PlantNodes(graph, new_nodes_vec, current_year - start_year + 1);
        this->LogTime(current_year, "Assign fitness values to new nodes");
        /* this->AgeAuthors(author_to_publication_map, number_published_to_author_map, author_remaining_years_map); */
        new_nodes_vec.clear();
        new_edges_vec.clear();
        /* this->WriteToLogFile("writing temp graph", Log::info); */
        /* graph->WriteGraph(this->output_file + "_" + std::to_string(current_year)); */
        /* this->WriteToLogFile("writing temp graph and aux", Log::info); */
        /* this->UpdateGraphAttributesWeights(graph, initial_next_node_id, pa_weight_arr, rec_weight_arr, fit_weight_arr, final_graph_size - initial_graph_size); */
        /* this->UpdateGraphAttributesAlphas(graph, initial_next_node_id, alpha_arr, final_graph_size - initial_graph_size); */
        /* this->UpdateGraphAttributesOutDegrees(graph, initial_next_node_id, out_degree_arr, final_graph_size - initial_graph_size); */
        /* for(auto const& node_id : graph->GetNodeSet()) { */
        /*     graph->SetIntAttribute("in_degree", node_id, graph->GetInDegree(node_id)); */
        /*     graph->SetIntAttribute("out_degree", node_id, graph->GetOutDegree(node_id)); */
        /* } */
        /* graph->WriteGraph(this->output_file + "_" + std::to_string(current_year)); */
        /* graph->WriteAttributes(this->auxiliary_information_file + "_" + std::to_string(current_year)); */
        /* this->LogTime(current_year, "Write temporary graph and aux files"); */
    }
    /* this->WriteTimingFile(start_year, start_year + this->num_cycles); */

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
