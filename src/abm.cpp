#include "abm.h"

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

/* void ABM::AssignPeakFitnessValues(Graph* graph, const std::set<int>& nodeset) { */
/*     for(auto const& node : nodeset) { */
/*         double fitness_uniform = this->fitness_value_uniform_distribution(this->generator); */
/*         double adjusted_alpha = this->fitness_alpha + 1; */
/*         double base_left = (pow(this->fitness_value_max, adjusted_alpha) - pow(this->fitness_value_min, adjusted_alpha)) * fitness_uniform; */
/*         double base_right = pow(this->fitness_value_min, adjusted_alpha); */
/*         double exponent = 1.0/adjusted_alpha; */
/*         int fitness_power = pow(base_left + base_right ,exponent); */
/*         graph->SetAttribute("peak_fitness_value", node, fitness_power); */
/*     } */
/* } */

/* void ABM::AssignFitnessLagDuration(Graph* graph, const std::set<int>& nodeset) { */
/*     for(auto const& node : nodeset) { */
/*         int fitness_lag_uniform = this->fitness_lag_duration_uniform_distribution(this->generator); */
/*         graph->SetAttribute("fitness_lag_duration", node, fitness_lag_uniform); */
/*     } */
/* } */

/* void ABM::AssignFitnessPeakDuration(Graph* graph, const std::set<int>& nodeset) { */
/*     for(auto const& node : nodeset) { */
/*         int fitness_peak_uniform = this->fitness_peak_duration_uniform_distribution(this->generator); */
/*         graph->SetAttribute("fitness_peak_duration", node, fitness_peak_uniform); */
/*     } */
/* } */

void ABM::InitializeFitness(Graph* graph) {
    this->AssignPeakFitnessValues(graph, graph->GetNodeSet());
    this->AssignFitnessLagDuration(graph, graph->GetNodeSet());
    this->AssignFitnessPeakDuration(graph, graph->GetNodeSet());
}

void ABM::FillFitnessArr(Graph* graph, const std::map<int, int>& continuous_node_mapping, int current_year, int* fitness_arr) {
    for(auto const& node : graph->GetNodeSet()) {
        int peak_fitness_value = graph->GetAttribute("peak_fitness_value", node);
        int fitness_lag_duration = graph->GetAttribute("fitness_lag_duration", node);
        int fitness_peak_duration = graph->GetAttribute("fitness_peak_duration", node);
        int published_year = graph->GetAttribute("year", node);
        int continuous_index = continuous_node_mapping.at(node);
        if (published_year + fitness_lag_duration >= current_year) {
            fitness_arr[continuous_index] = 1;
        } else if (published_year + fitness_lag_duration + fitness_peak_duration >= current_year) {
            fitness_arr[continuous_index] =peak_fitness_value;
        } else {
            double decayed_fitness_value = peak_fitness_value / pow(current_year - published_year - fitness_lag_duration - fitness_peak_duration + 1, this->fitness_decay_alpha);
            fitness_arr[continuous_index] = decayed_fitness_value;
        }
    }
}


void ABM::FillRecencyArr(Graph* graph, const std::map<int, int>& continuous_node_mapping, int current_year, double* recency_arr) {
    for(auto const& node : graph->GetNodeSet()) {
        int current_published_year = graph->GetAttribute("year", node);
        int continuous_index = continuous_node_mapping.at(node);
        int year_diff = current_year - current_published_year;
        recency_arr[continuous_index] = this->recency_probabilities_map[year_diff];
    }
    double recency_arr_sum = 0.0;
    #pragma omp parallel for reduction(+:recency_arr_sum)
    for(size_t i = 0; i < graph->GetNodeSet().size(); i ++) {
        recency_arr_sum += recency_arr[i];
    }
    #pragma omp parallel for
    for(size_t i = 0; i < graph->GetNodeSet().size(); i ++) {
        recency_arr[i] /= recency_arr_sum;
    }
}


int ABM::GetMaxYear(Graph* graph) {
    int max_year = -1;
    bool is_first = true;
    for(auto const& node : graph->GetNodeSet()) {
        int current_node_year = graph->GetAttribute("year", node);
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
    if(this->alpha == -1) {
        for(int i = 0; i < len; i ++) {
            int alpha_uniform = this->alpha_uniform_distribution(generator);
            alpha_arr[i] = alpha_uniform;
        }
    } else {
        for(int i = 0; i < len; i ++) {
            alpha_arr[i] = this->alpha;
        }
    }
}

void ABM::PopulateWeightArrs(double* pa_weight_arr, double* rec_weight_arr, double* fit_weight_arr, int len) {
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

void ABM::PopulateOutDegreeArr(int* out_degree_arr, int len) {
    std::uniform_int_distribution<int> outdegree_index_uniform_distribution{0, this->out_degree_bag_vec.size() - 1};
    /* std::cout << "randomly sampling from 0 to " << this->out_degree_bag_vec.size() << " to get out-degree" << std::endl; */
    for(int i = 0; i < len; i ++) {
        int index_uniform = outdegree_index_uniform_distribution(generator);
        /* std::cout << "chose index " << index_uniform << std::endl; */
        out_degree_arr[i] = this->out_degree_bag_vec[index_uniform];
    }
}

void ABM::CalculateScores(int* src_arr, double* dst_arr, int len) {
    double sum = 0;
    #pragma omp parallel for
    for(int i = 0; i < len; i ++) {
        dst_arr[i] = pow(src_arr[i], this->gamma) + 1;
    }
    #pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < len; i ++) {
        sum += dst_arr[i];
    }
    #pragma omp parallel for
    for(int i = 0; i < len; i ++) {
        dst_arr[i] /= sum;
    }
}

int ABM::MakeCitations(int* citations, double* score_arr, double* random_weight_arr, int len, int num_citations) {
    int actual_num_cited = num_citations;
    if (len < num_citations) {
        actual_num_cited = len;
    }
    std::random_device rand_dev;
    std::mt19937 generator{rand_dev()};
    for(int i = 0; i < actual_num_cited; i ++) {
        std::discrete_distribution<int> int_discrete_distribution(score_arr, score_arr + len);
        int current_citation = int_discrete_distribution(generator);
        citations[i] = current_citation;
        score_arr[i] = 0.0;
    }
    return actual_num_cited;

    /* double* random_weight_arr = new double[len]; */
    // BEGIN code that scales with size of array
    /* #pragma omp parallel for */
    /* for(int i = 0; i < len; i ++) { */
    /*     double wrs_uniform = this->wrs_uniform_distribution(this->generator); */
    /*     if (score_arr[i] != 0) { */
    /*         random_weight_arr[i] = pow(wrs_uniform, 1.0/score_arr[i]); */
    /*     } else { */
    /*         random_weight_arr[i] = 0; */
    /*     } */
    /* } */

    /* std::vector<std::pair<int, int>> element_index_vec; */
    /* for (size_t i = 0; i < len; ++i) { */
    /*     element_index_vec.push_back({random_weight_arr[i], i}); */
    /* } */
    /* int actual_num_cited = num_citations; */
    /* if (len < num_citations) { */
    /*     actual_num_cited = len; */
    /* } */
    /* /1* std::partial_sort(element_index_vec.begin(), element_index_vec.begin() + actual_num_cited, element_index_vec.end(), std::greater<std::pair<int, int>>()); *1/ */
    /* std::sort(element_index_vec.begin(), element_index_vec.end(), [](auto& left, auto& right) { */
    /*     return left.first < right.first; */
    /* }); */
    /* for (int i = 0; i < actual_num_cited; i ++) { */
    /*     citations[i] = element_index_vec[i].second; */
    /* } */
    /* return actual_num_cited; */
    // END code that scales with size of array
}


int ABM::ZeroOutNonNeighbors(Graph* graph, const std::map<int, int>& continuous_node_mapping, const std::map<int, int>& reverse_continuous_node_mapping, double* current_score_arr) {
    std::uniform_int_distribution<int> generator_uniform_distribution{0, graph->GetNodeSet().size() - 1};
    int continuous_generator_node = generator_uniform_distribution(this->generator);
    int generator_node = reverse_continuous_node_mapping.at(continuous_generator_node);
    /* std::cout << "generator node currentnly is: " << "(continuous) " << continuous_generator_node << " and (actual node id) " << generator_node << std::endl; */
    if(graph->GetOutDegree(generator_node) > 0) {
        for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(generator_node)) {
            int continuous_node_id = continuous_node_mapping.at(outgoing_neighbor);
            current_score_arr[continuous_node_id] = 0.0;
        }
    }
    if(graph->GetInDegree(generator_node) > 0) {
        for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(generator_node)) {
            /* std::cout << "has incoming neighbor " << incoming_neighbor << std::endl; */
            int continuous_node_id = continuous_node_mapping.at(incoming_neighbor);
            /* std::cout << "incoming neighbor continuous id is " << continuous_node_id << std::endl; */
            current_score_arr[continuous_node_id] = 0.0;
        }
    }
    double score_arr_sum = 0.0;
    #pragma omp parallel for reduction(+:score_arr_sum)
    for(size_t i = 0; i < graph->GetNodeSet().size(); i ++) {
        score_arr_sum += current_score_arr[i];
    }
    #pragma omp parallel for
    for(size_t i = 0; i < graph->GetNodeSet().size(); i ++) {
       current_score_arr[i] /= score_arr_sum;
    }
    return generator_node;
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

    int start_year = this->GetMaxYear(graph);
    int next_node_id = this->GetMaxNode(graph) + 1;

    /* get input to score arrays based on continuous_node_mapping */
    int initial_graph_size = graph->GetNodeSet().size();
    int final_graph_size = this->GetFinalGraphSize(graph);
    this->WriteToLogFile("final graph size is " + std::to_string(final_graph_size), Log::info);
    /* std::cout << "final graph size is: " << final_graph_size << std::endl; */
    int* in_degree_arr = new int[final_graph_size];
    this->WriteToLogFile("alloc 1", Log::info);
    int* fitness_arr = new int[final_graph_size];
    this->WriteToLogFile("alloc 2", Log::info);
    double* pa_arr = new double[final_graph_size];
    this->WriteToLogFile("alloc 3", Log::info);
    double* fit_arr = new double[final_graph_size];
    this->WriteToLogFile("alloc 4", Log::info);
    double* recency_arr = new double[final_graph_size];
    this->WriteToLogFile("alloc 5", Log::info);
    double* random_weight_arr = new double[final_graph_size];
    this->WriteToLogFile("alloc 6", Log::info);
    double* current_score_arr = new double[final_graph_size];
    this->WriteToLogFile("alloc 7", Log::info);
    int* citations = new int[250]; // maximum size container with 250 assumption for now
    this->WriteToLogFile("alloc 8", Log::info);

    // the first new agent node has index 0 but is actually index initial_graph_size in the continuous mapping
    double* pa_weight_arr = new double[final_graph_size - initial_graph_size];
    this->WriteToLogFile("alloc 9", Log::info);
    double* rec_weight_arr = new double[final_graph_size - initial_graph_size];
    this->WriteToLogFile("alloc 10", Log::info);
    double* fit_weight_arr = new double[final_graph_size - initial_graph_size];
    this->WriteToLogFile("alloc 11", Log::info);
    double* alpha_arr = new double[final_graph_size - initial_graph_size];
    this->WriteToLogFile("alloc 12", Log::info);
    int* out_degree_arr = new int[final_graph_size - initial_graph_size];
    this->WriteToLogFile("alloc 13", Log::info);
    this->WriteToLogFile("allocated arrays", Log::debug);

    this->PopulateWeightArrs(pa_weight_arr, rec_weight_arr, fit_weight_arr, final_graph_size - initial_graph_size);
    this->WriteToLogFile("populated weight arrays", Log::debug);
    this->PopulateAlphaArr(alpha_arr, final_graph_size - initial_graph_size);
    this->WriteToLogFile("populated alpha array", Log::debug);
    this->PopulateOutDegreeArr(out_degree_arr, final_graph_size - initial_graph_size);
    this->WriteToLogFile("populated out degree array", Log::debug);
    /* std::cout << "outdegree bag contains:"; */
    /* for(int i = 0; i < final_graph_size - initial_graph_size; i ++) { */
    /*     std::cout << " " << out_degree_arr[i]; */
    /* } */
    /* std::cout << std::endl; */
    std::map<int, std::set<int>> author_to_publication_map; // map[1] = {98} meaning author 1 published the node 98
    std::map<int, std::set<int>> number_published_to_author_map; // map[1] = {1, 3, 4} meaning authors 1, 3, and 3 have 1 publications each
    std::map<int, int> author_remaining_years_map; // map[1] = 2 meaning  author 1 has 2 more years left to live
    /* this->InitializeAuthors(graph, continuous_node_mapping, author_to_publication_map, number_published_to_author_map, author_remaining_years_map); */

    std::vector<int> new_nodes_vec;
    std::vector<std::pair<int, int>> new_edges_vec;
    for (int current_year = start_year; current_year < start_year + this->num_cycles; current_year ++) {
        /* std::cout << "current year: " << current_year << std::endl; */
        int current_graph_size = graph->GetNodeSet().size();
        this->WriteToLogFile("current year is: " + std::to_string(current_year) + " and the graph is " + std::to_string(current_graph_size) + " nodes large", Log::info);
        this->FillInDegreeArr(graph, continuous_node_mapping, in_degree_arr);
        this->WriteToLogFile("indegree for current year filled", Log::debug);
        this->FillFitnessArr(graph, continuous_node_mapping, current_year, fitness_arr);
        this->WriteToLogFile("fitness for current year filled", Log::debug);
        this->FillRecencyArr(graph, continuous_node_mapping, current_year, recency_arr);
        this->WriteToLogFile("recency for current year calculated", Log::debug);
        this->CalculateScores(in_degree_arr, pa_arr, current_graph_size);
        this->WriteToLogFile("indegree gammad", Log::debug);
        this->CalculateScores(fitness_arr, fit_arr, current_graph_size);
        this->WriteToLogFile("fitness gammad", Log::debug);

        /* print scores */
        /* std::cout << "["; */
        /* for(auto const& [node_id, continuous_id] : continuous_node_mapping) { */
        /*     std::cout << " (" << pa_arr[continuous_id] << ", " << recency_arr[continuous_id] << ", " << fit_arr[continuous_id] << ") "; */
        /* } */
        /* std::cout << "]" << std::endl; */

        /* initialize new nodes */
        int num_new_nodes = std::ceil(current_graph_size * this->growth_rate);
        this->WriteToLogFile("making " + std::to_string(num_new_nodes) + " nodes this year", Log::info);
        /* std::cout << "making " << num_new_nodes << " new nodes" << std::endl; */
        for(int i = 0; i < num_new_nodes; i ++) {
            continuous_node_mapping[next_node_id] = current_graph_size + i;
            reverse_continuous_node_mapping[current_graph_size + i] = next_node_id;
            new_nodes_vec.push_back(next_node_id);
            graph->SetAttribute("year", next_node_id, current_year);
            next_node_id ++;
        }
        this->WriteToLogFile("all new nodes initialized with years and mapped", Log::debug);

        /* for(auto const& new_node : new_nodes_set) { */
        /* for(std::set<int>::iterator it = new_nodes_set.begin(); it != new_nodes_set.end(); it++) { */
        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            this->WriteToLogFile("starting node " + std::to_string(i) + "/" + std::to_string(new_nodes_vec.size()), Log::debug);
            int new_node = new_nodes_vec[i];
            /* this->AssignAuthor(graph, new_node, author_to_publication_map, number_published_to_author_map, author_remaining_years_map); */
            int weight_arr_index = continuous_node_mapping[new_node] - initial_graph_size;
            /* std::cout << "new node made with: " << "(continuous) " << continuous_node_mapping[new_node] << " and actual node id: " << new_node << std::endl; */
            /* std::cout << "weight arr index of: " << weight_arr_index << std::endl; */
            double pa_weight = pa_weight_arr[weight_arr_index];
            double rec_weight = rec_weight_arr[weight_arr_index];
            double fit_weight = fit_weight_arr[weight_arr_index];
            double alpha = alpha_arr[weight_arr_index];
            #pragma omp parallel for
            for(int i = 0; i < current_graph_size; i ++) {
                current_score_arr[i] = pa_weight * pa_arr[i] + rec_weight * recency_arr[i] + fit_weight * fit_arr[i];
            }
            this->WriteToLogFile("weighted sum array done", Log::debug);
            // cite from the entire graph for outdegree * alpha
            /* std::cout << "assigned outdegree: " << out_degree_arr[weight_arr_index] << std::endl; */
            /* std::cout << "requsted in-neighbor outdegree: " << (int)(out_degree_arr[weight_arr_index] * (1 - alpha)) << std::endl; */
            int num_citations_outside = out_degree_arr[weight_arr_index] * (1 - alpha);
            int num_actually_cited = this->MakeCitations(citations, current_score_arr, random_weight_arr, current_graph_size, num_citations_outside);
            this->WriteToLogFile(std::to_string(num_actually_cited) + " cited within neighborhood", Log::debug);

            // cite within generator for outdegree * alpha
            int num_citations_inside = out_degree_arr[weight_arr_index] - num_actually_cited;
            int generator_node = this->ZeroOutNonNeighbors(graph, continuous_node_mapping, reverse_continuous_node_mapping, current_score_arr);
            this->WriteToLogFile("generator 1-hop zerod out", Log::debug);
            num_actually_cited += this->MakeCitations(citations + num_actually_cited, current_score_arr, random_weight_arr, current_graph_size, num_citations_inside);
            this->WriteToLogFile("outside cited, total citations: " + std::to_string(num_actually_cited), Log::debug);

            new_edges_vec.push_back({new_node, generator_node});
            for(int i = 0; i < num_actually_cited; i ++) {
                int destination_id = reverse_continuous_node_mapping[citations[i]];
                new_edges_vec.push_back({new_node, destination_id});
                /* std::cout << "adding" << new_node << "-" << destination_id << std::endl; */
            }
            this->WriteToLogFile("edge vec now has " + std::to_string(new_edges_vec.size()) + " edges", Log::debug);
        }

        this->WriteToLogFile("edges saved to vector", Log::debug);
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
        /* this->AgeAuthors(author_to_publication_map, number_published_to_author_map, author_remaining_years_map); */
        new_nodes_vec.clear();
        new_edges_vec.clear();
    }

    this->WriteToLogFile("finished sim", Log::info);
    graph->WriteGraph(this->output_file);
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
    delete[] citations;
    delete[] random_weight_arr;
    delete[] current_score_arr;

    delete graph;
    return 0;
}
