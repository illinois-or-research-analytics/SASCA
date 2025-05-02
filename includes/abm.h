#ifndef ABM_H
#define ABM_H

#include <cmath>
#include <chrono>
#include <condition_variable>
#include <random>
#include <thread>
#include <map>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <queue>
#include <omp.h>
#include "graph.h"


enum Log {info, debug, error = -1};

class ABM {
    public:
        ABM(std::string edgelist, std::string nodelist, std::string out_degree_bag, std::string recency_probabilities, double alpha, double growth_rate, int num_cycles, std::string output_file, std::string log_file, int num_processors, int log_level) : edgelist(edgelist), nodelist(nodelist), out_degree_bag(out_degree_bag), recency_probabilities(recency_probabilities), alpha(alpha), growth_rate(growth_rate), num_cycles(num_cycles), output_file(output_file), log_file(log_file), num_processors(num_processors), log_level(log_level) {
            if(this->log_level > -1) {
                this->start_time = std::chrono::steady_clock::now();
                this->log_file_handle.open(this->log_file);
            }
            this->num_calls_to_log_write = 0;
            this->ReadOutDegreeBag();
            this->ReadRecencyProbabilities();
            omp_set_num_threads(this->num_processors);
        };

        ~ABM() {
            if(this->log_level > -1) {
                this->log_file_handle.close();
            }
        }

        int main();
        int WriteToLogFile(std::string message, Log message_type);
        void ReadOutDegreeBag();
        void ReadRecencyProbabilities();
        std::map<int, int> BuildContinuousNodeMapping(Graph* graph);
        std::map<int, int> ReverseMapping(std::map<int, int> mapping);
        int GetFinalGraphSize(Graph* graph);
        void FillInDegreeArr(Graph* graph, const std::map<int, int>& continuous_node_mapping, int* in_degree_arr);
        /* void AssignPeakFitnessValues(Graph* graph, const std::set<int>& nodeset); */
        /* void AssignFitnessLagDuration(Graph* graph, const std::set<int>& nodeset); */
        /* void AssignFitnessPeakDuration(Graph* graph, const std::set<int>& nodeset); */
        void InitializeFitness(Graph* graph);
        void FillFitnessArr(Graph* graph, const std::map<int, int>& continuous_node_mapping, int current_year, int* fitness_arr);
        void FillRecencyArr(Graph* graph, const std::map<int, int>& continuous_node_mapping, int current_year, double* recency_arr);
        void PopulateWeightArrs(double* pa_weight_arr, double* rec_weight_arr, double* fit_weight_arr, int len);
        void PopulateAlphaArr(double* alpha_arr, int len);
        int GetMaxYear(Graph* graph);
        int GetMaxNode(Graph* graph);
        void PopulateOutDegreeArr(int* out_degree_arr, int len);
        void CalculateScores(int* src_arr, double* dst_arr, int len);
        int MakeCitations(int* citations, double* score_arr, double* random_weight_arr, int len, int num_citations);
        int ZeroOutNonNeighbors(Graph* graph, const std::map<int, int>& continuous_node_mapping, const std::map<int, int>& reverse_continuous_node_mapping, double* current_score_arr);
        /* void InitializeAuthors(Graph* graph, const std::map<int, int>& continuous_node_mapping, std::map<int, std::vector<int>>& author_to_publication_map, std::map<int, std::vector<int>>& number_published_to_author_map, std::map<int, int>& author_remaining_years_map); */
        /* void AgeAuthors(std::map<int, std::vector<int>>& author_to_publication_map, std::map<int, std::vector<int>>& number_published_to_author_map, std::map<int, int>& author_remaining_years_map); */

        template<typename T>
        void AssignFitnessLagDuration(Graph* graph, const T& container) {
            for(auto const& node : container) {
                int fitness_lag_uniform = this->fitness_lag_duration_uniform_distribution(this->generator);
                graph->SetAttribute("fitness_lag_duration", node, fitness_lag_uniform);
            }
        }

        template<typename T>
        void AssignFitnessPeakDuration(Graph* graph, const T& container) {
            for(auto const& node : container) {
                int fitness_peak_uniform = this->fitness_peak_duration_uniform_distribution(this->generator);
                graph->SetAttribute("fitness_peak_duration", node, fitness_peak_uniform);
            }
        }

        template<typename T>
        void AssignPeakFitnessValues(Graph* graph, const T& container) {
            for(auto const& node : container) {
                double fitness_uniform = this->fitness_value_uniform_distribution(this->generator);
                double adjusted_alpha = this->fitness_alpha + 1;
                double base_left = (pow(this->fitness_value_max, adjusted_alpha) - pow(this->fitness_value_min, adjusted_alpha)) * fitness_uniform;
                double base_right = pow(this->fitness_value_min, adjusted_alpha);
                double exponent = 1.0/adjusted_alpha;
                int fitness_power = pow(base_left + base_right ,exponent);
                graph->SetAttribute("peak_fitness_value", node, fitness_power);
            }
        }


    protected:
        std::string edgelist;
        std::string nodelist;
        std::string out_degree_bag;
        std::string recency_probabilities;
        double alpha;
        double growth_rate;
        int num_cycles;
        std::string output_file;
        std::string log_file;
        int num_processors;
        int log_level;
        std::chrono::steady_clock::time_point start_time;
        std::ofstream log_file_handle;
        int num_calls_to_log_write;
        std::random_device rand_dev;
        std::mt19937 generator{rand_dev()};
        const int fitness_value_min = 1;
        const int fitness_value_max = 1000;
        const int fitness_lag_duration_min = 1;
        const int fitness_lag_duration_max = 7;
        const int fitness_peak_duration_min = 1;
        const int fitness_peak_duration_max = 7;
        const int fitness_alpha = -3;
        const int fitness_decay_alpha = 3;
        const int gamma = 3;
        const int max_author_lifetime = 30;
        const int k = 2;
        int next_author_id = 0;
        std::uniform_real_distribution<double> fitness_value_uniform_distribution{0, 1};
        std::uniform_real_distribution<double> weights_uniform_distribution{0, 1};
        std::uniform_real_distribution<double> wrs_uniform_distribution{0, 1};
        std::uniform_real_distribution<double> alpha_uniform_distribution{0, 1};
        std::uniform_int_distribution<int> fitness_lag_duration_uniform_distribution{fitness_lag_duration_min, fitness_lag_duration_max};
        std::uniform_int_distribution<int> fitness_peak_duration_uniform_distribution{fitness_peak_duration_min, fitness_peak_duration_max};
        std::vector<int> out_degree_bag_vec;
        std::map<int, double> recency_probabilities_map;
};

#endif
