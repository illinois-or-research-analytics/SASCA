#include "graph.h"


Graph::Graph(std::string edgelist, std::string nodelist): edgelist(edgelist), nodelist(nodelist) {
    this->ParseEdgelist();
    this->ParseNodelist();
}

void Graph::ParseEdgelist() {
    char delimiter = Graph::get_delimiter(this->edgelist);
    std::ifstream input_edgelist(this->edgelist);
    std::string line;
    while(std::getline(input_edgelist, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        std::string citing = current_line[0];
        if(citing[0] == '#') {
            continue;
        }
        int integer_citing = std::stoi(citing);
        int integer_cited = std::stoi(current_line[1]);
        this->AddEdge({integer_citing, integer_cited});
    }
}

void Graph::ParseNodelist() {
    char delimiter = Graph::get_delimiter(this->nodelist);
    std::ifstream input_nodelist(this->nodelist);
    std::string line;
    while(std::getline(input_nodelist, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        std::string node = current_line[0];
        if(node[0] == '#') {
            continue;
        }
        int integer_node = std::stoi(node);
        int integer_year = std::stoi(current_line[1]);
        this->SetIntAttribute("year", integer_node, integer_year);
        this->SetStringAttribute("type", integer_node, "seed");
    }
}

void Graph::SetIntAttribute(std::string attribute_key, int node, int attribute_value) {
    this->int_attribute_map[attribute_key][node] = attribute_value;
}

int Graph::GetIntAttribute(std::string attribute_key, int node) const {
    return this->int_attribute_map.at(attribute_key).at(node);
}

void Graph::SetStringAttribute(std::string attribute_key, int node, std::string attribute_value) {
    this->string_attribute_map[attribute_key][node] = attribute_value;
}

std::string Graph::GetStringAttribute(std::string attribute_key, int node) const {
    return this->string_attribute_map.at(attribute_key).at(node);
}

void Graph::SetDoubleAttribute(std::string attribute_key, int node, double attribute_value) {
    this->double_attribute_map[attribute_key][node] = attribute_value;
}

double Graph::GetDoubleAttribute(std::string attribute_key, int node) const {
    return this->double_attribute_map.at(attribute_key).at(node);
}

bool Graph::HasIntAttribute(std::string attribute_key, int node) const {
    return this->int_attribute_map.contains(attribute_key) && this->int_attribute_map.at(attribute_key).contains(node);
}

void Graph::AddEdge(std::pair<int, int> edge) {
    this->forward_adj_map[edge.first].insert(edge.second);
    this->backward_adj_map[edge.second].insert(edge.first);
    this->AddNode(edge.first);
    this->AddNode(edge.second);
}

int Graph::GetInDegree(int node) const {
    if (this->backward_adj_map.contains(node)) {
        return this->backward_adj_map.at(node).size();
    }
    return 0;
}

int Graph::GetOutDegree(int node) const {
    if (this->forward_adj_map.contains(node)) {
        return this->forward_adj_map.at(node).size();
    }
    return 0;
}


void Graph::AddNode(int u) {
    this->node_set.insert(u);
}

const std::set<int>& Graph::GetNodeSet() const {
    return this->node_set;
}
const std::unordered_map<int, std::set<int>>& Graph::GetForwardAdjMap() const {
    return this->forward_adj_map;
}

const std::unordered_map<int, std::set<int>>& Graph::GetBackwardAdjMap() const {
    return this->backward_adj_map;
}

void Graph::PrintGraph() const {
    for(auto const& [u,u_neighbors] : this->GetForwardAdjMap()) {
        for(const int& v : u_neighbors) {
            /* if (u < v) { */
            std::cout << u << "-" << v << std::endl;
            /* } */
        }
    }
}

void Graph::WriteGraph(std::string output_file) const {
    std::ofstream output_filehandle(output_file);
    output_filehandle << "#source,target" << std::endl;
    for(auto const& [u,u_neighbors] : this->GetForwardAdjMap()) {
        for(const int& v : u_neighbors) {
            /* if (u < v) { */
            output_filehandle << u << "," << v << std::endl;
            /* } */
        }
    }
    output_filehandle.close();
}

void Graph::WriteAttributes(std::string auxiliary_information_file) const {
    std::ofstream auxiliary_information_filehandle(auxiliary_information_file);
    auxiliary_information_filehandle << "node_id,type,year,pa_weight,rec_weight,fit_weight,fit_lag_duration,fit_peak_value,fit_peak_duration,alpha,in_degree,out_degree,assigned_out_degree,planted_nodes_line_number,generator_node_string,one_hop_neighborhood_size,two_hop_neighborhood_size\n";
    for(const auto& node_id : this->GetNodeSet()) {
        std::string node_type = this->GetStringAttribute("type", node_id);
        int year = this->GetIntAttribute("year", node_id);
        double pa_weight = -1;
        double rec_weight = -1;
        double fit_weight = -1;
        int fit_lag_duration = this->GetIntAttribute("fitness_lag_duration", node_id);
        int fit_peak_value = this->GetIntAttribute("fitness_peak_value", node_id);
        int fit_peak_duration = this->GetIntAttribute("fitness_peak_duration", node_id);
        double alpha  = -1;
        int out_degree = this->GetIntAttribute("out_degree", node_id);
        int assigned_out_degree  = -1;
        int in_degree = this->GetIntAttribute("in_degree", node_id);
        int planted_nodes_line_number = -1;
        std::string generator_node_string  = "no_generators";
        int one_hop_neighborhood_size = -1;
        int two_hop_neighborhood_size = -1;
        if(node_type == "agent") {
            pa_weight = this->GetDoubleAttribute("preferential_attachment_weight", node_id);
            rec_weight = this->GetDoubleAttribute("recency_weight", node_id);
            fit_weight = this->GetDoubleAttribute("fitness_weight", node_id);
            alpha = this->GetDoubleAttribute("alpha", node_id);
            assigned_out_degree = this->GetIntAttribute("assigned_out_degree", node_id);
            generator_node_string = this->GetStringAttribute("generator_node_string", node_id);
            if(this->HasIntAttribute("planted_nodes_line_number", node_id)) {
                planted_nodes_line_number = this->GetIntAttribute("planted_nodes_line_number", node_id);
            }
            one_hop_neighborhood_size = this->GetIntAttribute("one_hop_neighborhood_size", node_id);
            two_hop_neighborhood_size = this->GetIntAttribute("two_hop_neighborhood_size", node_id);
        }
        auxiliary_information_filehandle << node_id << "," << node_type << "," << year << "," << pa_weight << "," << rec_weight << "," << fit_weight << "," << fit_lag_duration << "," << fit_peak_value << "," << fit_peak_duration << "," << alpha << "," << in_degree << "," << out_degree << "," << assigned_out_degree << "," << planted_nodes_line_number << "," << generator_node_string << "," << one_hop_neighborhood_size << "," << two_hop_neighborhood_size << "\n";
    }
    auxiliary_information_filehandle.close();
}
