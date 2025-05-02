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
        this->SetAttribute("year", integer_node, integer_year);
    }
}

void Graph::SetAttribute(std::string attribute_key, int node, int attribute_value) {
    this->attribute_map[attribute_key][node] = attribute_value;
}

int Graph::GetAttribute(std::string attribute_key, int node) const {
    return this->attribute_map.at(attribute_key).at(node);
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
const std::map<int, std::set<int>>& Graph::GetForwardAdjMap() const {
    return this->forward_adj_map;
}

const std::map<int, std::set<int>>& Graph::GetBackwardAdjMap() const {
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
