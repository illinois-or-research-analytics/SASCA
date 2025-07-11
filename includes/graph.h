#ifndef GRAPH_H
#define GRAPH_H

#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <queue>
#include <iostream>

class Graph {
    public:
        Graph(std::string edgelist, std::string nodelist);
        void AddEdge(std::pair<int, int> edge);

        static inline char get_delimiter(std::string filepath) {
            std::ifstream edgelist(filepath);
            std::string line;
            getline(edgelist, line);
            if (line.find(',') != std::string::npos) {
                return ',';
            } else if (line.find('\t') != std::string::npos) {
                return '\t';
            } else if (line.find(' ') != std::string::npos) {
                return ' ';
            }
            throw std::invalid_argument("Could not detect filetype for " + filepath);
        }


        const std::set<int>& GetNodeSet() const;
        const std::unordered_map<int, std::set<int>>& GetForwardAdjMap() const;
        const std::unordered_map<int, std::set<int>>& GetBackwardAdjMap() const;
        void SetIntAttribute(std::string attribute_key, int node, int attribute_value);
        int GetIntAttribute(std::string attribute_key, int node) const;
        void SetStringAttribute(std::string attribute_key, int node, std::string attribute_value);
        std::string GetStringAttribute(std::string attribute_key, int node) const;
        void SetDoubleAttribute(std::string attribute_key, int node, double attribute_value);
        double GetDoubleAttribute(std::string attribute_key, int node) const;
        bool HasIntAttribute(std::string attribute_key, int node) const;
        void ParseNodelist();
        void ParseEdgelist();
        int GetInDegree(int node) const;
        int GetOutDegree(int node) const;
        void AddNode(int u);
        void PrintGraph() const;
        void WriteGraph(std::string output_file) const;
        void WriteAttributes(std::string auxiliary_information_file) const;

    private:

        std::set<int> node_set;
        std::string edgelist;
        std::string nodelist;

    protected:
        std::unordered_map<int, std::set<int>> forward_adj_map;
        std::unordered_map<int, std::set<int>> backward_adj_map;
        std::unordered_map<std::string, std::unordered_map<int, int>> int_attribute_map;
        std::unordered_map<std::string, std::unordered_map<int, std::string>> string_attribute_map;
        std::unordered_map<std::string, std::unordered_map<int, double>> double_attribute_map;
};

#endif
