#ifndef GRAPH_H
#define GRAPH_H

#include <set>
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
        const std::map<int, std::set<int>>& GetForwardAdjMap() const;
        const std::map<int, std::set<int>>& GetBackwardAdjMap() const;
        void SetAttribute(std::string attribute_key, int node, int attribute_value);
        int GetAttribute(std::string attribute_key, int node) const;
        void ParseNodelist();
        void ParseEdgelist();
        int GetInDegree(int node) const;
        int GetOutDegree(int node) const;
        void AddNode(int u);
        void PrintGraph() const;
        void WriteGraph(std::string output_file) const;

    private:

        std::set<int> node_set;
        std::string edgelist;
        std::string nodelist;

    protected:
        std::map<int, std::set<int>> forward_adj_map;
        std::map<int, std::set<int>> backward_adj_map;
        std::map<std::string, std::map<int, int>> attribute_map;
};

#endif
