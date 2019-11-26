#ifndef GRAPH_PRINTER_H
#define GRAPH_PRINTER_H

#include <fstream>
#include <string>
#include <vector>

class Graph {
    public:
    Graph() = default;

    struct Edge {
        Edge(int, int);

        int start;
        int end;
    };

    typedef std::vector<float> Features;
    typedef std::vector<int> Labels;

    private:
    std::vector<Edge> edges;
    std::vector<Features> features;
    std::vector<Labels> labels;
};

class GraphPrinter {
    public:
    GraphPrinter(std::string const&, bool, bool);
    GraphPrinter(GraphPrinter const&) = delete;
    GraphPrinter& operator=(GraphPrinter const&) = delete;

    GraphPrinter& operator<<(Graph const&);

    void close();

    private:
    void create_directory();
    void print_indicator();
    void print_adj();
    void print_node_labels();
    void print_node_features();

    std::string m_name;
    std::ofstream m_indicator_file;
    std::ofstream m_adj_file;
    std::ofstream m_node_labels_file;
    std::ofstream m_node_features_file;
    int m_offset;
};


void print_graphs();

#endif