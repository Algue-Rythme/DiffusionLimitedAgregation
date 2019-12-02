#ifndef GRAPH_PRINTER_H
#define GRAPH_PRINTER_H

#include <fstream>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/progress.hpp>

extern const int printer_rank;

template<typename T>
void print_vector(std::ofstream& file, std::vector<T> const& v) {
    for (int id = 0; id < (int)v.size(); ++id) {
        file << v[id];
        if (id+1 == (int)v.size())
            file << "\n";
        else
            file << ", ";
    }
}

class Graph {
    public:
    Graph();

    struct Edge {
        Edge() = default;
        Edge(int, int);

        int start;
        int end;
    };

    typedef std::vector<float> Features;
    typedef std::vector<int> Labels;
    
    void add_edge(Edge const&);
    int add_node();
    int add_node(Features const&);

    void set_label(int);

    int num_nodes;
    std::vector<Edge> edges;
    std::vector<Features> features;
    std::vector<Labels> labels;
    int label;
};

template<typename Archive>
void serialize(Archive& ar, Graph::Edge& edge, unsigned int version) {
    ar & edge.start;
    ar & edge.end;
}

template<typename Archive>
void serialize(Archive& ar, Graph& graph, unsigned int version) {
    ar & graph.num_nodes;
    ar & graph.edges;
    ar & graph.features;
    ar & graph.labels;
    ar & graph.label;
}

class GraphPrinter {
    public:
    GraphPrinter(std::string const&, bool, bool, int);
    GraphPrinter(GraphPrinter const&) = delete;
    GraphPrinter& operator=(GraphPrinter const&) = delete;

    void open();

    GraphPrinter& operator<<(Graph const&);

    void close();
    int num_graphs() const;

    private:
    boost::filesystem::path create_directory();
    void print_indicator(Graph const& graph);
    void print_adj(Graph const& graph);
    void print_node_labels(Graph const& graph);
    void print_node_features(Graph const& graph);

    std::string m_name;
    std::ofstream m_indicator_file;
    std::ofstream m_adj_file;
    std::ofstream m_graph_labels;
    std::ofstream m_node_labels_file;
    std::ofstream m_node_features_file;
    int m_offset;
    int m_graph_id;
    bool m_has_node_labels;
    bool m_has_node_features;
    int m_total_num_graphs;
    std::unique_ptr<boost::progress_display> m_show_progress;
};

void print_graphs(GraphPrinter& graph_printer, std::vector<Graph> const& graphs);

#endif