#include <string>
#include <boost/filesystem.hpp>

#include "graph_printer.h"

using namespace std;
namespace mpi = boost::mpi;
namespace filesystem = boost::filesystem;

const int printer_rank = 0;

Graph::Edge::Edge(int _start, int _end): start(_start), end(_end) {}

void Graph::add_edge(Edge const& edge) {
    edges.push_back(edge);
}

int Graph::add_node() {
    num_nodes += 1;
    return num_nodes-1;
}

GraphPrinter::GraphPrinter(std::string const& name, bool has_node_features, bool has_node_labels)
: m_name(name), m_offset(1), m_graph_id(1) {
    auto path = create_directory();
    m_indicator_file.open(path.string() + "_graph_indicator.txt");
    m_adj_file.open(path.string() + "_graph_A.txt");
    if (has_node_labels)
        m_node_labels_file.open(path.string() + "_node_labels.txt");
    if (has_node_features)
        m_node_features_file.open(path.string() + "_node_features.txt");
}

GraphPrinter& GraphPrinter::operator<<(Graph const& graph) {
    print_indicator(graph);
    print_adj(graph);
    print_node_labels(graph);
    print_node_features(graph);
    m_offset += graph.num_nodes;
    m_graph_id += 1;
    return *this;
}

void GraphPrinter::print_indicator(Graph const& graph) {
    for (int node = 0; node < graph.num_nodes; ++node)
        m_indicator_file << m_graph_id << '\n';
}

void GraphPrinter::print_adj(Graph const& graph) {
    for (const auto& edge : graph.edges) {
        int a = edge.start + m_offset;
        int b = edge.end + m_offset;
        m_adj_file << a << ", " << b << '\n';
        m_adj_file << b << ", " << a << '\n';
    }
}

void GraphPrinter::print_node_labels(Graph const& graph) {
    if (!m_node_labels_file.is_open())
        return ;
    for (auto const& node_labels : graph.labels) {
        print_vector(m_node_labels_file, node_labels);
    }
}

void GraphPrinter::print_node_features(Graph const& graph) {
    if (!m_node_features_file.is_open())
        return ;
    for (auto const& node_features : graph.features) {
        print_vector(m_node_features_file, node_features);
    }
}

void GraphPrinter::close() {
    m_indicator_file.close();
    m_adj_file.close();
    if (m_node_labels_file.is_open())
        m_node_labels_file.close();
    if (m_node_features_file.is_open())
        m_node_features_file.close();
}

int GraphPrinter::num_graphs() const {
    return m_graph_id - 1;
}

filesystem::path GraphPrinter::create_directory() {
    filesystem::path dir(m_name + "/");
    filesystem::create_directory(dir);
    return dir / m_name;
}

void print_graphs(mpi::communicator const& com, GraphPrinter& graph_printer) {
    Graph graph;
    mpi::request request = com.irecv(mpi::any_source, mpi::any_tag, graph);
    while (request.test()) {
        graph_printer << graph;
        request = com.irecv(mpi::any_source, mpi::any_tag, graph);
    }
}