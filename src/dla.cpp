#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include "dla.h"
#include "graph_printer.h"

using namespace std;
namespace mpi = boost::mpi;

const int default_tag = 42;
const DLA_Graph::Particle DLA_Graph::m_center{0., 0.};
default_random_engine rgen;

DLA_Graph::Particle::Particle(double _x, double _y): x(_x), y(_y) {}

DLA_params::DLA_params(int rank, int num_args, char* args[]): m_rank(rank) {
    vector<string> params;
    for (int i_arg = 0; i_arg < num_args; ++i_arg)
        params.emplace_back(args[i_arg]);
    assert(params.size() == 5);
    m_num_particles = stoi(params[0]);
    m_particle_radius = stod(params[1]);
    m_gaussian_var = stod(params[2]);
    m_spawn_radius = stod(params[3]);
    if (m_rank%2 == 0) {
        m_stick_prob = stod(params[4]);
        m_label = 2;
    } else {
        m_stick_prob = 1.;
        m_label = 1;
    }
}

int DLA_params::get_num_particles() const {
    return m_num_particles;
}

double DLA_params::get_particle_radius() const {
    return m_particle_radius;
}

double DLA_params::get_gaussian_var() const {
    return m_gaussian_var;
}

double DLA_params::get_spawn_radius() const {
    return m_spawn_radius;
}

double DLA_params::get_stick_prob() const {
    return m_stick_prob;
}

int DLA_params::get_rank() const {
    return m_rank;
}

int DLA_params::get_label() const {
    return m_label;
}

std::ostream& operator<<(std::ostream& out, DLA_params const& params) {
    out << "Rank=" << params.get_rank() << "\tstick_prob=" << params.get_stick_prob() << endl;
    return out;
}

void produce_graph(
    mpi::communicator const& com,
    GraphPrinter& graph_printer,
    DLA_params const& params) {
    DLA_Graph dla_graph(params);
    dla_graph.aggregate_particles(params.get_num_particles()-1);
    auto graph = dla_graph.get_graph();
    vector<Graph> graphs;
    mpi::gather(com, graph, graphs, printer_rank);
    if (com.rank() == printer_rank)
        print_graphs(graph_printer, graphs);
}

DLA_Graph::DLA_Graph(DLA_params const& params):
    m_particle_radius(params.get_particle_radius()),
    m_spawn_radius(params.get_spawn_radius()),
    m_stick_prob(params.get_stick_prob()),
    m_label(params.get_label()),
    m_gaussian(0., params.get_gaussian_var()),
    m_unif(0., 1.),
    m_particles{m_center} {
        Graph::Features feature{0., 0.};
        m_graph.add_node(feature);
}

Graph DLA_Graph::get_graph() const {
    return m_graph;
}

void DLA_Graph::aggregate_particles(int num_particles) {
    m_graph.set_label(m_label);
    for (int particle = 0; particle < num_particles; ++particle) {
        aggregate_particle();
    }
}

void DLA_Graph::aggregate_particle() {
    Particle particle = get_random_particle();
    while (true) {
        constrained_brownian_motion(particle);
        int neighbor_id = get_nearest_particle(particle);
        if (is_collision(m_particles[neighbor_id], particle)) {
            double coin = m_unif(rgen);
            if (coin > m_stick_prob)
                continue ;
            add_particle(neighbor_id, particle);
            break ;
        }
    }
}

void DLA_Graph::add_particle(int neighbor_id, Particle const& particle) {
    Graph::Features feature{float(particle.x), float(particle.y)};
    int particle_id = m_graph.add_node(feature);
    m_graph.add_edge(Graph::Edge(neighbor_id, particle_id));
    m_particles.push_back(particle);
}

DLA_Graph::Particle DLA_Graph::get_random_particle() const {
    double x = m_gaussian(rgen);
    double y = m_gaussian(rgen);
    Particle particle(x, y);
    double dsquare = particles_square_distance(particle, m_center);
    double d = sqrt(dsquare);
    particle.x *= (m_spawn_radius / d);
    particle.y *= (m_spawn_radius / d);
    return particle;
}

int DLA_Graph::get_nearest_particle(Particle const& particle) const {
    int closest_id = 0;
    double closest_dst = numeric_limits<double>::infinity();
    for (int i = 0; i < (int)m_particles.size(); ++i) {
        double dst = particles_square_distance(particle, m_particles[i]);
        if (dst < closest_dst) {
            closest_dst = dst;
            closest_id = i;
        }
    }
    return closest_id;
}

void DLA_Graph::constrained_brownian_motion(Particle& particle) const {
    Particle candidate;
    do {
        candidate = particle;
        pure_brownian_motion(candidate);
    } while(is_outside(candidate));
    particle = candidate;
}

bool DLA_Graph::is_farther(Particle const& candidate, Particle const& ref_pos) const {
    double dcand = particles_square_distance(candidate, m_center);
    double dref = particles_square_distance(ref_pos, m_center);
    return dcand > dref;
}

void DLA_Graph::pure_brownian_motion(Particle& particle) const {
    double dx = m_gaussian(rgen);
    double dy = m_gaussian(rgen);
    particle.x += dx;
    particle.y += dy;
}

bool DLA_Graph::is_outside(Particle const& particle) const {
    double d = particles_square_distance(m_center, particle);
    return d > m_spawn_radius * m_spawn_radius;
}

bool DLA_Graph::is_collision(Particle const& a, Particle const& b) const {
    double d = particles_square_distance(a, b);
    return d <= 4 * m_particle_radius * m_particle_radius;
}

double DLA_Graph::particles_square_distance(Particle const& a, Particle const& b) const {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return dx*dx + dy*dy;
}