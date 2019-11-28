#include <cmath>
#include <limits>
#include <random>
#include <string>
#include "dla.h"
#include "graph_printer.h"

using namespace std;
namespace mpi = boost::mpi;

const int default_tag = 0;
const DLA_Graph::Particle DLA_Graph::m_center{0., 0.};
const int seed = 577;
default_random_engine rgen(seed);

DLA_Graph::Particle::Particle(double _x, double _y): x(_x), y(_y) {}

DLA_params::DLA_params(int num_args, char* args[]) {
    vector<string> params;
    for (int i_arg = 0; i_arg < num_args; ++i_arg)
        params.emplace_back(args[i_arg]);
}

int DLA_params::get_num_particles() const {
    return 42;
}

std::ostream& operator<<(std::ostream& out, DLA_params const&) {
    return out;
}

void produce_graph(boost::mpi::communicator const& com, DLA_params const& params) {
    DLA_Graph dla_graph;
    dla_graph.aggregate_particles(params.get_num_particles());
    auto const& graph = dla_graph.get_graph();
    com.isend(printer_rank, default_tag, graph);
}

DLA_Graph::DLA_Graph():
    m_particle_radius(1.),
    m_spawn_radius(10.),
    m_gaussian(0., 1.),
    m_particles{m_center} {}

Graph const& DLA_Graph::get_graph() const {
    return m_graph;
}

void DLA_Graph::aggregate_particles(int num_particles) {
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
            add_particle(neighbor_id, particle);
        }
    }
}

void DLA_Graph::add_particle(int neighbor_id, Particle const& particle) {
    int particle_id = m_graph.add_node();
    m_graph.add_edge(Graph::Edge(neighbor_id, particle_id));
    m_particles.push_back(particle);
}

DLA_Graph::Particle DLA_Graph::get_random_particle() const {
    double x = m_gaussian(rgen);
    double y = m_gaussian(rgen);
    Particle particle(x, y);
    double dsquare = particles_square_distance(particle, m_center);
    double d = sqrt(dsquare);
    double norm = m_spawn_radius / d;
    particle.x *= norm;
    particle.y *= norm;
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
    Particle candidate = particle;
    do {
        pure_brownian_motion(candidate);
    } while(is_outside(candidate));
    particle = candidate;
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
    return d <= m_particle_radius * m_particle_radius;
}

double DLA_Graph::particles_square_distance(Particle const& a, Particle const& b) const {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return dx*dx + dy*dy;
}