#ifndef DLA_H
#define DLA_H

#include <iostream>
#include <random>
#include <vector>
#include <boost/mpi.hpp>
#include "graph_printer.h"

class DLA_params {
    public:
    DLA_params(int, char* []);

    int get_num_particles() const;
    double get_particle_radius() const;
    double get_gaussian_var() const;
    double get_spawn_radius() const;

    private:

    int m_num_particles;
    double m_particle_radius;
    double m_gaussian_var;
    double m_spawn_radius;
};

std::ostream& operator<<(std::ostream&, DLA_params const&);

class DLA_Graph {
    public:

    DLA_Graph(DLA_params const&);

    void aggregate_particles(int);
    Graph get_graph() const;

    struct Particle {
        Particle() = default;
        Particle(double, double);
        double x;
        double y;
    };

    private:

    static const Particle m_center;

    void aggregate_particle();
    Particle get_random_particle() const;
    void constrained_brownian_motion(Particle&) const;
    void pure_brownian_motion(Particle&) const;
    bool is_outside(Particle const&) const;
    bool is_farther(Particle const&, Particle const&) const;
    int get_nearest_particle(Particle const&) const;
    void add_particle(int, Particle const& particle);
    bool is_collision(Particle const&, Particle const&) const;
    double particles_square_distance(Particle const& a, Particle const& b) const;

    double m_particle_radius;
    double m_spawn_radius;
    mutable std::normal_distribution<double> m_gaussian;
    Graph m_graph;
    std::vector<Particle> m_particles;
};

void produce_graph(boost::mpi::communicator const&, GraphPrinter&, DLA_params const& params);

#endif