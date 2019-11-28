#include <cmath>
#include <iostream>
#include <string>
#include <boost/mpi.hpp>

#include "dla.h"
#include "graph_printer.h"

using namespace std;

int main(int argc, char *argv[])
{
  boost::mpi::environment env{argc, argv};
  boost::mpi::communicator world;

  int total_num_graphs = stoi(string(argv[1]));
  int graph_per_process = int(ceil(double(total_num_graphs) / double(world.size())));
  total_num_graphs = graph_per_process * world.size();
  cout << "Producing " << total_num_graphs << " on " << world.size() << " processes" << endl;

  DLA_params params(argc-2, argv+2);
  cout << params;

  if (world.rank() == printer_rank) {
    GraphPrinter graph_printer("DLA", false, false);
    for (int graph_id = 0; ; ++graph_id) {
      if (graph_id < graph_per_process)
        produce_graph(world, params);
      print_graphs(world, graph_printer);
      if (graph_printer.num_graphs() >= total_num_graphs)
        break ;
    }
  } else {
    for (int graph_id = 0; graph_id < graph_per_process; ++graph_id) {
        produce_graph(world, params);
    }
  }
}