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

  DLA_params params(world.rank(), argc-2, argv+2);

  GraphPrinter graph_printer(string("DLA"), false, true, total_num_graphs);
  if (world.rank() == printer_rank) {
    cout << params;
    cout << "Producing " << total_num_graphs << " on " << world.size() << " processes" << endl;
    graph_printer.open();
  }
  for (int graph_id = 0; graph_id < graph_per_process; ++graph_id) {
    produce_graph(world, graph_printer, params);
  }
}