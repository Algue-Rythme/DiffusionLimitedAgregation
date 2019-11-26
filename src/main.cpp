#include <boost/mpi.hpp>
#include <iostream>
#include <string>

#include "dla.h"
#include "graph_printer.h"

using namespace std;
namespace mpi = boost::mpi;

const int print_rank = 0;


int main(int argc, char *argv[])
{
  boost::mpi::environment env{argc, argv};
  boost::mpi::communicator world;
  if (world.rank() == print_rank) {
      print_graphs();
  } else {
      produce_graphs();
  }
}