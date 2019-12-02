# DiffusionLimitedAgregation
Produce graphs with DLA model.  
Efficient parallel implementation using Boost::MPI.

## Usage

```
cmake .
./launch_n.sh 6 1000 100 1. 2. 40. 0.05
```

The cmake command is requireed only once.
This command will launch the program with 6 threads.
It will produce 1000 graphs, with 100 nodes per graph.
Each node is of diameter 1. The Brownian motion has a standard deviation equal to 2.
The radius of the initialization sphere is 40
The probability that a new particle stick is 0.05