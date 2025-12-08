# Shortest path, Minimum Spanning Tree in a graph

The Dijkstra algorithm for the shortest path in a graph, and the prim's algotithm for the Minimum Spanning Tree (MST) are implemented. These implementations use the Min-Heap structure for updating the set of new vertices and extracting the minimum in an efficient way.

## Data

The graph is expected to be provided as the number of vertices, followed by 1 row per edge. The precise format is as follows: 

```
#number_of_vertices
#extremity1 #extremity2 #weight
#extremity1 #extremity2 #weight
...
```

An example is provided as <code>SampleTestData.txt</code> .

## How to use

Instructions for compiling and running in debug mode using the compiler <code>g++</code> are listed in <code>Makefile</code>. They can be called using:
* Compile:
```bash
make build
```
* Run:
```bash
make run
```
* Remove the executable:
```bash
make clean
```

## Feedback

All <b>suggestions</b> should be adressed to its author Alexandre Aksenov:
* GitHub: Alexandre-aksenov
* Email: alexander1aksenov@gmail.com