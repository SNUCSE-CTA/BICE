## Introduction

Executable files for "BICE: Exploring Compact Search Space by Using Bipartite Matching and Cell-Wide Verification".

BICE_S: BICE for subgraph search.
BICE_M_100k: BICE for subgraph matching, which finds up to 100,000 embeddings.

## Compile

The executable files were compiled on CentOS-7 using g++ (GCC) 11.2.1 20220127 (Red Hat 11.2.1-9).

## Usages
```zsh
./BICE_S -dg [data graph file] -qg [query graph file] -o [output file]
./BICE_M_100k -dg [data graph file] -qg [query graph file]
```
## Input File Format:
[graph ID]
[the number of vertices (n)]
[the label of v_1]
...
[the label of v_n]
[the number of edges (m)]
[vertex ID of e_1] [vertex ID of e_1]
...
[vertex ID of e_m] [vertex ID of e_m]

Example Input:
```zsh
#0
3
1
0
0
2
0 1
1 2
```

## Datasets

The real-world datasets and the query sets used in our paper can be downloaded at the following link.

[download](https://drive.google.com/drive/folders/1UTCWtXQyvFOx6kp8N2RB_fezuK3ulgnj?usp=sharing)

You can find files for subgraph matching query sets generation in the query_graphs_generator_subgraphmatching directory.

The graph file format for data graph and query graph is a text format to store a vertex-labeled undirected graph.

- The first line of the file should be "t #vertices #graph-ID".
- Following lines of "v vertex-ID vertex-label" indicate the vertices in the graph.
- Following lines of "e vertex-ID1 vertex-ID2 edge-label" after the vertices indicate the undirected edges in the graph. 

For example:
```zsh 
t # 0
v 0 0
v 1 1
v 2 2
v 3 3
e 0 1 0
e 0 2 0
e 1 2 0
e 1 3 0
e 2 3 0
```

As an example, the main.py is given in the query_graphs_generator_subgraphmatching directory.

## Contact
yunyoungchoi96@gmail.com