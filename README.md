## Introduction

Exectuable files for "BICE: Exploring Compact Search Space by Using Bipartite Matching and Cell-Wide Verification".

BICE_S: BICE for subgraph search.
BICE_M_100k: BICE for subgraph matching, which finds up to 100,000 embeddings.

## Compile

Executable files are compiled on CentOS-7.

## Usages
```zsh
./BICE_S -dg [data graph file] -qg [query graph file] -o [output file]
./BICE_M_100k -dg [data graph file] -qg [query graph file]
```

## Datasets

link to datasets on the paper(https://drive.google.com/drive/folders/1UTCWtXQyvFOx6kp8N2RB_fezuK3ulgnj?usp=sharing).

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
