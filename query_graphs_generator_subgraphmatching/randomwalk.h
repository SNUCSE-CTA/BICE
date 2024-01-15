#pragma once

#include <math.h>
#include <stdlib.h>  //srand, rand, atoi
#include <time.h>    //time

#include <algorithm>  //sort
#include <chrono>
#include <fstream>  //ofstream
#include <iostream>
#include <queue>    //queue
#include <set>      //set
#include <sstream>  //istringstream
#include <string>   //string
#include <utility>  //pair, make_pair
#include <vector>   //vector

using namespace std;

class RandomWalk {
 private:
  int mCurrentVertex = -1;  // v_t
  int mWalkCount = -1;

  vector<int> datagraph_vertices;  // datagraph_vertices[i] contains label of
                                   // v_i
  vector<set<int> >
      datagraph_edges;  // datagraph_edges[i] contains adjacent vertices of v_i
 public:
  RandomWalk(char*);      // load data graph from file
  void PrintDataGraph();  // print in igraph format
  int Initialize();       // pick v_0 and return v0
  int Walk();             // pick v_t from v_{t-1} and return v_t
  int Random(int, int);   // return uniform random number between [arg1,arg2]
  int Label(int);         // return label of datagraph vertex
  int WalkCount();
  set<int> AdjacentSet(int);
};
