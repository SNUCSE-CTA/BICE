#include "randomwalk.h"

// load data graph from file
RandomWalk::RandomWalk(char* a_data_graph_file) {
  ifstream infile(a_data_graph_file);
  int datagraph_num_vertices;
  char tag;
  int id;
  int label;

  string line;
  bool first_time = true;
  while (getline(infile, line)) {
    istringstream iss(line);
    iss >> tag;
    switch (tag) {
      case 't':
        break;
      case 'v':
        iss >> id;
        iss >> label;
        datagraph_vertices.push_back(label);
        break;
      case 'e':
        if (first_time) {
          first_time = false;
          datagraph_num_vertices = datagraph_vertices.size();
          for (int i = 0; i < datagraph_num_vertices; i++) {
            set<int> s;
            datagraph_edges.push_back(s);
          }
        }
        int x, y;
        iss >> x >> y;
        datagraph_edges[x].insert(y);
        datagraph_edges[y].insert(x);
        break;
      default:
        break;
    }
  }
  cout << "|V|: " << datagraph_vertices.size() << endl;
}

void RandomWalk::PrintDataGraph() {
  cout << "t 0" << datagraph_vertices.size() << endl;
  for (int i = 0; i < datagraph_vertices.size(); i++) {
    cout << "v " << i << " " << datagraph_vertices[i] << endl;
  }

  for (int i = 0; i < datagraph_edges.size(); i++) {
    for (set<int>::iterator j = datagraph_edges[i].begin();
         j != datagraph_edges[i].end(); j++) {
      if (*j > i)  // caution: we assume there is no loop (i.e. *j==i)
        cout << "e " << i << " " << *j << " 0" << endl;
    }
  }
}

// select v_0 with uniform random probability
int RandomWalk::Initialize() {
  mWalkCount = 1;

  mCurrentVertex = Random(0, datagraph_vertices.size() - 1);
  return mCurrentVertex;
}

// select v_t by selecting one of v_{t-1}'s neighbor with uniform probability
int RandomWalk::Walk() {
  if (datagraph_edges[mCurrentVertex].size() == 0) {
    mWalkCount = -1;
    return -1;
  } else
    mWalkCount++;

  int neighborIndex = Random(0, datagraph_edges[mCurrentVertex].size() - 1);
  set<int>::iterator neighbor = datagraph_edges[mCurrentVertex].begin();
  std::advance(neighbor, neighborIndex);
  mCurrentVertex = *neighbor;
  return mCurrentVertex;
}

// returns a random number among [a_left, a_right]
int RandomWalk::Random(int a_left, int a_right) {
  return (rand() % (a_right - a_left + 1) + a_left);
}

int RandomWalk::Label(int a_vertex_id) {
  return datagraph_vertices[a_vertex_id];
}

int RandomWalk::WalkCount() { return mWalkCount; }

set<int> RandomWalk::AdjacentSet(int a_vertex_id) {
  return datagraph_edges[a_vertex_id];
}
