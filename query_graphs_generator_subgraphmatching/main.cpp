#include <algorithm>
#include <fstream>

#include "randomwalk.h"

using namespace std;

double fRand(double fMin, double fMax);
void PrintQuery(int, RandomWalk*, vector<int>, vector<set<int> >, ofstream&);

// usage ./a.out datagraph_filename #querygraphs #queryvertex
// sparse_query_filename non-sparse_query_file_name
int main(int argc, char* argv[]) {
  srand(time(NULL));

  RandomWalk walker(argv[1]);
  // cerr << "walker instance ready" << '\n';
  int queryCount = atoi(argv[2]);  // the number of query graphss to generate
  int queryVertexCount =
      atoi(argv[3]);                  // the number of vertices in a query graph
  ofstream sparseQueryFile(argv[4]);  // sparse query file name
  ofstream nonsparseQueryFile(argv[5]);  // non sparse query file name

  for (int i = 0; i < queryCount; i++) {
    int currNumEdge = 0;

    vector<int> vertices;
    vector<set<int> > edges;
    int prevVertex = walker.Initialize();
    // cerr << "walker.Initialize() done" << '\n';
    vertices.push_back(prevVertex);
    set<int> edgeSet;
    edges.push_back(edgeSet);

    while (vertices.size() < queryVertexCount) {
      int currVertex = walker.Walk();
      if (currVertex == -1) break;

      if (std::find(vertices.begin(), vertices.end(), currVertex) ==
          vertices.end()) {
        vertices.push_back(currVertex);
        set<int> s;
        edges.push_back(s);
      }

      int prevIndex = -1;
      int currIndex = -1;
      for (int j = 0; j < vertices.size(); j++) {
        if (vertices[j] == prevVertex) prevIndex = j;
        if (vertices[j] == currVertex) currIndex = j;
      }

      if (prevIndex != currIndex) {
        if (find(edges[prevIndex].begin(), edges[prevIndex].end(), currIndex) ==
            edges[prevIndex].end()) {
          edges[prevIndex].insert(currIndex);
          edges[currIndex].insert(prevIndex);

          currNumEdge += 2;
        }
      }
      prevVertex = currVertex;

      if (walker.WalkCount() == -1 || walker.WalkCount() > 1000000) {
        break;
      }
    }
    if (walker.WalkCount() == -1 || walker.WalkCount() > 1000000) {
      i--;
      continue;
    }
    double avgDegreeMin = currNumEdge / (double)queryVertexCount;
    vector<int> backup_vertices;
    backup_vertices.assign(vertices.begin(), vertices.end());
    vector<set<int> > backup_edges;
    backup_edges.assign(edges.begin(), edges.end());

    int backupCurrNumEdge = currNumEdge;
    vector<pair<int, int> > edge_pool;

    for (int j = 0; j < vertices.size(); j++) {
      set<int> orin_edge = walker.AdjacentSet(vertices[j]);
      set<int>::iterator orin_it = orin_edge.begin();
      while (orin_it != orin_edge.end()) {
        int index = -1;
        for (int k = 0; k < vertices.size(); k++) {
          if (vertices[k] == *orin_it) {
            index = k;
            break;
          }
        }
        if (index != -1) {
          if (std::find(edges[j].begin(), edges[j].end(), index) ==
              edges[j].end()) {
            pair<int, int> p(j, index);
            edge_pool.push_back(p);
          }
        }
        orin_it++;
      }
    }
    vector<pair<int, int> > backup_edge_pool;
    backup_edge_pool.assign(edge_pool.begin(), edge_pool.end());

    while (edge_pool.size() > 0) {
      int x = edge_pool[0].first;
      int y = edge_pool[0].second;
      if (x != y and
          std::find(edges[x].begin(), edges[x].end(), y) == edges[x].end()) {
        edges[x].insert(y);
        edges[y].insert(x);
        currNumEdge += 2;
      }
      edge_pool.erase(edge_pool.begin());
    }

    double avgDegreeMax = currNumEdge / (double)queryVertexCount;
    edge_pool.clear();
    edge_pool.assign(backup_edge_pool.begin(), backup_edge_pool.end());
    if (avgDegreeMax <= 3 || avgDegreeMin > 3) {
      --i;
      continue;
    }

    double sparseNumEdge = (int)(fRand(avgDegreeMin, 3) * queryVertexCount);
    double nonsparseNumEdge = (int)(fRand(3, avgDegreeMax) * queryVertexCount);

    vertices.clear();
    vertices.assign(backup_vertices.begin(), backup_vertices.end());
    edges.clear();
    edges.assign(backup_edges.begin(), backup_edges.end());

    currNumEdge = backupCurrNumEdge;
    while (edge_pool.size() > 0 && sparseNumEdge > currNumEdge) {
      int rand_index = walker.Random(0, edge_pool.size() - 1);
      int x = edge_pool[rand_index].first;
      int y = edge_pool[rand_index].second;
      if (x != y and
          std::find(edges[x].begin(), edges[x].end(), y) == edges[x].end()) {
        edges[x].insert(y);
        edges[y].insert(x);
        currNumEdge += 2;
      }
      edge_pool.erase(edge_pool.begin() + rand_index);
    }
    PrintQuery(i, &walker, vertices, edges, sparseQueryFile);
    while (edge_pool.size() > 0 and nonsparseNumEdge >= currNumEdge) {
      int rand_index = walker.Random(0, edge_pool.size() - 1);
      int x = edge_pool[rand_index].first;
      int y = edge_pool[rand_index].second;
      if (x != y and
          std::find(edges[x].begin(), edges[x].end(), y) == edges[x].end()) {
        edges[x].insert(y);
        edges[y].insert(x);
        currNumEdge += 2;
      }
      edge_pool.erase(edge_pool.begin() + rand_index);
    }
    PrintQuery(i, &walker, vertices, edges, nonsparseQueryFile);
  }

  sparseQueryFile.close();
  nonsparseQueryFile.close();
  return 0;
}

void PrintQuery(int a_graph_index, RandomWalk* a_walker, vector<int> a_vertices,
                vector<set<int> > a_edges, ofstream& a_file) {
  a_file << "t # " << a_graph_index << endl;
  for (int j = 0; j < a_vertices.size(); j++) {
    a_file << "v " << j << " " << a_walker->Label(a_vertices[j]) << endl;
  }
  for (int k = 0; k < a_edges.size(); k++) {
    for (set<int>::iterator j = a_edges[k].begin(); j != a_edges[k].end();
         j++) {
      if (*j > k)  // caution: we assume there is no loop (i.e. *j==k)
        a_file << "e " << k << " " << *j << " 0" << endl;
    }
  }
}

double fRand(double fMin, double fMax) {
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}
