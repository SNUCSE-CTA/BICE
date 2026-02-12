#ifndef DAF
#define DAF

#define FILTERING_BY_NEIGHBOR_SAFETY 
#define LEAF_ADAPTIVE_MATCHING 
#define N_OPTIMIZATION 
#define PRUNING_BY_EQUIVALENCE_SETS

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <iomanip>
#include <cfloat>
#include <climits>
#include <set>
#include <map>
#include <queue>
#include <limits>
#include <cassert>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <unistd.h>
#include <limits.h> //for LLONG_MAX which is 9223372036854775807 
#include <unordered_map>
#include <string>
#include <cstring> 
#ifdef MAPPING_FUNCTION_LOG
#include <bitset> 
#endif
#include "process.h"
using namespace std;

#ifdef PRUNING_BY_EQUIVALENCE_SETS
class Array
{
public: 
  int* b;
  int size;
  Array(){};
  Array(int* arr, int* pos, int start, int end){
    size = end - start;
    b = new int[size];
    for(int i = 0; i < size; ++i)
      b[i] = arr[pos[start + i]];
    sort(b, b + size);
  }

  bool operator==(const Array& other) const{
    if(size != other.size) 
      return false;
    else
    {
      bool ret = true;
      for(int x = 0; x < size; ++x)
	ret &= (b[x] == other.b[x]);
      return ret;
    }
  }
};

struct ArrayHash{
  size_t operator()(const Array& array) const{
    size_t ret = 0;
    for(int i = 0; i < array.size; ++i)
      ret = ret ^ array.b[i];
    return ret + array.b[0];
  }
};

unordered_map<Array, int, ArrayHash> cellToID;
#ifdef SUBGRAPH_MATCHING
vector<int> cellPos;
vector<int> cellList;
#else
int* cellPos;
int* cellList;
#endif
unordered_map<Array, int, ArrayHash>::iterator uIter;
unordered_map<int, Array>::iterator iIter;
int uI;
int nCandOffset;
int globalCellID = 1;
int* posArray;
int* candArray;
int* candOffset;
#endif

// Main result variables
bool* answer;
int nCandidate = 0;
double nMatch, nCurrMatch, nRemainingMatch;
#ifdef SUBGRAPH_MATCHING
const double nMaxMatch = 100000;
#else
const double nMaxMatch = 1;
#endif
long long int recursiveCallCount = 0;
long long int recursiveCallCountPerGraph;
// Elapsed time variable used in the main function
bool exceedTimeLimit;
chrono::high_resolution_clock::time_point startPoint, endPoint;
long long int quotient, lastQuotient;
double elapsedTime;
const double timeLimit = 600000;
double filteringTime;
double verificationTime;
double totalTime;
double sumFilteringTime = 0;
double sumVerificationTime = 0;
double sumTotalTime = 0;
int auxDataStructureSize; 
// File paths
string dataGraphFile;
string queryGraphFile;
string answerFile;

inline void swapValue(int & v1, int & v2) 
{
  int temp = v1;
  v1 = v2;
  v2 = temp;
}

inline void swapValue(pair<int, int> & v1, pair<int, int> & v2) 
{
  pair<int, int> temp = v1;
  v1 = v2;
  v2 = temp;
}
// Variables for buildling CS
struct CandidateSpace 
{
  int size; // the size for both path and candidates
  int* candidates; //candidate set
  int*** adjacent = NULL; //adjacent[i][j] = candidates of this unit when the i-th parent, regarding DAG, mapped to the j-th candidate of the parent.
  int** nAdjacent = NULL; //nAdjacent[i][j] = size of back_trak_index[i][j]. That is, the number of candidates of this unit when the i-th parent mapped to the j-th candidate of the parent.
  
  #ifdef N_OPTIMIZATION
  int** capacity = NULL;
  #else
  int* nParentCand = NULL; //nParentCand[i] = the number of candidates of the i-th parent.
  #endif
  long long* weight = NULL;
  #ifdef PRUNING_BY_EQUIVALENCE_SETS
  int* cell;
  int* cellVertex;
  int nCellVertex;
  #endif
  int* eqSetId;
  int szEqSetId;
  
  vector < pair < int, int > > intervals_cell;
  int * components_cell = nullptr;
};
int uSequence[MAX_NUM_VERTEX];
int uSequenceSize = 0;

int orderDAG[MAX_NUM_VERTEX];

CandidateSpace candSpace[MAX_NUM_VERTEX];
int * arrayToClean;
int toCleanIndex = 0;
char * isVisited; //data count
int * candToPos; //data count

// Variables for maximal matching
int nodeIndex = 0;
int dist[MAX_NUM_VERTEX + 1];
int * uCand;
int * necMapping;
int * pairU;
int * pairV;
number_type leafCandSize = 0;
int leafCandIndex = 0;
vector<int> leafCand;
pair<int, int> * leafCandInfo;
int * * sumNECCand; //added for the optimization of performing all MM first
int * nSumNECCand;   //added for the optimization of performing all MM first
int * vCandIndex; //length: maxNumDataVertex. its each entry stores the number of query nodes contain the corresponding candidate.
int posSumNECCand;
vector < pair <int, int> > * vCandInfo; // a combination of v_cands and v_cands_pos
pair <int, int>  * uCandInfo; //two-d array:
//[i,j] means the j-th candidate pair of the i-th nec node(in the nec-node categorized by label)
//each candidate pair stores a candidate and the position of "the i-th nec node" in this candidate's u_cand set.
int * flagSumNECCand; //indicate whether or not a candidate has been added into the candidate already.

// Variables for search
int uCurr;
int labelID = -1;
int NECSetSize = -1;
pair<int, pair<int, double> >  NECRank[MAX_NUM_VERTEX];
int NECRegionSize[MAX_NUM_VERTEX];
int NECCountSet[MAX_NUM_VERTEX];
int NECCountSetSize;
pair <int, int>  vNECCount[MAX_NUM_VERTEX]; // used in the map special tree function
double PRE_COMPUTED_PERMUTATION;
int candPos[MAX_NUM_VERTEX];
int currMapping[MAX_NUM_VERTEX];
int nMappedParent[MAX_NUM_VERTEX];
int MappedParent[MAX_NUM_VERTEX][MAX_QUERY_DEGREE];
#ifdef QUERY_VERTEX_WITH_MAX
weight_type WEIGHT_MIN = LLONG_MIN;
#else
weight_type WEIGHT_MAX = LLONG_MAX;
#endif
int posExtendable = 0;
int minPosition = -1;
long long optWeight = LLONG_MAX;
int* iec[MAX_NUM_VERTEX][MAX_QUERY_DEGREE];
int iecSize[MAX_NUM_VERTEX][MAX_QUERY_DEGREE];
weight_type weightArray[MAX_NUM_VERTEX];
weight_type currWeight[MAX_NUM_VERTEX];
int position[MAX_NUM_VERTEX];
int nInserted[MAX_NUM_VERTEX];
int numMatchedChild[MAX_NUM_VERTEX];
struct Stack 
{
  int* address = NULL;
  int addressSize;
  int addressPos;
  int vertex;
  void update(
    int* add, int addSize,
    int addPos, int ver
  )
  {
    this->address = add;
    this->addressSize = addSize;
    this->addressPos = addPos;
    this->vertex = ver;
  }
  //for failing set
  bool propagated = false;
  uint64_t* failingSet;
  #ifdef PRUNING_BY_EQUIVALENCE_SETS
  uint64_t* buffer;
  uint64_t* conflictCell;
  bool* isPruned;
  int nPruned;
  int* pruned;
  #ifdef SUBGRAPH_MATCHING
  uint64_t* auxBuffer;
  uint64_t* nonConflictCell;
  bool* neverVisited;
  double* nFoundMatch;
  #endif
  #endif
  #ifdef LEAF_ADAPTIVE_MATCHING
  int* problemChild;
  int nProblemChild;
  #endif
  #ifdef CANDIDATE_VERTEX_WITH_MAX
  bool* isVisited;
  #endif

  long long* visEqSetId;
  long long visEqSetIdCnt = 0;
};
Stack element[MAX_NUM_VERTEX];
Stack *currE;

// Variables for pruning by failing sets
int FAILING_SET_SIZE = ((MAX_NUM_VERTEX + 63) >> 6);
uint64_t* ancestors[MAX_NUM_VERTEX];
bool isRedundant;
// Variables for pruning by equivalence sets
#ifdef PRUNING_BY_EQUIVALENCE_SETS
int CELL_SIZE;
int nMaxCell;
#endif

// Variables for filtering by neighbor safety
#ifdef FILTERING_BY_NEIGHBOR_SAFETY
using lint = long long;  
bool on = false;
int index_vis_color = 0;
int index_set_color = 0;
int index_vis_adj = 0;
int index_ngb_existence = 0;
double s1 = 0;
double s2 = 0;
#ifdef HUGE_GRAPH
int DATAV;
int LABEL;
int **color;//label 2  
int *set_color;//label
int *vis_adj;//datav
int **num_ngb_existence;//datav 2
lint *cnt_included_cs;//datav
#else
const int DATAV = 4300000;
const int LABEL = 50000;
int color[LABEL][2] = { 0, };
int set_color[LABEL] = { 0, };
int vis_adj[DATAV] = { 0, };
int num_ngb_existence[DATAV][2] = { 0, }; 
lint cnt_included_cs[DATAV] = { 0, };
#endif
#endif

// Temporary variables
uint64_t* b;
uint64_t arr;

//======================= sorting functions ========================================
inline bool sortLLCNode(pair<int, pair<int, double> > a, pair<int, pair<int, double> > b) {// the second pair is <node sum, cand sum>
  if (a.second.first == b.second.first)
    return a.second.second < b.second.second;
  else
    return a.second.first < b.second.first;
}
inline bool sortByEffectiveCandSize(pair<int, int> p1, pair<int, int> p2) {
  return NECRegionSize[p1.first] < NECRegionSize[p2.first];
}
inline bool sortBySecondElement(pair<int, int> p1, pair<int, int> p2) {
  return (p1.second < p2.second);
}
inline bool sortQueryByDegree(const int id1, const int id2) {
  return (currQ->degree[id1] > currQ->degree[id2]);
}
inline bool sortQueryByLabelFrequency(const int id1, const int id2) {
  return (currG->labelFrequency[currQ->label[id1]] < currG->labelFrequency[currQ->label[id2]]);
}
inline bool sortPairsByWeight(const pair<long long, long long>& p1, const pair<long long, long long>& p2){
  return p1.second < p2.second;
}

#ifdef SUBGRAPH_MATCHING
inline string doubleToString(double value)
{
  string strVal = to_string(value);
  int pos = strVal.find_first_of(".");
  return strVal.substr(0, pos);
}
#endif

inline double factorization(int x) 
{
  switch (x) {
  case 1:
    return 1;
  case 2:
    return 2;
  case 3:
    return 6;
  case 4:
    return 24;
  case 5:
    return 120;
  case 6:
    return 720;
  case 7:
    return 5040;
  case 8:
    return 40320;
  case 9:
    return 362880;
  case 10:
    return 3628800;
  case 11:
    return 39916800;
  case 12:
    return 479001600;
  case 13:
    return 6227020800;
  case 14:
    return 87178291200;
  case 15:
    return 1307674368000;
  case 16:
    return 20922789888000;
  case 17:
    return 355687428096000;
  case 18:
    return 6402373705728000;
  case 19:
    return 121645100408832000;
  case 20:
    return 2432902008176640000;
  }

  double result = 2432902008176640000;
  for (double i = 21; i <= x; i++)
    result *= i;

  return result;
}

inline void AllocateForDataGraph()
{
  #ifdef FILTERING_BY_NEIGHBOR_SAFETY
  #ifdef HUGE_GRAPH
  DATAV = maxNumDataVertex + 10;
  LABEL = nUniqueLabel + 10;
  color = new int*[LABEL];
  for(int i=0;i<LABEL;i++){
    color[i] = new int[2];
    for(int j=0;j<2;j++){
      color[i][j] = 0;
    }
  }
  set_color = new int[LABEL];
  for(int i=0;i<LABEL;i++){
    set_color[i] = 0;
  }
  vis_adj = new int[DATAV];
  for(int i=0;i<DATAV;i++){
    vis_adj[i] = 0;
  }
  num_ngb_existence = new int*[DATAV];
  for(int i=0;i<DATAV;i++){
    num_ngb_existence[i] = new int[2];
    for(int j=0;j<2;j++){
      num_ngb_existence[i][j] = 0;
    }
  }
  cnt_included_cs = new lint[DATAV];
  for(int i=0;i<DATAV;i++){
    cnt_included_cs[i] = (lint)(0);
  }
  #endif
  #endif
  pairV = new int[maxNumDataVertex];
  
  vCandIndex = new int[maxNumDataVertex];
  memset(vCandIndex, 0, sizeof(int) * maxNumDataVertex);
  
  flagSumNECCand = new int[maxNumDataVertex];
  memset(flagSumNECCand, -1, sizeof(int) * maxNumDataVertex);

  vCandInfo = new vector< pair<int, int> >[maxNumDataVertex];
  
  sumNECCand = new int *[nUniqueLabel + 1];
  for (int i = 0; i < (nUniqueLabel + 1); i++)
    sumNECCand[i] = new int[maxNumCandidate];
  nSumNECCand = new int[nUniqueLabel + 1]();
  
  mapTo = new int[maxNumDataVertex];
  memset(mapTo, -1, sizeof(int) * maxNumDataVertex);

  arrayToClean = new int[maxNumDataVertex];
  isVisited = new char[maxNumDataVertex];
  memset(isVisited, 0, sizeof(char) * maxNumDataVertex);
  candToPos = new int[maxNumDataVertex];
  memset(candToPos, -1, sizeof(int) * maxNumDataVertex);
  
  answer = new bool[nGraph];
  memset(answer, false, sizeof(bool) * nGraph);
  #ifdef PRUNING_BY_EQUIVALENCE_SETS
  #ifdef SUBGRAPH_MATCHING
  cellPos.reserve(maxNumCandidate);
  cellList.reserve(maxNumDataVertex);
  #else
  cellPos = new int[maxNumCandidate];
  cellList = new int[maxNumDataVertex];
  #endif
  #endif
}

inline void Deallocate(Graph& query)
{
  delete[] query.NECMapping;
  delete[] query.NECElement;
  delete[] query.NECMap;
  #ifdef LEAF_ADAPTIVE_MATCHING
  delete[] query.isProblemLeaf;
  #endif
  leafCand.clear();
  delete[] leafCandInfo;
  delete[] necMapping;
  delete[] uCand;
  delete[] pairU;
  delete[] uCandInfo;
  for(int i = 0; i < maxNumDataVertex; ++i)
    vCandInfo[i].clear();
 
  #ifdef PRUNING_BY_EQUIVALENCE_SETS
  #ifdef SUBGRAPH_MATCHING
  delete[] nActive;
  for(int i = 0; i < maxNumDataVertex; ++i)
    delete[] active[i];
  delete[] active;
  #endif
  #endif
  for(int i = 0; i < nQueryVertex; i++) 
  {
    delete[] element[i].failingSet;
    #ifdef PRUNING_BY_EQUIVALENCE_SETS
    delete[] element[i].conflictCell;
    delete[] element[i].isPruned;
    delete[] element[i].pruned;
    #ifdef SUBGRAPH_MATCHING
    delete[] element[i].nonConflictCell;
    delete[] element[i].neverVisited;
    delete[] element[i].nFoundMatch;
    #endif
    #endif
    #ifdef LEAF_ADAPTIVE_MATCHING
    delete[] element[i].problemChild;
    #endif
	  #ifdef CANDIDATE_VERTEX_WITH_MAX
    delete[] element[i].isVisited;
	  #endif
    for(int j = 0; j < MAX_QUERY_DEGREE; j++) 
      delete[] iec[i][j];
    //if(useFailingSet)
      delete[] ancestors[i];
    delete[] element[i].visEqSetId;
  }

  // For BuildDAG
  for(int i = 0; i < nQueryVertex; i++)
  {
    #ifdef NEIGHBOR_LABEL_FREQUENCY
    delete[] NLFArray[i];
    #endif
    delete[] DAG_child_query[i];
    delete[] DAG_parent_query[i];
    delete[] DAG_child_query_parent_index[i];
    delete[] DAG_ancestor[i];
    delete[] DAG_ngb_query_ngb_index[i];
    delete[] DAG_ngb_query[i];
  }

  // For ConstructAdjacencyList
  for(int u = 0; u < nQueryVertex; ++u)
  {
    CandidateSpace& currSet = candSpace[u];
    delete[] currSet.candidates;
    delete[] currSet.weight;
    #ifdef N_OPTIMIZATION
    for(int k = 0; k < query.maxDegree; ++k)
    {
      delete[] currSet.nAdjacent[k];
      delete[] currSet.adjacent[k];
      delete[] currSet.capacity[k];
    }
    delete[] currSet.capacity;
    #else
    delete[] currSet.nParentCand;
    #endif
    delete[] currSet.nAdjacent;
    delete[] currSet.adjacent;  	
    #ifdef PRUNING_BY_EQUIVALENCE_SETS
    delete[] currSet.cell;
    delete[] currSet.cellVertex;
    #endif
    
    delete[] currSet.components_cell;
    currSet.intervals_cell.clear();
  }
  #ifdef PRUNING_BY_EQUIVALENCE_SETS
  delete[] posArray;
  delete[] candArray;
  delete[] candOffset;
  #endif
  for(int u = 0; u < nQueryVertex; ++u)
  {
    CandidateSpace& currSet = candSpace[u];
    delete[] currSet.eqSetId;
  }
}

inline void AllocateForQueryGraph(Graph& query)
{
  nQueryVertex = query.nVertex;
  query.NECMapping = new int[nUniqueLabel * nQueryVertex];
  query.NECElement = new Element[nUniqueLabel * nQueryVertex];
  query.NECMap = new int[nQueryVertex];
  #ifdef LEAF_ADAPTIVE_MATCHING
  query.isProblemLeaf = new bool[nQueryVertex];
  #endif
  FAILING_SET_SIZE = (nQueryVertex + 63) >> 6;
  
  leafCandSize = nQueryVertex * maxNumDataVertex;
  leafCand.resize(leafCandSize); // tight upper bound
  leafCandInfo = new pair<int, int>[(nUniqueLabel + 1) * nQueryVertex];
  necMapping = new int[nQueryVertex + 1];
  uCand = new int[nQueryVertex * maxNumCandidate]; //definite upper bound

  pairU = new int[nQueryVertex + 1];
  
  uCandInfo = new pair<int, int>[nQueryVertex * maxNumCandidate]; //definite upper bound

  memset(localFlagQuery, 0, sizeof(char) * nQueryVertex);

  for(int i = 0; i < maxNumDataVertex; i++)
    vCandInfo[i].resize(nQueryVertex);
 
  #ifdef PRUNING_BY_EQUIVALENCE_SETS
  nMaxCell = maxNumCandidate;
  CELL_SIZE = ((maxNumCandidate + 63) >> 6);
  #ifdef SUBGRAPH_MATCHING
  nActive = new int[maxNumDataVertex];
  active = new pair<int, int>*[maxNumDataVertex];
  unordered_map<int, int>::const_iterator lit;
  query.maxNumSameLabelVertex = 0;
  for(lit = query.labelFrequency.begin(); lit != query.labelFrequency.end(); ++lit)
  {
    if(lit->second > query.maxNumSameLabelVertex)
      query.maxNumSameLabelVertex = lit->second;
  }
  for(int i = 0; i < maxNumDataVertex; ++i)
    active[i] = new pair<int, int>[query.maxNumSameLabelVertex];
  #endif
  #endif
  for(int i = 0; i < nQueryVertex; i++) 
  {
    element[i].failingSet = new uint64_t[FAILING_SET_SIZE];
    #ifdef PRUNING_BY_EQUIVALENCE_SETS
    element[i].conflictCell = new uint64_t[CELL_SIZE]();
    element[i].isPruned = new bool[maxNumCandidate](); //TODO. to maxDegree
    element[i].nPruned = 0;
    element[i].pruned = new int[maxNumCandidate];
    #ifdef SUBGRAPH_MATCHING
    element[i].nonConflictCell = new uint64_t[CELL_SIZE]();
    element[i].neverVisited = new bool[maxNumCandidate](); //TODO. to maxDegree
    element[i].nFoundMatch = new double[maxNumCandidate]();
    #endif
    #endif
    #ifdef LEAF_ADAPTIVE_MATCHING
    element[i].problemChild = new int[query.maxDegree];
    element[i].nProblemChild = 0;
    #endif
	  #ifdef CANDIDATE_VERTEX_WITH_MAX
    element[i].isVisited = new bool[maxNumCandidate](); //TODO. to maxDegree
    #endif
	  for(int j = 0; j < MAX_QUERY_DEGREE; j++) 
      iec[i][j] = new int[maxDegree];
    //if(useFailingSet)
      ancestors[i] = new uint64_t[FAILING_SET_SIZE];
    element[i].visEqSetId = new long long[maxNumCandidate];
  }

  // For BuildDAG
  for(int i = 0; i < nQueryVertex; i++)
  {
    #ifdef NEIGHBOR_LABEL_FREQUENCY
    NLFArray[i] = new pair<int, int>[query.degree[i]];
    #endif
    DAG_child_query[i] = new int[query.maxDegree];
    DAG_parent_query[i] = new int[query.maxDegree];
    DAG_child_query_parent_index[i] = new int[query.maxDegree];
    DAG_ancestor[i] = new uint64_t[FAILING_SET_SIZE];
    DAG_ngb_query[i] = new int[query.maxDegree];
    DAG_ngb_query_ngb_index[i] = new int[query.maxDegree];
  }

  // For ConstructAdjacencyList
  for(int u = 0; u < nQueryVertex; ++u)
  {
    CandidateSpace& currSet = candSpace[u];
    currSet.candidates = new int[maxNumCandidate];
    currSet.weight = new long long[maxNumCandidate];
    currSet.nAdjacent = new int*[query.maxDegree];
    currSet.adjacent = new int**[query.maxDegree];  	
    #ifdef N_OPTIMIZATION
    currSet.capacity = new int*[query.maxDegree];
    for(int k = 0; k < query.maxDegree; ++k)
    {
      currSet.nAdjacent[k] = new int[maxNumCandidate];
      currSet.adjacent[k] = new int*[maxNumCandidate];
      currSet.capacity[k] = new int[maxNumCandidate];
      memset(currSet.capacity[k], 0, sizeof(int) * maxNumCandidate);
    }
    #else
    currSet.nParentCand = new int[query.maxDegree];
    memset(currSet.nParentCand, 0, sizeof(int) * query.maxDegree);
    #endif
    #ifdef PRUNING_BY_EQUIVALENCE_SETS
    currSet.cell = new int[maxNumCandidate];
    currSet.cellVertex = new int[maxNumCandidate];
    #endif
    currSet.intervals_cell.resize(maxNumCandidate);
    currSet.components_cell = new int[maxNumCandidate];
  }
  #ifdef PRUNING_BY_EQUIVALENCE_SETS
  posArray = new int[maxNumCandidate];
  candArray = new int[maxNumCandidate];
  candOffset = new int[maxNumCandidate];
  #endif
  for(int u = 0; u < nQueryVertex; ++u)
  {
    CandidateSpace& currSet = candSpace[u];
    currSet.eqSetId = new int[maxNumCandidate];
  }
}

inline void SetQueryGraphResource(Graph& query)
{
  #ifdef PRUNING_BY_EQUIVALENCE_SETS
  globalCellID = 1;
  cellToID.clear();
  #ifdef SUBGRAPH_MATCHING
  nMatch = 0;
  memset(nActive, 0, sizeof(int) * maxNumDataVertex);
  #ifdef MULTIPLE_QUERIES
  lastQuotient = 0;
  cellPos.clear();
  cellList.clear();
  #endif
  cellPos.push_back(0);
  #else
  cellPos[0] = 0;
  #endif
  #endif
  memset(visitedForQuery, 0, sizeof(char) * nQueryVertex);
  memset(cntLabel, 0, sizeof(int) * nUniqueLabel);
  for(int u = 0; u < nQueryVertex; ++u)
  {
    candSpace[u].size = 0;
    #ifdef PRUNING_BY_EQUIVALENCE_SETS
    candSpace[u].nCellVertex = 0;
    memset(candSpace[u].cell, -1, sizeof(int) * maxNumCandidate);
    #endif
  }
}

inline bool FilterByCount(const Graph& query, const Graph& data)
{
  if(query.nVertex > data.nVertex 
    || query.nEdge > data.nEdge 
    || query.maxDegree > data.maxDegree 
    || query.labelFrequency.size() > data.labelFrequency.size())
    return false;
  unordered_map<int, int>::const_iterator qIter, dIter;
  for(qIter = query.labelFrequency.begin(); qIter != query.labelFrequency.end(); ++qIter)
  {
    dIter = data.labelFrequency.find(qIter->first);
    if(dIter == data.labelFrequency.end()) return false;
    else if(qIter->second > dIter->second) return false;
  }
  return true;
}

inline int SelectRootVertex(const Graph& query, const Graph& data)
{
  double leastRank = DBL_MAX, rank;
  int root = -1, label, degree;
  for (int i = 0; i < nQueryVertex; i++)
  {
    if(!isTree && query.core[i] < 2) continue;
    label = query.label[i];
    degree = query.degree[i];
    int start = data.vertexOffset[label];
    int end = data.vertexOffset[label + 1];
    while(start < end && degree > data.degree[data.vertex[start]]) start++;
    rank = (double)(end - start) /(double)degree;
    if(rank < leastRank)
    {
      leastRank = rank;
      root = i;
    }
  }
  return root;
}

inline void FindNEC(Graph& query, int u)
{
  pair<int, int> p;
  unordered_map<pair<int, int>, int, pairHash>::iterator sizeIter;
  #ifdef LEAF_ADAPTIVE_MATCHING
  unordered_map<int, int>::const_iterator freqIter;
  #endif
  //for all of node i's children
  for (int j = query.nbrOffset[u]; j < query.nbrOffset[u+1]; j++)
  {
    int child = query.nbr[j];
    // If child is not in the two-core
    if (query.core[child] < 2) 
    { 
      // If child has degree one (i.e., NEC vertex)
      if (query.degree[child] == 1)
      {
        int label = query.label[child];
        p = make_pair(label, u);
        sizeIter = NECSize.find(p);
        if(sizeIter == NECSize.end()) 
        {
        // child is the representative node
          query.NECElement[query.numNECMapping++] = Element(label, u, child);
          query.NECMap[child] = child;
          NECChild[p] = child;
          NECSize[p] = 1;
          #ifdef LEAF_ADAPTIVE_MATCHING
          freqIter = query.labelFrequency.find(label);
          if(freqIter != query.labelFrequency.end() && freqIter->second > 1)
            query.isProblemLeaf[child] = true;
          #endif
        }
        else
        {
          query.NECMap[child] = NECChild[p];
          sizeIter->second += 1; 
        }
        // the label with parent being i, necCount++
      }
      // If child is in tree
      else
      {
        // extract the query tree for extra region candidate extraction, based on DFS
        // also give a DFS-based query sequence at the same time
        int* stack = arrForQuery;
        int posStack = 0;
        visitedForQuery[u] = 1; // this is the start node's parent node (a marked node)
        visitedForQuery[child] = 1; // this is the start node
        stack[posStack ++] = child;
        while(posStack != 0){
          int curr = stack[posStack - 1];
          posStack--;
          for(int m = query.nbrOffset[curr]; m < query.nbrOffset[curr+1]; m++)
          {
            int child = query.nbr[m];
            if(!visitedForQuery[child]) 
            {
              visitedForQuery[child] = 1;
              //======== special treatment here: if a node is a leaf (degree being 1), then put it into nec node set
              if(query.degree[child] == 1)
              {
  	        int label = query.label[child];
  	        p = make_pair(label, curr);
  	        sizeIter = NECSize.find(p);
  	        if(sizeIter == NECSize.end())
  	        {
  	          // child is the repesentive node
  	          query.NECElement[query.numNECMapping++] = Element(label, curr, child);
  	          query.NECMap[child] = child;
  	          NECChild[p] = child;
  	          NECSize[p] = 1;
              #ifdef LEAF_ADAPTIVE_MATCHING
              freqIter = query.labelFrequency.find(label);
              if(freqIter != query.labelFrequency.end() && freqIter->second > 1)
                query.isProblemLeaf[child] = true;
              #endif
  	        }
  	        else
                {
                  query.NECMap[child] = NECChild[p];
                  sizeIter->second += 1; 
                  // increase the NEC size of vertices with 'label' and parent 'curr'
  	            }
                continue;
              }
              stack[posStack ++] = child;
            }
          }
        }
      }
    }
  }
}

inline bool sortByNECLabel(const Element& p1, const Element& p2) 
{
  return (p1.label < p2.label);
}

inline void ExtractResidualStructures(Graph& query)
{
  NECSize.clear();
  NECChild.clear();
  NECElementOffset.clear();
  query.numNECMapping = 0; 
  memset(query.NECMap, -1, sizeof(int) * nQueryVertex);
  #ifdef LEAF_ADAPTIVE_MATCHING
  memset(query.isProblemLeaf, false, sizeof(bool) * nQueryVertex);
  #endif
  if(isTree)
  {
    FindNEC(query, root);
  }
  else
  {
    for(int i = 0; i < nQueryVertex; i++)
  	{
      if(query.core[i] < 2) continue; //not in two-core
      FindNEC(query, i); 
    }
  }
  
  sort(query.NECElement, query.NECElement + query.numNECMapping, sortByNECLabel);
  // NECElement: a list of (label: lID, parentID: uP, representative: u, sum: |NEC(u)|) 
  // sorted in the ascending order of label ID
  
  int sum;
  int lastLabel;
  int NECSetSize = 0;
  if(query.numNECMapping == 1)
  {
    Element & elem = query.NECElement[0];
    int label = elem.label;
    int parentID = elem.parentID;
    int representChild = elem.representative;
    sum = NECSize[make_pair(label, parentID)];
    NECElementOffset.push_back(make_pair(label, NECSetSize));
    query.NECElement[0].sum = sum;
    NECElementOffset.push_back(make_pair(-1, query.numNECMapping)); // redundant currE to set the end
  } 
  else
  {
    for (int i = 0; i < query.numNECMapping; i++) 
    {
      Element & elem = query.NECElement[i];
      int label = elem.label;
      int parentID = elem.parentID;
      int representChild = elem.representative;
      pair<int, int> p = make_pair(label, parentID);
      sum = NECSize[make_pair(label, parentID)];
      if (i == 0) 
      {
        NECElementOffset.push_back(make_pair(label, NECSetSize));
        query.NECElement[NECSetSize++].sum = sum;
        lastLabel = label;
        continue;
      } 
      else if(i == query.numNECMapping - 1)
      {
        if(label != lastLabel)
          NECElementOffset.push_back(make_pair(label, NECSetSize));
        query.NECElement[NECSetSize ++].sum = sum;
        NECElementOffset.push_back(make_pair(-1, query.numNECMapping)); // redunant currE to set the end
        continue;
      }
      if (label != lastLabel) 
      {
        NECElementOffset.push_back(make_pair(label, NECSetSize));
        lastLabel = label;
      }
      query.NECElement[NECSetSize++].sum = sum;
    }
  }
  // NECElemtntOffset[i] = (l, pos) 
  // l: label ID
  // pos: the starting position of NEC vertices with label l in NECElement if l != -1,
  // the last position of NECElement otherwise
}
#ifdef NEIGHBOR_LABEL_FREQUENCY
inline void GenerateNLF(Graph& query)
{
  for(int ui = 0; ui < nQueryVertex; ++ui)
  {
    for(int j = query.nbrOffset[ui]; j < query.nbrOffset[ui+1]; ++j)
    {
      int un = query.nbr[j];
      int label = query.label[un];
      bool newLabel = true;
      for(int k = 0; k < numNLF[ui]; ++k)
      {
        if(label == NLFArray[ui][k].first)
        {
          newLabel = false;
          NLFArray[ui][k].second += 1;
          break;
        }
      }
      if(newLabel)
      {
        NLFArray[ui][numNLF[ui]++] = make_pair(label, 1);
      }
    }
  }
}
#endif
inline void BuildDAG(Graph& query) 
{
  char * uCurr = visitedForQuery;  		//visitedForQuery[i] = 1 if i-node have uCurr from queue.
  memset(uCurr, 0, sizeof(char) * nQueryVertex);
  int * visited = visitedForQueryInt;  		//visited[i] = level if i-node had pushed into queue, where level is it's BFS-level.
  memset(visited, 0, sizeof(int) * nQueryVertex);
  #ifdef NEIGHBOR_LABEL_FREQUENCY
  memset(numNLF, 0, sizeof(int) * nQueryVertex);
  #endif
  memset(DAG_child_query_size, 0, sizeof(int) * nQueryVertex);
  memset(DAG_parent_query_size, 0, sizeof(int) * nQueryVertex);
  memset(DAG_ngb_query_size, 0, sizeof(int) * nQueryVertex);
  for(int i = 0; i < nQueryVertex; i++)  
    memset(DAG_ancestor[i], 0, sizeof(uint64_t) * FAILING_SET_SIZE);
  arrForQuery[0] = root;
  int startPtr = 0, endPtr = 1, nextEndPtr = 1;
  uSequenceSize = 0;
  visited[root] = 1;
  // Sort labels in the ascending order of label frequency

  while (true) 
  {
    //<sort by label_freq first, and sort by degree for each label group>
    sort(arrForQuery + startPtr, arrForQuery + endPtr, sortQueryByDegree);
    stable_sort(arrForQuery + startPtr, arrForQuery + endPtr, sortQueryByLabelFrequency); 
    while (startPtr != endPtr) 
    { 
      int currNode = arrForQuery[startPtr++];
      uCurr[currNode] = 1;
      for (int i = query.nbrOffset[currNode]; i < query.nbrOffset[currNode + 1]; ++i)
      {
        int childNode = query.nbr[i];
        if (uCurr[childNode] == 0) //childNode is not currNode's parent
        {
          {
            DAG_child_query[currNode][ DAG_child_query_size[currNode] ] = childNode;
            DAG_child_query_parent_index[currNode][ DAG_child_query_size[currNode] ] = DAG_parent_query_size[childNode];
            DAG_parent_query[childNode][ DAG_parent_query_size[childNode] ] = currNode;	
            DAG_ancestor[childNode][currNode >> 6] |= (1ULL << (currNode & 0x3f));
            for(int x = 0; x < FAILING_SET_SIZE; ++x)
              DAG_ancestor[childNode][x] |= DAG_ancestor[currNode][x];
            DAG_child_query_size[currNode]++;
            DAG_parent_query_size[childNode]++;				
          }
          {
            int& currIdx = DAG_ngb_query_size[currNode];
            int& childIdx = DAG_ngb_query_size[childNode];

            DAG_ngb_query[currNode][currIdx] = childNode;
            DAG_ngb_query[childNode][childIdx] = currNode;

            DAG_ngb_query_ngb_index[currNode][currIdx] = childIdx;
            DAG_ngb_query_ngb_index[childNode][childIdx] = currIdx;
            
            ++currIdx;
            ++childIdx;
          }
        }
        if (visited[childNode] == 0) //this child node has not been visited.
        {
          visited[childNode] = 1;
          arrForQuery[nextEndPtr] = childNode;
          nextEndPtr++;
          if(query.core[childNode] < 2)
            continue;
        }
      }
      orderDAG[currNode] = uSequenceSize;
      uSequence[uSequenceSize++] = currNode;
    }
    if(nextEndPtr == endPtr) break;
    startPtr = endPtr;
    endPtr = nextEndPtr;
  }
  query.nNonLeaf = 0;
  for(int i = 0; i < nQueryVertex; ++i)
  {
    if(query.NECMap[i] == -1) ++query.nNonLeaf;
  }
  #ifdef NEIGHBOR_LABEL_FREQUENCY
  GenerateNLF(query);
  #endif
}

enum direction{
  topDown, 
  bottomUp
};

enum ngbType{
  parentNgb,
  childNgb,
  allNgb
};

inline bool filteringWithDAG(
  const Graph& query, 
  const Graph& data,
  direction filterDir,
  ngbType filterNgb,
  bool useNbrSafety
)
{
  unordered_map<pair<int, int>, pair<number_type, number_type>, pairHash>::const_iterator iter;
  
  int orderSt, orderEd, orderAdd;
  
  if(filterDir == topDown)
  {
    orderSt = 0;
    orderEd = uSequenceSize;
    orderAdd = 1;
  } 
  else
  {
    orderSt = uSequenceSize - 1;
    orderEd = -1;
    orderAdd = -1;
  }
  for(
    int i=orderSt;
    i!=orderEd;
    i+=orderAdd
  ){
    ++index_vis_adj;      
    ++index_vis_color;      
    const int query_node_cur = uSequence[i];
    if(query.NECMap[query_node_cur] != -1 
    && query.NECMap[query_node_cur] != query_node_cur)
      continue;
    const int color_cur = query.label[query_node_cur];
    index_set_color = 0;
    int s = query.nbrOffset[query_node_cur];
    int d = query.degree[query_node_cur];

    //1.1 cnt the number of color      
    for(int idx_node_ngb = s; idx_node_ngb < s + d; idx_node_ngb++){
      const int query_node_ngb = query.nbr[idx_node_ngb];
      const int color_ngb = query.label[query_node_ngb];

      if(color[color_ngb][0] == index_vis_color){
        color[color_ngb][1]++;
      }
      else{
        color[color_ngb][0] = index_vis_color;
        color[color_ngb][1] = 1;
        set_color[index_set_color++] = color_ngb;
      }

      if(query.NECMap[query_node_ngb] != -1 && query.NECMap[query_node_ngb] != query_node_ngb)
      {
        continue;
      }

      //1.1.2 see the whole cs of neighbor
      const CandidateSpace &node_unit_ngb = candSpace[query_node_ngb];
      for(int y=0;y<node_unit_ngb.size;y++){
        const int cand_ngb = node_unit_ngb.candidates[y];
        if(vis_adj[cand_ngb] != index_vis_adj){
          vis_adj[cand_ngb] = index_vis_adj;
        }
      }

      //1.1.3 filtering
      
    }
   
    //1.2 calculate parent existence part
    ++index_ngb_existence;
    int nVisitedNeighbor = 0;
    for(int y = 0; y < query.degree[query_node_cur]; y++) // Filtering rule 3.1
    {
      const int query_node_nbr = query.nbr[query.nbrOffset[query_node_cur] + y]; // Filtering rule 3.1
      
      if(
        filterNgb == parentNgb && 
        orderDAG[query_node_cur] < orderDAG[query_node_nbr]   
      )
      {
        continue;
      }
      if(
        filterNgb == childNgb && 
        orderDAG[query_node_cur] > orderDAG[query_node_nbr]   
      )
      {
        continue;
      }
      
      if(query.NECMap[query_node_nbr] != -1 && query.NECMap[query_node_nbr] != query_node_nbr)
        continue;

      CandidateSpace & nbr_unit = candSpace[query_node_nbr];
      for(int z = 0; z < nbr_unit.size; z++){
	    int cand_ngb = nbr_unit.candidates[z];
  #ifdef HUGE_GRAPH
  	  for(long long w = data.nbrOffset[cand_ngb]; w < data.nbrOffset[cand_ngb + 1]; w++)
  #else
        iter = data.labelToNbrOffset.find(make_pair(cand_ngb,color_cur));
        if(iter == data.labelToNbrOffset.end()) continue;
        for(int w=iter->second.first;w<iter->second.first+iter->second.second;w++)
  #endif	
        {
          int cand_around = data.nbr[w];
  #ifdef HUGE_GRAPH
          if(data.label[cand_around] != color_cur) continue;
  #endif
          if(nVisitedNeighbor == 0){
            num_ngb_existence[cand_around][0] = index_ngb_existence;
            num_ngb_existence[cand_around][1] = 1;
          }
          else{
            if(num_ngb_existence[cand_around][0] == index_ngb_existence
            && num_ngb_existence[cand_around][1] == nVisitedNeighbor){
              num_ngb_existence[cand_around][1]++;
            }
          }
        } 
      }
      ++nVisitedNeighbor;
    }

    //1.3 filtering      
    CandidateSpace &node_unit_cur = candSpace[query_node_cur];
    int& new_size = node_unit_cur.size;
    int x = 0;
    while(x < new_size){
      const int cand_cur = node_unit_cur.candidates[x];
      bool iSVaildCand = true;
      //1.3.1 : filtering by parent

	    // Filtering rule 3.1 (in the codition below after ||)
      if(
        (
          num_ngb_existence[cand_cur][0] != index_ngb_existence ||
          num_ngb_existence[cand_cur][1] != nVisitedNeighbor
        ) &&
        (
          nVisitedNeighbor != 0
        )
      )
      {
        iSVaildCand = false;
      }
    
      //1.3.2 : filtering by cardinality 
      if(
        iSVaildCand &&
        useNbrSafety
      ){
  #ifdef HUGE_GRAPH
        int t = 0;
        for(long long z = data.nbrOffset[cand_cur]; z < data.nbrOffset[cand_cur+1]; ++z)
        {
          int ngb = data.nbr[z];
          if(vis_adj[ngb] == index_vis_adj){
              ++t;
          }
        }
	      if(t < d){
          iSVaildCand = false;
        }
  #else
        for(int y = 0;y < index_set_color;y++){
          int t = 0;
          iter = data.labelToNbrOffset.find(make_pair(cand_cur, set_color[y]));
	        if(iter == data.labelToNbrOffset.end()) {
            iSVaildCand = false;
            break;
          }
          for(int z = iter->second.first; z < iter->second.first + iter->second.second; ++z)
          {
            int ngb = data.nbr[z];
            
            if(vis_adj[ngb] == index_vis_adj){
              ++t;
            }
          }
	        if(t < color[set_color[y]][1]){
            iSVaildCand = false;
            break;
          } 
        }
  #endif
      }
      //1.3.3 : remain or remove
      if(!iSVaildCand){
        int y = node_unit_cur.candidates[new_size - 1];
        node_unit_cur.candidates[x] = y;
        --new_size;
      }
      else{
        ++x;
      }
    }
    if(new_size == 0)return false;
    node_unit_cur.size = new_size;
  }
  return true;
}

inline bool TopDownInitial(const Graph& query, const Graph& data) 
{
  // Initialize candidates of root vertex
  CandidateSpace & rootUnit = candSpace[root];
  int label = query.label[root];
  int start = data.vertexOffset[label];
  int end = data.vertexOffset[label + 1];
  unordered_map<pair<int, int>, pair<number_type, number_type>, pairHash>::const_iterator iter;
  while(start < end && query.degree[root] > data.degree[data.vertex[start]]) start++;
  // For each v in G with label
  for (int j = start; j < end; j++) 
  {
    int v = data.vertex[j]; //v
    if(data.core[v] < query.core[root] 
    || data.maxNbrDegree[v] < query.maxNbrDegree[root])
      continue;
    char isAdded = 1;
    #ifdef NEIGHBOR_LABEL_FREQUENCY
    const unordered_map<int, int>& dataNLF = data.NLF[v];
    if(dataNLF.size() < numNLF[root]) continue;
    pair<int, int> currNLF;
    unordered_map<int, int>::const_iterator iter;
    for(int k = 0; k < numNLF[root]; ++k)
    {
      currNLF = NLFArray[root][k];
      iter = dataNLF.find(currNLF.first);
      if(iter == dataNLF.end() || iter->second < currNLF.second)
      {
        isAdded = 0; break;
      }
    }  
    #else
    for (int pos = NLFSize - 1; pos >= 0; pos--)
    {
      #ifdef HUGE_GRAPH
      if(data.NLF[(long long)v * (long long)NLFSize + (long long)pos] != (query.NLF[root * NLFSize + pos] | data.NLF[(long long)v * (long long)NLFSize + (long long)pos]))
      #else
      if(data.NLF[v * NLFSize + pos] != (query.NLF[root * NLFSize + pos] | data.NLF[v * NLFSize + pos]))
      #endif
      {
        isAdded = 0; break;
      }
    }
    #endif
    if (isAdded)
      rootUnit.candidates[rootUnit.size++] = v;
  }
  if(rootUnit.size == 0) return false;
  // For each query vertex u of q in a top-down fashion
  for (int i = 1; i < uSequenceSize; i++)
  {
    // i = 0 is root, which is already processed above.
    int currVertex = uSequence[i]; //u
    if(query.NECMap[currVertex] != -1 
    && query.NECMap[currVertex] != currVertex)
      continue;
    int currLabel = query.label[currVertex];
    char checkVal = 0;  //Cnt = 0

    // For each parent p of u in dag
    for(int parentPos = 0; parentPos < DAG_parent_query_size[currVertex]; parentPos++)
    { 
      int parent = DAG_parent_query[currVertex][parentPos]; //p
      CandidateSpace & parentUnit = candSpace[parent];
      // For each candidate vertex v' in p.C
      for (int y = 0; y < parentUnit.size; y++) 
      {
        int parentCand = parentUnit.candidates[y]; //v'
        // For each vertex v adjacent to v' with label
	      #ifdef HUGE_GRAPH
        for(long long z = data.nbrOffset[parentCand]; z < data.nbrOffset[parentCand + 1]; ++z)
        {
          int canID = data.nbr[z];
	        if(data.label[canID] != currLabel) continue;
          if(isVisited[canID] == checkVal)
          {
            isVisited[canID] ++; //v.cnt++;
            if(checkVal == 0)
              arrayToClean[toCleanIndex++] = canID;
          }
        }
	      #else
        iter = data.labelToNbrOffset.find(make_pair(parentCand, currLabel));
        if(iter == data.labelToNbrOffset.end()) continue;
        for(int z = iter->second.first; z < iter->second.first + iter->second.second; ++z)
        {
          int canID = data.nbr[z];
          if(isVisited[canID] == checkVal)
          {
            isVisited[canID] ++; //v.cnt++;
            if (checkVal == 0)
              arrayToClean[toCleanIndex++] = canID;
          }
        }
	      #endif
      }
      checkVal++; //update the check value by one
    }
    // For each vertex v with v.cnt = Cnt, if v passes filters, add candidates
    CandidateSpace & currSet = candSpace[currVertex];
    int candIndex = 0;
    for(int checkIndex = 0; checkIndex < toCleanIndex; checkIndex++)
    { 
      int canID = arrayToClean[checkIndex]; //v
      if(isVisited[canID] == checkVal)
      {
        //check degree, core, and max_neighbor degree together
        if(data.degree[canID] < query.degree[currVertex] 
        || data.core[canID] < query.core[currVertex] 
        || data.maxNbrDegree[canID] < query.maxNbrDegree[currVertex])
          continue;
        //check lightweight NLF
        char isAdded = 1;
        #ifdef NEIGHBOR_LABEL_FREQUENCY
        const unordered_map<int, int>& dataNLF = data.NLF[canID];
        if(dataNLF.size() < numNLF[currVertex]) continue;
        pair<int, int> currNLF;
        unordered_map<int, int>::const_iterator iter;
        for(int k = 0; k < numNLF[currVertex]; ++k)
        {
          currNLF = NLFArray[currVertex][k];
          iter = dataNLF.find(currNLF.first);
          if(iter == dataNLF.end() || iter->second < currNLF.second)
          {
            isAdded = 0; break;
          }
        }	
        #else
        for (int pos = NLFSize - 1; pos >= 0; pos--)
        {
	        #ifdef HUGE_GRAPH
          if(data.NLF[(long long)canID * (long long)NLFSize + (long long)pos] != (query.NLF[currVertex * NLFSize + pos] | data.NLF[(long long)canID * (long long)NLFSize + (long long)pos] ))
	        #else
          if(data.NLF[canID * NLFSize + pos] != (query.NLF[currVertex * NLFSize + pos] | data.NLF[canID * NLFSize + pos] ))
	        #endif
          {
            isAdded = 0; break;
          }
        }
        #endif
        // lightweight NLF OK
        if (isAdded){
          currSet.candidates[candIndex] = canID;
          candIndex++;
        }
      }
    }
    currSet.size = candIndex;
    // Reset v.cnt for all v in G such that v.cnt > 0  
    while(toCleanIndex != 0)
      isVisited[arrayToClean[--toCleanIndex]] = 0; 
    if(currSet.size == 0) return false;
  } // In a top-down manner
  return true;
}

inline bool BottomUp(const Graph& query, const Graph& data)
{
  unordered_map<pair<int, int>, pair<number_type, number_type>, pairHash>::const_iterator iter;
  // For each query vertex u of dag(q) in a bottom-up fashion
  for (int i = uSequenceSize-1; i >= 0; --i)
  { 
    int currVertex = uSequence[i]; //u
    //No child, No prunning.
    if(DAG_child_query_size[currVertex] == 0) 
      continue;
    int currLabel = query.label[currVertex];
    char checkVal = 0;  //Cnt = 0
    // For each child uc of u in dag(q)
    for(int childPos = 0; childPos < DAG_child_query_size[currVertex]; childPos++)
    { 
      int child = DAG_child_query[currVertex][childPos]; //uc
      CandidateSpace& childSet = candSpace[child];
      //process only represent node when the child is in NEC
      if(query.NECMap[child] != -1 && query.NECMap[child] != child) 
	      continue;
      // For each candidate vc in uc.C
      for (int y = 0; y < childSet.size; y++) 
      {
        int childCand = childSet.candidates[y]; //vc
        // For each vertex v adjacent to vc
	      #ifdef HUGE_GRAPH
        for(long long z = data.nbrOffset[childCand]; z < data.nbrOffset[childCand + 1]; ++z)
	      #else
        iter = data.labelToNbrOffset.find(make_pair(childCand, currLabel));
        if(iter == data.labelToNbrOffset.end()) continue;
        for(int z = iter->second.first; z < iter->second.first + iter->second.second; ++z)
	      #endif
        {
          int canID = data.nbr[z]; //v
	        #ifdef HUGE_GRAPH
	        if(data.label[canID] != currLabel) continue;
	        #endif
          if(isVisited[canID] == checkVal)
          {
            isVisited[canID]++;
            if(checkVal == 0)
	            arrayToClean[toCleanIndex++] = canID;
          }
        }
      }
      checkVal++; //update the check value by one
    }
    // For each vertex v in u.C  
    CandidateSpace & currSet = candSpace[currVertex];
    for(int candIndex = 0; candIndex < currSet.size; candIndex++)
    { 
      int canID = currSet.candidates[candIndex]; //v
      //if v.cnt != Cnt then Remove v from u.C
      if(isVisited[canID] != checkVal)
      { 
        currSet.candidates[candIndex] = currSet.candidates[--currSet.size];
        --candIndex;
      }
    }
    // Reset v.cnt = 0 for all v in G such that v.cnt > 0
    while(toCleanIndex != 0)
      isVisited [arrayToClean[--toCleanIndex]] = 0;
    if(currSet.size == 0) return false;
  } // In a bottom-up manner
  return true;
}

inline bool TopDown(const Graph& query, const Graph& data)
{
  unordered_map<pair<int, int>, pair<number_type, number_type>, pairHash>::const_iterator iter;
  // For each query vertex u of dag(q) in a top-down fashion
  for (int i = 1; i < uSequenceSize; i++)
  { 
    // root(0) has no parent
    int currVertex = uSequence[i]; //u
    if(query.NECMap[currVertex] != -1 && query.NECMap[currVertex] != currVertex)
      continue;
    int currLabel = query.label[currVertex];
    char checkVal = 0;  //Cnt = 0
    // For each parent up of u in dag(q)
    for(int parentPos = 0; parentPos < DAG_parent_query_size[currVertex]; parentPos++)
    { 
      int parent = DAG_parent_query[currVertex][parentPos]; //up
      CandidateSpace & parentUnit = candSpace[parent];
      // For each candidate vp in up.C
      for (int y = 0; y < parentUnit.size; y++) 
      {
        int parentCand = parentUnit.candidates[y]; //vp
        // For each vertex v adjacent to vp
        #ifdef HUGE_GRAPH
        for(long long z = data.nbrOffset[parentCand]; z < data.nbrOffset[parentCand + 1]; ++z)
        {
          int canID = data.nbr[z]; //v
	        if(data.label[canID] != currLabel) continue;
          if (isVisited[canID] == checkVal)
          {
            isVisited[canID]++;
            if(checkVal == 0)
              arrayToClean[toCleanIndex++] = canID;
          }
        }
        #else
        iter = data.labelToNbrOffset.find(make_pair(parentCand, currLabel));
        if(iter == data.labelToNbrOffset.end()) continue;
        for(int z = iter->second.first; z < iter->second.first + iter->second.second; ++z)
        {
          int canID = data.nbr[z]; //v
          if (isVisited[canID] == checkVal)
          {
            isVisited[canID]++;
            if(checkVal == 0)
              arrayToClean[toCleanIndex++] = canID;
          }
        }
        #endif
      }
      checkVal++; //update the check value by one
    }
    // For each vertex v in u.C  
    CandidateSpace & currSet = candSpace[currVertex];
    for(int candIndex = 0; candIndex < currSet.size; candIndex++)
    { 
      int canID = currSet.candidates[candIndex]; //v
      //if v.cnt != Cnt then Remove v from u.C
      if(isVisited[canID] != checkVal) 
      {
        currSet.candidates[candIndex] = currSet.candidates[ --currSet.size];
        --candIndex;
      }
    }
    // Reset v.cnt = 0 for all v in G such that v.cnt > 0
    while(toCleanIndex != 0) 
      isVisited[arrayToClean[--toCleanIndex]] = 0;
    if(currSet.size == 0) return false;
  } // In a top-down manner
  return true;
}

inline bool ConstructAdjacencyList(const Graph& query, const Graph& data, bool filter)
{
  unordered_map<pair<int, int>, pair<number_type, number_type>, pairHash>::const_iterator iter;
  // For each query vertex u of dag(q) in a bottom-up fashion
  for(int i = uSequenceSize - 1; i >= 0; --i)
  { 
    int u = uSequence[i]; //u
    int currLabel = query.label[u];
    if(query.NECMap[u] != -1 && query.NECMap[u] != u) continue;
    CandidateSpace& currSet = candSpace[u];
    //No child, No adjacency list (adjacency list is stored in the child)
    if(filter)
    {
      // <GlobalLabelFrequency>
      for(int vPos = 0; vPos < currSet.size; ++vPos)
      {    
        int currCand = currSet.candidates[vPos];
        if(isVisited[currCand] == 0)
        {
          isVisited[currCand] = 1;
          arrayToClean[toCleanIndex++] = currCand;
          ++cntLabel[currLabel];
        }
      }
      // </GlobalLabelFrequency>
    }
    //if (DAG_child_query_size[u] == 0) continue;
    
    // Preprocess u.C 1) to check (if v in u.C) in constant time 2) to know v's posCand
    for(int vPos = 0; vPos < currSet.size;++vPos)
    { 
      candToPos[currSet.candidates[vPos]] = vPos;
    }
    //1.1. for each neighbor un of u in q
    for(int unPos = 0; unPos < DAG_ngb_query_size[u]; ++unPos)
    {
      int un = DAG_ngb_query[u][unPos]; //uc
      CandidateSpace& nbrSet = candSpace[un];
      int uPos = DAG_ngb_query_ngb_index[u][unPos];
      
      // Process only represent node when the child is in NEC
      if (query.NECMap[un] != -1 && query.NECMap[un] != un) continue;
      // Initialize candidate space
      // Make sure it won't allocate new array every time
      #ifdef N_OPTIMIZATION
      #else
      if (nbrSet.nParentCand[uPos] < currSet.size)
      {
        nbrSet.nAdjacent[uPos] = new int[currSet.size];
        nbrSet.adjacent[uPos] = new int*[currSet.size];
        nbrSet.nParentCand[uPos] = currSet.size;
      }
      #endif
      memset(nbrSet.nAdjacent[uPos], 0 , sizeof(int) * currSet.size);
      #ifdef N_OPTIMIZATION
      for (int vPos = 0; vPos < currSet.size; ++vPos)
      {
        int currSize = data.degree[currSet.candidates[vPos]];
        if(nbrSet.capacity[uPos][vPos] < currSize)
        {
          nbrSet.capacity[uPos][vPos] = currSize;
          nbrSet.adjacent[uPos][vPos] = new int[currSize];
        }
      }
      #else
      for (int vPos = 0; vPos < currSet.size; ++vPos)
      {
        nbrSet.adjacent[uPos][vPos] = new int[data.degree[currSet.candidates[vPos]]];
      }
      #endif
      //1.1.1. for each vertex vn in C(un) Most time-consuming in ConstructAdjacencyList
      for (int vnPos = 0; vnPos < nbrSet.size; ++vnPos)
      { 
        int vn = nbrSet.candidates[vnPos];
        // For each v adjacent to vc in G with label
	      #ifdef HUGE_GRAPH
        for(long long z = data.nbrOffset[vn]; z < data.nbrOffset[vn + 1]; ++z)
        {
          // For each vertex v in N_G(v_p) with label l_q(u)
          int v = data.nbr[z]; //v
	        if(data.label[v] != currLabel) continue;
          int vPos = candToPos[v];
	        if(vPos < 0) continue;
          nbrSet.adjacent[uPos][vPos][nbrSet.nAdjacent[uPos][vPos]++] = vnPos;
        }
	      #else
        iter = data.labelToNbrOffset.find(make_pair(vn, currLabel));
        if(iter == data.labelToNbrOffset.end()) continue;
        for(int z = iter->second.first; z < iter->second.first + iter->second.second; ++z)
        {
          // For each vertex v in N_G(v_p) with label l_q(u)
          int v = data.nbr[z]; //v
          int vPos = candToPos[v];
	        if(vPos < 0) continue;
          nbrSet.adjacent[uPos][vPos][nbrSet.nAdjacent[uPos][vPos]++] = vnPos;
        }
        #endif
      }
    }
    // Reset v in G such that v in u.C and v.index = v's index in u.C
    for(int vPos = 0; vPos < currSet.size; ++vPos)
    {
      candToPos[currSet.candidates[vPos]] = -1;
    }
  }
  if(filter)
  {
    //<GlobalLabelFrequency>
    while(toCleanIndex != 0)
      isVisited[arrayToClean[--toCleanIndex]] = 0;
    for(unordered_map<int, int>::const_iterator it = query.labelFrequency.begin(); it != query.labelFrequency.end(); ++it)
    {
      if(cntLabel[it->first] < it->second)
        return false;
    }
    //</GlobalLabelFrequency>
  }
  return true;
}

#ifdef SUBGRAPH_MATCHING
inline void ComputeAuxiliaryDataStructureSize()
{
  auxDataStructureSize = 0;
  for(int i = 0; i < uSequenceSize; ++i)
  {
    int u = uSequence[i];
    if(currQ->NECMap[u] != -1 && currQ->NECMap[u] != u) 
      continue;
    auxDataStructureSize += candSpace[u].size;
  }
}
#endif

inline int partitionByParent(int arr[], int index, int nbr, int nbrPos, int low, int high)
{
  int pos, i = low - 1, nVisited = 0;
  for(int j = 0; j < candSpace[uI].nAdjacent[index][nbrPos]; ++j)
  {
    pos = candSpace[uI].adjacent[index][nbrPos][j];
    arrayToClean[nVisited++] = pos; 
    candToPos[pos] = j;
  }
  for(int j = low; j < high; ++j)
  {
    if(candToPos[arr[j]] < 0) continue;
    ++i;
    swap(arr[i], arr[j]);
  }
  while(nVisited != 0)
    candToPos[arrayToClean[--nVisited]] = -1;
  return (i + 1);
}

inline int partitionByChild(int arr[], int uPos, int nbr, int nbrPos, int low, int high)
{
  int i = low - 1;
  for(int j = low; j < high; ++j)
  {
    bool exist = false;
    while(posArray[arr[j]] < candSpace[nbr].nAdjacent[uPos][arr[j]] &&
    candSpace[nbr].adjacent[uPos][arr[j]][posArray[arr[j]]] == nbrPos)
    {
      exist = true;
      ++posArray[arr[j]];
    }
    if(exist) swap(arr[++i], arr[j]);
  }
  return (i + 1);
}

inline void sortByParent(int arr[], int index, int nbr, int nbrPos, int low, int high)
{
  if(nbrPos >= candSpace[nbr].size)
  {
    nbrPos = 0; ++index; 
    if(index >= DAG_parent_query_size[uI]) return;
    nbr = DAG_parent_query[uI][index];
  }
  int pivot = partitionByParent(arr, index, nbr, nbrPos, low, high);
  if(pivot > low && pivot < high) candOffset[nCandOffset++] = pivot;
  //cout<<"Partition by v"<<candSpace[nbr].candidates[nbrPos]<<"("<<nbrPos<<") in C(u"<<nbr<<"). low: "<<low<<". pivot: "<<pivot<<". high: "<<high<<endl;  
  if(pivot - low > 1)
    sortByParent(arr, index, nbr, nbrPos + 1, low, pivot);
  if(high - pivot > 1)
    sortByParent(arr, index, nbr, nbrPos + 1, pivot, high);
}

inline void sortByChild(int arr[], int index, int uPos, int nbr, int nbrPos, int low, int high)
{
  if(nbrPos >= candSpace[nbr].size)
  {
    nbrPos = 0; ++index; 
    if(index >= DAG_child_query_size[uI]) return;
    nbr = DAG_child_query[uI][index];
    while(currQ->NECMap[nbr] != -1 && currQ->NECMap[nbr] != nbr)
    {
      ++index;
      if(index >= DAG_child_query_size[uI]) return;
      nbr = DAG_child_query[uI][index];
    }
    uPos = DAG_child_query_parent_index[uI][index];
    for(int x = low; x < high; ++x) posArray[arr[x]] = 0;
  }
  int pivot = partitionByChild(arr, uPos, nbr, nbrPos, low, high);
  if(pivot > low && pivot < high) candOffset[nCandOffset++] = pivot;
  //cout<<"!Partition by v"<<candSpace[nbr].candidates[nbrPos]<<"("<<nbrPos<<") in C(u"<<nbr<<"). low: "<<low<<". pivot: "<<pivot<<". high: "<<high<<endl;  
  if(pivot - low > 1)
    sortByChild(arr, index, uPos, nbr, nbrPos + 1, low, pivot);
  if(high - pivot > 1)
    sortByChild(arr, index, uPos, nbr, nbrPos + 1, pivot, high);
}

inline void GetCell(const Graph& query, int u, int d, int* extCand, int nExtCand)
{
  int uPos, nbr, nCand = 0;
  CandidateSpace& currSet = candSpace[u];
  if(u == root)
  {
    for(int j = 0; j < currSet.size; ++j) candArray[j] = j;
    nCand = currSet.size;
  }
  else
  {
    int curr = -1;
    for(int j = 0; j < nExtCand; ++j)
    {
      if(currSet.cell[extCand[j]] < 0 && curr != extCand[j]) 
	    candArray[nCand++] = extCand[j];
      curr = extCand[j];
    }
  }
  if(nCand == 0) return;
  //cout<<"Bundling up vertices in C(u"<<u<<")"<<endl; 
  uI = u;
  nCandOffset = 0;
  candOffset[nCandOffset++] = 0;
  candOffset[nCandOffset++] = nCand;
  int nParent = DAG_parent_query_size[u];
  int nChild = DAG_child_query_size[u];
  int index = 0;
  #ifdef CELL_OPTIMIZATION
  #else
  if(nParent > 0)
  {
    nbr = DAG_parent_query[u][index];  
    // Compute cells in C(u)
    sortByParent(candArray, index, nbr, 0, 0, nCand);
    // Sort the offset of cells in candArray in the ascending order
    sort(candOffset, candOffset + nCandOffset);
  }
  #endif
  if(nChild > 0)
  {
    nbr = DAG_child_query[u][index];
    while(query.NECMap[nbr] != -1 && query.NECMap[nbr] != nbr)
    {
      ++index;
      if(index >= nChild) break;
      nbr = DAG_child_query[u][index];
    }
    if(index < nChild)
    {
      uPos = DAG_child_query_parent_index[u][index];
      for(int x = 0; x < nCand; ++x) posArray[candArray[x]] = 0;
      int nPrevCandOffset = nCandOffset;
      for(int j = 0; j < nPrevCandOffset - 1; ++j)
      {
        if(candOffset[j+1] - candOffset[j] == 1)
        {
          currSet.cell[candArray[candOffset[j]]] = 0;
          currSet.cellVertex[currSet.nCellVertex++] = candArray[candOffset[j]];
          continue;
        }
        // Compute cells in C(u)
        sortByChild(candArray, index, uPos, nbr, 0, candOffset[j], candOffset[j + 1]);
      }
      // Sort the offset of cells in candArray in the ascending order
      sort(candOffset, candOffset + nCandOffset);
    }
  }
  // Save cell info in candidate space 
  for(int j = 0 ; j < nCandOffset - 1; ++j)
  {
    if(candOffset[j+1] - candOffset[j] == 1)
    {
      currSet.cell[candArray[candOffset[j]]] = 0;
      currSet.cellVertex[currSet.nCellVertex++] = candArray[candOffset[j]];
      continue;
    }
    Array* a = new Array(currSet.candidates, candArray, candOffset[j], candOffset[j+1]);
    uIter = cellToID.find(*a);
    int currID;
    if(uIter == cellToID.end())
    {
      currID = cellToID[*a] = globalCellID++;
      if(globalCellID >= nMaxCell)
      {
        int prevSize = CELL_SIZE;
        nMaxCell *= 2;
        CELL_SIZE = ((nMaxCell + 63) >> 6);
        for(int i = 0; i < nQueryVertex; ++i)
        {
          element[i].buffer = new uint64_t[CELL_SIZE]();
          #ifdef SUBGRAPH_MATCHING
          element[i].auxBuffer = new uint64_t[CELL_SIZE](); 
          #endif
        }
        for(int i = 0; i < d; ++i)
        {
          for(int j = 0; j < prevSize; ++j)
          {
            element[i].buffer[j] = element[i].conflictCell[j];
            #ifdef SUBGRAPH_MATCHING
            element[i].auxBuffer[j] = element[i].nonConflictCell[j];
            #endif
          }
        }
        for(int i = 0; i < nQueryVertex; ++i) 
        {
          element[i].conflictCell = element[i].buffer;
          #ifdef SUBGRAPH_MATCHING
          element[i].nonConflictCell[j] = element[i].auxBuffer[j];
          #endif
        }
      }
      #ifdef SUBGRAPH_MATCHING
      cellPos.push_back(cellPos[currID - 1]);
      #else
      cellPos[currID] = cellPos[currID - 1];
      #endif
      #ifdef MAPPING_FUNCTION_LOG
      cout<<"Cell "<<currID<<":";
      #endif
      for(int k = candOffset[j]; k < candOffset[j+1]; ++k)
      {
        #ifdef MAPPING_FUNCTION_LOG
        cout<<" v"<<currSet.candidates[candArray[k]]<<"("<<cellPos[currID]<<")"; 
        #endif
        #ifdef SUBGRAPH_MATCHING
        cellList.push_back(currSet.candidates[candArray[k]]);
        ++cellPos[currID];
        #else
        cellList[cellPos[currID]++] = currSet.candidates[candArray[k]];
        #endif
      }
      #ifdef MAPPING_FUNCTION_LOG
      cout<<endl;
      #endif
    }
    else
    {
      currID = uIter->second;
    }
    for(int k = candOffset[j]; k < candOffset[j+1]; ++k)
    {
      currSet.cell[candArray[k]] = currID;
      currSet.cellVertex[currSet.nCellVertex++] = candArray[k];
    }	
  }
}

#ifdef SUBGRAPH_MATCHING
inline void Deactivate(Stack& retE, int depth)
{
  bool isRoot = (retE.vertex == root);
  CandidateSpace& retSet = candSpace[retE.vertex];
  for(int aPos = 0; aPos < retE.addressSize; ++aPos)
  {
    int idx = isRoot?aPos:retE.address[aPos];
    int v = retSet.candidates[idx];
    if(nActive[v] > 0 && active[v][nActive[v] - 1].first == depth) 
      --nActive[v];
  }
}
#endif

inline void GetEquivalence(Stack& retE)
{
  bool isRoot = (retE.vertex == root);
  CandidateSpace& retSet = candSpace[retE.vertex];
  int initPos = retE.addressPos;
  int idx = isRoot?initPos:retE.address[initPos];
  int cID = retSet.cell[idx];
  if((retE.conflictCell[0] & 1ULL) != 0 || cID == 0) return;
  int nCell = 1, cPos = 0, aPos = cellPos[cID - 1];
  int v = retSet.candidates[idx];
  #ifdef MAPPING_FUNCTION_LOG
  cout<<"GetEquivalence. ConflictCell(u"<<retE.vertex<<"): "<<bitset<64>(retE.conflictCell[0])<<endl;
  #ifdef SUBGRAPH_MATCHING
  cout<<"GetEquivalence. NonconflictCell(u"<<retE.vertex<<"): "<<bitset<64>(retE.nonConflictCell[0])<<endl;
  #endif
  #endif
  for(int y = cellPos[cID - 1]; y < cellPos[cID]; ++y)
  {
    arrayToClean[toCleanIndex++] = cellList[y];
    ++isVisited[cellList[y]];
  }
  for(int z = 0; z < CELL_SIZE; ++z)
  {
    arr = retE.conflictCell[z]; 
    nCell += __builtin_popcountll(arr);
    while(arr)
    {
      int id = (z << 6) + __builtin_ctzl(arr);
      for(int y = cellPos[id - 1]; y < cellPos[id]; ++y)
      {
        if(isVisited[cellList[y]] == 1)
          ++isVisited[cellList[y]];
      }
      arr ^= (arr & -arr);
    }
  }
  #ifdef SUBGRAPH_MATCHING
  for(int z = 0; z < CELL_SIZE; ++z)
  {
    arr = retE.nonConflictCell[z] & (~retE.conflictCell[z]); 
    while(arr)
    {
      int id = (z << 6) + __builtin_ctzl(arr);
      for(int y = cellPos[id - 1]; y < cellPos[id]; ++y)
	      isVisited[cellList[y]] = 0;
      arr ^= (arr & -arr);
    }
  }
  #endif
  for(aPos = initPos; aPos < retE.addressSize; ++aPos)
  {
    if(retE.isPruned[aPos]) continue;
    idx = isRoot?aPos:retE.address[aPos];
    v = retSet.candidates[idx];
    #ifdef SUBGRAPH_MATCHING
    if((retE.nFoundMatch[initPos] > 0 && isVisited[v] == nCell && retE.neverVisited[aPos])
    || (retE.nFoundMatch[initPos] == 0 && isVisited[v] == nCell))
    #else
    if(isVisited[v] == nCell)
    #endif
    {
      #ifdef SUBGRAPH_MATCHING
      retE.neverVisited[aPos] = false;
      retE.nFoundMatch[aPos] = retE.nFoundMatch[initPos];
      #ifdef MAPPING_FUNCTION_LOG
      cout<<"Check v"<<v<<"("<<aPos<<") \\in C(u"<<retE.vertex<<"). ";
      cout<<"E[u"<<retE.vertex<<"].nFoundMatch["<<aPos<<"]: "<<"E[u"<<retE.vertex<<"].nFoundMatch["<<initPos<<"] ("<<retE.nFoundMatch[initPos]<<")"<<endl;
      #endif
      #endif
      retE.isPruned[aPos] = true;
      retE.pruned[retE.nPruned++] = aPos;
    }
  }
  while(toCleanIndex != 0)
    isVisited[arrayToClean[--toCleanIndex]] = 0;
  #ifdef MAPPING_FUNCTION_LOG
  cout<<"Pruned in C(u"<<retE.vertex<<"):";
  for(int i = 0; i < retE.nPruned; ++i)
  {
    idx = isRoot?retE.pruned[i]:retE.address[retE.pruned[i]];
    cout<<" v"<<retSet.candidates[idx];
  }
  cout<<endl;
  #endif
}

inline void InitializeWeight(const Graph& query)
{
  int v, vcIndex, u, uc;
  long long min, total;
  bool flag = false;
  for(int k = uSequenceSize - 1; k >= 0; --k)
  {
    u = uSequence[k];
    if(query.NECMap[u] != -1) continue;
    CandidateSpace& currSet = candSpace[u];
    for(int i = 0; i < currSet.size; ++i)
    { // for each vertex v in u.C
      v = currSet.candidates[i];     // v = u.C
      if(v < 0) continue;
      min = LLONG_MAX; flag = false;
      for(int j = 0; j < DAG_child_query_size[u]; ++j)
      { // for each child uc of u in dag(q)
        uc = DAG_child_query[u][j];
        if(query.NECMap[uc] != -1) continue;
        if(DAG_parent_query_size[uc] > 1) continue;
        flag = true; total = 0;
        CandidateSpace& childSet = candSpace[uc];
        for(int k = 0; k < childSet.nAdjacent[0][i]; ++k)
        { // for each vc in N^{u}_{uc}(v)
          vcIndex = childSet.adjacent[0][i][k];
          if(vcIndex < 0) continue;
          if(childSet.weight[vcIndex] > LLONG_MAX-1 - total)
          {
            total = LLONG_MAX-1; break;
          }
          total += childSet.weight[vcIndex];
        }
        if(total < min) min = total;
      }
      if(!flag) 
        currSet.weight[i] = 1;
      else
        currSet.weight[i] = min;
    }
  }  
}

inline bool BFS_MM()
{ 
  //BFS function for maximum matching
  /*
   * Note that in this function, we are NOT manipulating the actual query nodes and data nodes.
   * Instead, we are dealing with the index in LLC and the data node index in the sum NEC candidate set.
   */
  int* queue = arrayToClean;
  int queueStart = 0;//use for popping up
  int queueEnd = 0;//use for pushing back
  for(int u = 1; u < nodeIndex; u++)
  {
    if(pairU[u] == 0) //NULL
    {
      dist[u] = 0;
      queue[queueEnd ++] = u;
    }
    else
      dist[u] = INT_MAX;
  }
  dist[0] = INT_MAX; //NULL
  while(queueStart != queueEnd)
  {
    int u = queue[queueStart];
    queueStart++;
    if(dist[u] < dist[0]) //NULL
    {
      int nec = necMapping[u]; //get the actual index in the LLC
      for(int j = 0; j < NECRegionSize[nec]; j++)
      {
        int v = uCand[nec * maxNumCandidate + j];
        if(dist[pairV[v]] == INT_MAX){
          dist[pairV[v]] = dist[u] + 1; 
          //here use "u" instead of "nec", because "nec" is only used to get the candidate set.
          queue[queueEnd++] = pairV[v];
      	}
      }
    }
  }
  return (dist[0] != INT_MAX); //NULL
}

inline bool DFS_MM(int u)
{
  if (u != 0) //NULL
    {
    int nec = necMapping[u]; //get the actual index in the LLC
    for(int j = 0; j < NECRegionSize[nec]; j++)
    {
      int v =  uCand[nec * maxNumCandidate + j];
      //cout<<"uCand["<<nec<<"*"<<maxNumCandidate<<"+"<<j<<"<"<<nQueryVertex<<"*"<<maxNumCandidate<<"]: "<<v<<" add. "<<&uCand[nec*maxNumCandidate+j]<<endl;
      if(dist[pairV[v]] == dist[u] + 1)
        if(DFS_MM(pairV[v]))
        {
          pairV[v] = u;
          pairU[u] = v;
          return true;
        }
    }
    dist[u] = INT_MAX;
    return false;
  }
  return true;
}

inline void combine(double & cntEmbedding, double resultSoFar)
{
  #ifdef SUBGRAPH_MATCHING
  if(cntEmbedding * resultSoFar > nRemainingMatch) return; // nRest: the number of matches to find
  #else
  if(cntEmbedding * resultSoFar > 0) return;
  #endif
  //Step One: find a candidate that is with the max number of NEC nodes' candidate set containing it
  //locate the candidate with the last coverage number
  int maxCand = -1, maxCandSize = 1, maxCandIdx = -1;
  int * localSumNECCand = sumNECCand[labelID];
  int & localSizeSum = nSumNECCand[labelID];
  for(int i = 1; i < localSizeSum; i++)
  {
    int can = localSumNECCand[i];
    int temp = vCandIndex[can];
    if(mapTo[can] < 0 && temp > maxCandSize)
    {
      maxCandSize = temp;
      maxCandIdx = i;
    }
  }
  if(maxCandIdx == -1)
  {
    //in this case, there is no such a cand that has more than one query node having it as a candidate
    //hence, we can perform the combination and permutation to get the result here
    double res = 1;
    for(int i = 0; i < NECCountSetSize; i++)
    {
      if(NECCountSet[i] != localFlagQuery[i])
      {
        //this node is unmatched or still has unmatched NEC nodes
        int rest = NECCountSet[i] - localFlagQuery[i];
        int remainingNECRegionSize = NECRegionSize[i];
        if(remainingNECRegionSize < rest) //smaller ==> no result
          return;
        if(remainingNECRegionSize == rest) //equal ==> only one result
          res *= 1;
        else 
        {	
          //larger than rest, then perform the combination
          for(int x = remainingNECRegionSize - rest + 1; x <= remainingNECRegionSize; x++)
            res *= x;
          res /= factorization(rest);
        }
      }
    }
    res *= PRE_COMPUTED_PERMUTATION;
    cntEmbedding += res;
    return;
  } //end if maxCand = -1
  maxCand = localSumNECCand[maxCandIdx];//the case that maxCand is not -1
  swapValue(localSumNECCand[localSizeSum - 1], localSumNECCand[maxCandIdx]);//swap to the end
  localSizeSum--;  //reduce the size by one

  //then we examine the nec nodes that are covered by the max cand
  /*
   * The pruning rule here is that :
   *   for this cand, if there exists multiple query nodes covered by it and those nodes' cand size (the conditional cand size for the nec node) is 1,
   *   then there is no result. Because, one cand cannot be mapped to multiple query node at the same time.
   *   But if there only exists one such node, then it is OK to precede and this cand should be mapped to this query node.
   */

  //now we try to map this cand node, and remove it from all candidate regions that contain it
  //this for loop remove the maxCand from the candidate sets that contain it
  for (int i = 0; i < vCandIndex[maxCand]; i++)
  {
    //for each query node whose cand set containing it
    //qNodeIndex is the position of this maxCand in the nec node set
    int qNodeIndex = vCandInfo[maxCand][i].first;
    int posInRegion = vCandInfo[maxCand][i].second;
    int regionSize = NECRegionSize[qNodeIndex];
    if (posInRegion ==  regionSize - 1)
      NECRegionSize[qNodeIndex] --;
    else if(posInRegion < regionSize - 1){
      pair<int, int> & tempP = uCandInfo[qNodeIndex * maxNumCandidate + regionSize - 1];
      int lastCand = tempP.first;
      int lastPos = tempP.second;
      swapValue(tempP, uCandInfo[qNodeIndex * maxNumCandidate + posInRegion]);
      vCandInfo[maxCand][i].second = regionSize - 1;
      vCandInfo[lastCand][lastPos].second = posInRegion ;
      NECRegionSize[qNodeIndex] --;
    }
  } 

  //======multiple=====
  //sort accord to size of each query node's effective candidate size
  sort(vCandInfo[maxCand].begin(), vCandInfo[maxCand].begin() + vCandIndex[maxCand], sortByEffectiveCandSize);
  //the following loop maintains the position info after the sorting
  for (int i = 0; i < vCandIndex[maxCand]; i++)
  { 
    //for each query node containing maxCand
    pair<int, int> tempP =  vCandInfo[maxCand][i];
    int qNodeIndex = tempP.first;
    int posInRegion =  tempP.second;
    //set the position of "qNodeIndex" in the u_cand set of the candidate "maxCand"
    uCandInfo[qNodeIndex * maxNumCandidate + posInRegion].second = i;
  }
  for (int i = 0; i < vCandIndex[maxCand]; i++)
  { 
    //for each query node containing maxCand
    int qNodeIndex = vCandInfo[maxCand][i].first;
    if ((int)localFlagQuery[qNodeIndex] == NECCountSet[qNodeIndex])//check whether a nec node has been fully mapped or not
      continue;
    mapTo[maxCand] = nQueryVertex;
    localFlagQuery[qNodeIndex]++; //increase the flag value accordingly
    //here, if the query node "qNodeIndex" has been fully mapped, then we need to delete it from all the data node's candidate set which
    //containing it, such  that this query node will not be counted in the future round of finding the next maxCand
    if (localFlagQuery[qNodeIndex] == NECCountSet[qNodeIndex])
    {
      //this query node is fully mapped
      for (int j = 0; j < NECRegionSize[qNodeIndex]; j++)
      {
        //for each of its candidate, get the candidate from position "j"
        int vNode = uCandInfo[qNodeIndex * maxNumCandidate + j].first; 
        vector<pair<int, int>> & vcInfo = vCandInfo[vNode];
        int pos = uCandInfo[qNodeIndex * maxNumCandidate + j].second;
        if (pos == vCandIndex[vNode] - 1)
          vCandIndex[vNode]--;
        else if (pos < vCandIndex[vNode] - 1)
        {
          int size = vCandIndex[vNode];
          swapValue(uCandInfo[vcInfo [size - 1].first * maxNumCandidate +  vcInfo [size - 1].second],
          uCandInfo[vcInfo[pos].first * maxNumCandidate + vcInfo[pos].second]);	
          swapValue(vcInfo[size - 1], vcInfo[pos]);	
          vCandIndex[ vNode ]--;
        }
      }// end for j
    }// end if
    combine(cntEmbedding, resultSoFar); //recursive function
    #ifdef SUBGRAPH_MATCHING
    if(cntEmbedding * resultSoFar > nRemainingMatch)
    #else
    if(cntEmbedding * resultSoFar > 0)
    #endif
    {
      //clean up and return
      mapTo[maxCand] = -1;
      return;
    }
    //now we have to update the c(v) again, by putting back the last mapped query node (actually here: the qNodeIndex)
    if (localFlagQuery[qNodeIndex] == NECCountSet[qNodeIndex]){//NEC equal
      for (int j = 0; j < NECRegionSize[qNodeIndex]; j++)
      {
        int vNode = uCandInfo[qNodeIndex * maxNumCandidate + j].first;
        if (vNode == maxCand) //skip "maxCand"
          continue;
        vCandIndex[vNode]++; //extend the size by one
      }
    }
    localFlagQuery[qNodeIndex]--;//reduce this query node's flag, because it has been unmapped.
  }
  mapTo[maxCand] = -1;
  combine(cntEmbedding, resultSoFar); //recursive function
  #ifdef SUBGRAPH_MATCHING
  if(cntEmbedding * resultSoFar > nRemainingMatch) 
  #else
  if(cntEmbedding * resultSoFar > 0) 
  #endif
  {
    //clean up and return
    mapTo[maxCand] = -1;
    return;
  }
  // end =============== (multiple)

  //put back the max cand into each query node's candidate set
  for (int i = 0; i < vCandIndex[maxCand]; i++)
  {
    int qNodeIndex = vCandInfo[maxCand][i].first;
    NECRegionSize[qNodeIndex] ++;
  }
  localSizeSum ++;
  mapTo[maxCand] = -1;
}

inline double CountLeafMatch(const Graph& query, double matchFound) 
{
  #ifdef SUBGRAPH_MATCHING
  nRemainingMatch = nMaxMatch - matchFound;
  #endif
  double matchResult = 1;
  int NECSetSize = NECElementOffset.size() - 1;// -1 is because the last currE is actually redundant
  int can, uniqueUCandCounter, sumNEC; // the sum number of all nec node of a label
  int uCandCounter, sumCand;
  memset(necMapping, -1, sizeof (int) * (nQueryVertex + 1)); //for MM
  leafCandIndex = 0;
  //=== First scan === fastly identify whether all candidates satisfy the (cand < nec) and (sumCand < sum_nec). (the sum cands here are not unique cands)
  for(int i = 0; i < NECSetSize; i++) 
  { 
    //for each label
    int c, d, loc;
    int start = NECElementOffset[i].second; // start position
    int end = NECElementOffset[i + 1].second; // end position
    sumNEC = 0;
    sumCand = 0;
    for(int j = start; j < end; j++) 
    {
      //basic info about this nec node
      int parentID = query.NECElement[j].parentID;
      int NECNum = query.NECElement[j].sum;
      int uR = query.NECElement[j].representative;
      sumNEC += NECNum;
      uCandCounter = 0; //record the number of candidate for this nec node
      int parentPos = candPos[parentID];
      CandidateSpace & unit = candSpace[uR];
      #ifdef PRUNING_BY_EQUIVALENCE_SETS
      order[uR] = query.nNonLeaf;
      GetCell(query, uR, order[uR], unit.adjacent[0][parentPos], unit.nAdjacent[0][parentPos]);
      #endif
      for(int it = 0; it < unit.nAdjacent[0][parentPos]; it++)
      {
        int idx = unit.adjacent[0][parentPos][it];
        int v = unit.candidates[idx];
        if(mapTo[v] == uR) mapTo[v] = -1;
        if(mapTo[v] < 0) 
        {
          ++uCandCounter;
          #ifdef PRUNING_BY_EQUIVALENCE_SETS
          #ifdef SUBGRAPH_MATCHING
          c = unit.cell[idx];
          if(c > 0)
          {
            for(int x = 0; x < nActive[v]; ++x)
            {
              d = active[v][x].first;
              element[d].nonConflictCell[c >> 6] |= (1ULL << (c & 0x3f));
            }
          }
          else
          {
            for(int x = 0; x < nActive[v]; ++x)
            {
              d = active[v][x].first;
              loc = active[v][x].second;
              element[d].neverVisited[loc] = false;
              #ifdef MAPPING_FUNCTION_LOG
              cout<<"[L] Set (d: "<<d<<", u"<<element[d].vertex<<", v"<<v<<") as visited"<<endl;
              #endif
            }
          }
          #endif
          #endif
        }
        #ifdef PRUNING_BY_EQUIVALENCE_SETS
        else
        {
          Stack& ancElem = element[order[mapTo[v]]];
          c = unit.cell[idx];
          ancElem.conflictCell[c >> 6] |= (1ULL << (c & 0x3f));
        }
        #endif
      }
      if(uCandCounter < NECNum) return 0;
      sumCand += uCandCounter;
    }
    if(sumCand < sumNEC) return 0;
  }
  //=== Second scan: extract cands and perform the maximum matching
  for(int i = 0; i < NECSetSize; i++) 
  { //for each label
    int start = NECElementOffset[i].second; // start position
    int end = NECElementOffset[i + 1].second; // end position
    int label = NECElementOffset[i].first; //the label of this nec set
    //initialization
    posSumNECCand = 1;
    //flagSumNECCand must be all -1 here.  	
    uniqueUCandCounter = 0;
    sumNEC = 0;
    sumCand = 0;
    nodeIndex = 1; //for MM
    int * localSumNECCand = sumNECCand[label];
    for(int j = start; j < end; j++) 
    {
      //basic info about this nec node
      int parentID = query.NECElement[j].parentID;
      int NECNum = query.NECElement[j].sum;
      int uR = query.NECElement[j].representative;
      sumNEC += NECNum;
      for(int mm = 0; mm < NECNum; mm++) //for MM
      necMapping[nodeIndex ++] = j;
      uCandCounter = 0; //record the number of candidate for this nec node
      int parentPos = candPos[parentID];
      CandidateSpace & unit = candSpace[uR];
      int startPos = leafCandIndex;
      for(int it = 0; it < unit.nAdjacent[0][parentPos]; it++)
      {
        int can = unit.candidates[unit.adjacent[0][parentPos][it]];
        if(mapTo[can] < 0)
        {
          if(flagSumNECCand[can] == -1)
          { //if not already in the candidate set,
            flagSumNECCand[can] = posSumNECCand; //and mark it true and indicate its position in "sumNECCands"
            localSumNECCand[posSumNECCand ++] = can; //then add it into the candidate's "posSumNECCand" position,
            uniqueUCandCounter++;
          }
          uCand[j * maxNumCandidate + uCandCounter] = flagSumNECCand[can]; //for the maximum matching
          uCandCounter++;
          leafCand[leafCandIndex++] = can; //store the candidate for the second stage (combination)
        }
      }    
      leafCandInfo[j] = make_pair(startPos, uCandCounter);
      NECRegionSize[j] = uCandCounter;  //set the size of the j-th candidate set of this nec node
      sumCand += uCandCounter;
      nSumNECCand[label] = posSumNECCand; //added on 20:15 09/05/2016
    }
    if(uniqueUCandCounter < sumNEC)
    {
      //no result //reset "flagSumNECCand"
      for(int i = 1; i < posSumNECCand; i++) //it starts from 1 NOT 0
        flagSumNECCand[ localSumNECCand[i] ] = -1;
      return 0;
    }
    //after the  "sumNECCands" is initialized, we reset the "flagSumNECCand" for next time use
    for(int i = 1; i < posSumNECCand; i++) //it starts from 1 NOT 0
      flagSumNECCand[localSumNECCand[i]] = -1;
    //== BEGIN checking MAXIMUM MATCHING ==============================================
    //using nec_region as the adjacent list
    memset(pairU, NULL, sizeof (int) * nQueryVertex); //must reset before using
	  #ifdef SUBGRAPH_MATCHING
    memset(pairV, NULL, sizeof (int) * maxNumCandidate); //must reset before using
	  #else
    memset(pairV, NULL, sizeof (int) * currG->nVertex); //must reset before using
	  #endif
    int matching = 0;
    while(BFS_MM())
    {
      for(int u = 1; u < nodeIndex; u++)
        if(pairU[u] == NULL)
          if(DFS_MM(u))
            matching++;
    }
    if(matching != nodeIndex - 1)
    {
      return 0;
    }
    //== END checking MAXIMUM MATCHING ==========================================
    NECRank[i] = make_pair(i, make_pair(sumNEC, sumCand));
  }
  sort(NECRank, NECRank + NECSetSize, sortLLCNode);
  //Deal with the labels one by one
  for(int i = 0; i < NECSetSize; i++) 
  { //for each nec node label
    #ifdef SUBGRAPH_MATCHING
    if(matchResult > nRemainingMatch) continue;
    #else
    if(matchResult > 0) continue; //20190425. Subgraph query
    #endif
    double cntLocalMatchLabel = 0;
    int posElem = NECRank[i].first;  // the label index of this nec label set ==> this is not the actual label!!!!
    int nodeSum = NECRank[i].second.first; // the number of node in this nec label set
    labelID = NECElementOffset[posElem].first;
    if(nodeSum == 1) 
    { //==== CASE ONE : there is only one node in this nec set with this label =====
      cntLocalMatchLabel = (long long)NECRank[i].second.second; //we have computed this one in the last step
      if(cntLocalMatchLabel == 0) return 0;
    } 
    else
    { //==== CASE TWO : more than one node, and possible more than one start nodes (nec_units)
      memset(localFlagQuery, 0, sizeof(char) * nQueryVertex);
      posSumNECCand = 1;//omit zero, start from 1 ===> Leave 0 for NULL of Maximum Matching
      int start = NECElementOffset[posElem].second;
      int end = NECElementOffset[posElem + 1].second;
      int localNECSize = end - start; //number of nec with posElem
      //======= sorting here is to decide the processing sequence of the nec nodes in this nec label set
      for (int j = start, x = 0; j < end; j++, x++)
  	    vNECCount[x] = make_pair(j, query.NECElement[j].sum);
      sort(vNECCount, vNECCount + localNECSize, sortBySecondElement);
      int * toCleanVCandIndex = arrayToClean; //store the candidates that need to be cleaned for the array "vCandIndex"
      NECCountSetSize = localNECSize;
      int idxToCleanVCandIndex = 0;
      //in this loop, for this point beyond, before any return, "vCandIndex" must be cleaned up using the above two parameters
      for(int j = 0; j < localNECSize; j++) 
      {
        int necIndex = vNECCount[j].first;
        int necCount = vNECCount[j].second;
        NECCountSet[j] = necCount; // the sum of nodes that the nec representative stands for
        uCandCounter = 0;
        pair<int, int> tempP = leafCandInfo[necIndex];
        for(int x = tempP.first; x < tempP.first + tempP.second; x++)
        {
          can = leafCand[x];
          int temp = vCandIndex[can]; //get the position for storing this candidate
          //the temp-th position of vCandInfo[can] stores a candidate "can" which is for "j-th" nec node and
          //stored in the "uCandCounter-th" position in the candidate set
          vCandInfo[can][temp] = make_pair(j, uCandCounter); // pair<j-th nec node, can's pos in the candidate set of j-th nec node>
          //store the candidate and its corresponding info ==>  <the candidate can, the candidate's pos in vCandIndex[can]>
          //in the single-array: "j-th" nec node and the "uCandCounter-th" position.
          uCandInfo[j * maxNumCandidate + uCandCounter] = make_pair(can, temp);
          if(temp == 0)//record for the cleaning purpose
            toCleanVCandIndex[ idxToCleanVCandIndex ++] = can;
          vCandIndex[can]++; // update the counter(next available position) for "vCandIndex[can]"
          uCandCounter++; //update the candidate counter(next available position)
        }
        //if the number of candidate for this nec node, is smaller than the nec count of this nec node, then no result can be found.
        NECRegionSize[j] = uCandCounter;
      }//end for j

      //to reduce the computation cost, we pre-compute the value of the all factorization
      PRE_COMPUTED_PERMUTATION = 1;
      for(int i = 0; i < NECCountSetSize; i++)
  	    PRE_COMPUTED_PERMUTATION *= factorization(NECCountSet[i]);
      int * mappedDataVertex = flagSumNECCand;
      int posMappedDataVertex = 0;
      //in this loop, for this point beyond, before any return, "mapTo" must be cleaned up using the above two parameters.
      int found = 1;
      while (found > 0)
      {
        //=============preprocess before the combination
        found = 0;
        //first, we scan each nec node, to see whether there exists such a node that its candidate size is equal to its unmatched nec size
        for (int i = 0; i < NECCountSetSize; i++)
        {
	        if(NECRegionSize[i] != 0 && NECCountSet[i] > NECRegionSize[i])
          {
            while(posMappedDataVertex != 0)
	          {
              mapTo[mappedDataVertex[--posMappedDataVertex]] = -1;
	            mappedDataVertex[posMappedDataVertex] = -1;
            }
            while(idxToCleanVCandIndex != 0) //cleaning up vCandIndex
              vCandIndex[toCleanVCandIndex[--idxToCleanVCandIndex]] = 0;
            return 0;
          }
          //all can be matched one on one
          if(NECCountSet[i] == NECRegionSize[i])
          {
            found++;
            localFlagQuery[i] = NECCountSet[i];//mark this nec node all matched
            //update the c(v); there is no neccessary to update the c(v)
            for (int j = 0; j < NECRegionSize[i]; j++)
            {
              int vNode = uCandInfo[i * maxNumCandidate + j].first; //the j-th candidate
              mapTo[vNode] = nQueryVertex; //set this cand's flag to 1
              mappedDataVertex[posMappedDataVertex++] = vNode; //recording this cand
              //for each query node that contains this candidate "vNode"
              for(int x = 0; x < vCandIndex[vNode]; x++)
              { // x is the position in "vCandInfo[vNode]"
                int qNodeIndex = vCandInfo[vNode][x].first;//get the original index in the nec node set
                if(qNodeIndex == i) continue;
  	              //the position of vNode in the corresponding("qNodeIndex-th") candidate set.
  	            int posInRegion = vCandInfo[vNode][x].second;
                int regionSize = NECRegionSize[qNodeIndex];
                if (posInRegion == regionSize - 1) //this is the last node in the candidate set
  	            NECRegionSize[qNodeIndex]--; //remove this candidate by downsizing the NECRegionSize by one
  	            else if(posInRegion < regionSize - 1)
                {
                  //normal case: not the last one.
                  //swap the two candidates and their corresponding info, then downsize the NECRegionSize by one
                  pair<int, int> & tempP = uCandInfo[qNodeIndex * maxNumCandidate +  regionSize - 1];
                  int lastCand = tempP.first;
                  int lastPos = tempP.second;
                  swapValue ( tempP, uCandInfo [ qNodeIndex * maxNumCandidate + posInRegion ]); //swap the candidate
                  vCandInfo[vNode][x].second = regionSize - 1;
                  vCandInfo[lastCand][lastPos].second = posInRegion;
                  NECRegionSize[qNodeIndex]--;
  	            }
                if(NECRegionSize[qNodeIndex] < NECCountSet[qNodeIndex])
                {
                  while(posMappedDataVertex != 0)
                  {
                    mapTo[mappedDataVertex[--posMappedDataVertex]] = -1;
                    mappedDataVertex[posMappedDataVertex] = -1;
                  }
                  while(idxToCleanVCandIndex != 0) //cleaning up vCandIndex
                  vCandIndex[toCleanVCandIndex[--idxToCleanVCandIndex]] = 0;
                  return 0;
                } 
              } // end for x 
            }// end for j 
            NECRegionSize[i] = 0;
          }// end if
        }// end for i 
      }//end PREPROCESSING
      //================== END PREPROCESSING ==============================
      combine(cntLocalMatchLabel, matchResult);
      
      // Clean up
      while(posMappedDataVertex != 0)
      {
        mapTo[mappedDataVertex[--posMappedDataVertex]] = -1;
	      mappedDataVertex[posMappedDataVertex] = -1;
      }
      while(idxToCleanVCandIndex != 0) //cleaning up vCandIndex
        vCandIndex[ toCleanVCandIndex[--idxToCleanVCandIndex]] = 0;
      if(cntLocalMatchLabel == 0) return 0;
    } //end else : end of case TWO
    matchResult *= cntLocalMatchLabel;
  }
  #ifdef MAPPING_FUNCTION_LOG
  printf("[L] Found %lf embeddings.\n", matchResult);
  #endif
  return matchResult;
}
#endif