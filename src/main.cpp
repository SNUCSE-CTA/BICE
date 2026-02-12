#include "daf.h"
#include "bice.h"

void writeAnswer(string& answerFile)
{
  ofstream outFile(answerFile);
  for(int g = 0; g < nGraph; ++g){
    if(answer[g]) outFile<<g<<endl;
  }
  outFile.close();
}

inline bool BuildCS()
{
  if(!FilterByCount(*currQ, *currG)) return false;
  root = SelectRootVertex(*currQ, *currG);
  ExtractResidualStructures(*currQ);
  BuildDAG(*currQ);
  if(!TopDownInitial(*currQ, *currG)) return false;
  //if(!BottomUp(*currQ, *currG)) return false;

  #ifdef SUBGRAPH_MATCHING
  // Neighbor-safety in the 2nd refinement for subgraph matching
  if(!filteringWithDAG(*currQ, *currG, bottomUp, childNgb, true)) return false;
  #else
  // Neighbor-safety in the 2nd refinement for subgraph search
  if(!filteringWithDAG(*currQ, *currG, bottomUp, childNgb, false)) return false;
  #endif   
  // Neighbor-safety in the 3rd refinement
  if(!filteringWithDAG(*currQ, *currG, topDown, parentNgb, true)) return false;

  if(!ConstructAdjacencyList(*currQ, *currG, true)) return false;

  bice::eqs::makeEqSet(*currQ, *currG);

  return true;
}

inline void ProcessQuery()
{
  bool isCandidate;
  chrono::high_resolution_clock::time_point startClock, endClock;
  for(int x = 0; x < queryGraph.size(); ++x)
  {
    exceedTimeLimit = false;
    filteringTime = 0;
    verificationTime = 0;
    nMatch = 0;
    nCandidate = 0;
    recursiveCallCount = 0;
    currQ = queryGraph[x];
    AllocateForQueryGraph(*currQ);
    GetQueryDataStructure(*currQ);
    isTree = true;
    for(int i = 0; i < nQueryVertex; i ++){
      if (currQ->core[i] >= 2){
        isTree = false;
    	break;
      }
    }
    #ifdef MULTIPLE_QUERIES
    cout<<"*Processing q"<<x<<"(|V|: "<<currQ->nVertex<<")"<<endl;
    #endif
    for(int y = 0; y < dataGraph.size(); ++y)
    {
      printf("data_id : %d\n", y);
      SetQueryGraphResource(*currQ);
      currG = dataGraph[y];
      #ifdef MAPPING_FUNCTION_LOG
	    cout<<"Processing g"<<y<<"(|V|: "<<currG->nVertex<<")"<<endl;
      #endif
      startPoint = GetClock();

      startClock = GetClock();
      isCandidate = BuildCS();
      
      endClock = GetClock(); 
      filteringTime += TimeDiffInMilliseconds(startClock, endClock);
      if(!isCandidate) continue;
      
      bice::BipartiteMatching::allocate(candSpace,uSequenceSize,maxNumDataVertex,*currQ);
      
      ++nCandidate;
      #ifdef SUBGRAPH_MATCHING
      ComputeAuxiliaryDataStructureSize();
      cout<<"Sum of |C(u)|: "<<auxDataStructureSize<<endl;
      #endif
      #ifdef FILTERING_TEST
      #else
      startClock = GetClock();
      if(usePathSizeOrder) 
        InitializeWeight(*currQ); 

      #ifdef SUBGRAPH_MATCHING
      recursiveCallCountPerGraph = bice::backtrack(*currQ, *currG, y, 100000);
      recursiveCallCount += recursiveCallCountPerGraph;
      #else
      recursiveCallCountPerGraph = bice::backtrack(*currQ, *currG, y);
      recursiveCallCount += recursiveCallCountPerGraph;
      #endif
      
      #ifdef DEBUG
      if(answer[y]){
        puts("found");
      }
      else{
        puts("notfound");
      }
      #endif
      endClock = GetClock();
      verificationTime += TimeDiffInMilliseconds(startClock, endClock);
      #endif
      #ifdef SPLIT_RECURSIVE_CALLS
      cout<<"Recursive calls(g"<<y<<"): "<<recursiveCallCountPerGraph<<endl;
      #endif
      
      bice::BipartiteMatching::deallocate(uSequenceSize);
    }
    #ifdef MULTIPLE_QUERIES
    cout<<"Exceeded Time Limit: "<<(exceedTimeLimit?1:0)<<endl;
    Deallocate(*currQ);
    #endif
    totalTime = filteringTime + verificationTime;
    cout<<"Total Recursive Call Count: "<<recursiveCallCount<<endl;
    #ifdef SUBGRAPH_MATCHING
    cout<<"Number of Matches: "<<doubleToString(nMatch)<<endl;
    #else
    cout<<"Number of Candidate Graphs: "<<nCandidate<<endl;
    cout<<"Number of Answer Graphs: "<<nMatch<<endl;
    writeAnswer(answerFile);
    #endif
    cout<<"Filtering Time (ms): "<< filteringTime<<endl;
    cout<<"Verification Time (ms): "<<verificationTime<<endl;
    cout<<"Processing Time (ms): "<<totalTime<<endl;
  }
}

// Subgraph query
// Usage: ./Program -d[dataFormat] [dataFile] -q[queryFormat] [queryFile] -o [answerFile]

// Subgraph matching
// Usage: ./Program -d[dataFormat] [dataFile] -q[queryFormat] [queryFile]

// dataFormat:
// -i: Igraph format
// -g: GFU format

// queryFormat:
// -c: CFL format
// -i: Igraph format
// -g: GFU format

int main(int argc, char *argv[])
{
  setvbuf(stdout, NULL, _IOLBF, 0);
  int queryFormat = GFU_FORMAT;
  int dataFormat = GFU_FORMAT;
  for(int i = 1; i < argc; ++i)
  {
    if(argv[i][0] == '-')
    {
      switch(argv[i][1])
      {
        case 'd':
	  switch(argv[i][2])
	  {
	    case 'i':
	      dataFormat = IGRAPH_FORMAT;
	      break;
	    case 'g':
	      dataFormat = GFU_FORMAT;
	      break; 
	  }
	  dataGraphFile = argv[i+1];
          break;
        case 'q':
	  switch(argv[i][2])
	  {
	    case 'c':
	      queryFormat = CFL_FORMAT;
	      break;
	    case 'i':
	      queryFormat = IGRAPH_FORMAT;
	      break;
	    case 'g':
	      queryFormat = GFU_FORMAT;
	      break; 
	  }
	  queryGraphFile = argv[i+1];
          break;
	case 'o':	
	  answerFile = argv[i+1];
          break;
      }
    }
  }
  cout<<"Data file: "<<dataGraphFile<<endl;
  cout<<"Query file: "<<queryGraphFile<<endl;
  cout<<"Output file: "<<answerFile<<endl;
  if(dataFormat == IGRAPH_FORMAT)
    ReadIgraphFormat(dataGraphFile, dataGraph);
  else
    ReadGFUFormat(dataGraphFile, dataGraph);
 
  ProcessDataGraphs();
  AllocateForDataGraph();
  
  if(queryFormat == CFL_FORMAT)
    ReadCFLFormat(queryGraphFile, queryGraph);
  else if(queryFormat == IGRAPH_FORMAT)
    ReadIgraphFormat(queryGraphFile, queryGraph);
  else
    ReadGFUFormat(queryGraphFile, queryGraph);

  ProcessQuery();
  return 0;
}
