#ifndef ADV_BSF
#define ADV_BSF

namespace adv_bsf{
  void pT(int n)
  {
    for(int i=0;i<n;i++){
      printf("\t");
    }
  }

  namespace eqs{
    using ii = pair < int, int >;
    const int DATAV = 4300000;
    const int LABEL = 50000;
    int candArray[DATAV];
    int posArray[DATAV];
    int nextPivot[DATAV];

    void pT(int n)
    {
      for(int i=0;i<n;i++){
        printf("\t");
      }
    }

    void printCandidates
    (
      const Graph &query,  const Graph & data,
      int u, int uCand, int uCandIdx
    )
    {
      CandidateSpace& currSet = candSpace[u];
      set < ii > st;
      for(int idxChQuVer = 0; idxChQuVer < DAG_child_query_size[u]; idxChQuVer++)
      {
        int chQuVer = DAG_child_query[u][idxChQuVer];
        if(query.NECMap[chQuVer] == -1 || query.NECMap[chQuVer] == chQuVer){
          CandidateSpace& childSet = candSpace[chQuVer];
          int uPos = DAG_child_query_parent_index[u][idxChQuVer];
          for(int i=0;i<childSet.nAdjacent[uPos][uCandIdx];i++){
            int idx = childSet.adjacent[uPos][uCandIdx][i];
            st.insert(ii(chQuVer,childSet.candidates[idx]));
          }
        }
      }

      for(int idxParQuVer = 0; idxParQuVer < DAG_parent_query_size[u]; idxParQuVer++)
      {
        int parQuVer = DAG_parent_query[u][idxParQuVer];
        int label_ParQuVer = query.label[parQuVer];
       
        CandidateSpace & parSet = candSpace[parQuVer];
        for(int posCand=0;posCand<parSet.size;posCand++)
        {
          int candidate = parSet.candidates[posCand];
          candToPos[candidate] = posCand;
        }
        
        auto iter = data.labelToNbrOffset.find(make_pair(uCand,label_ParQuVer));	
	      if(iter == data.labelToNbrOffset.end()) continue;
	      for(int w=iter->second.first;w<iter->second.first+iter->second.second;w++)
        {
          int nbr_data = data.nbr[w];
          if(candToPos[nbr_data] >= 0)
          {
            st.insert(ii(-parQuVer, nbr_data));
          }
        }

        for(int posCand=0;posCand<parSet.size;posCand++)
        {
          int candidate = parSet.candidates[posCand];
          candToPos[candidate] = -1;
        }
      }
      for(const ii & e:st){
        printf("[%d,%d]", e.first, e.second);
      }
      puts("");
    }

    void printEqSetResult
    (
      const Graph & query, const Graph & data
    )
    {
      for(int i=0;i<uSequenceSize;i++){
        //[0]
        int u = uSequence[i];
        CandidateSpace& currSet = candSpace[u];
        if(query.NECMap[u] == -1 || query.NECMap[u] == u){
          
          printf("eqSetAlign u : %d ( sz : %d) \n", u, currSet.size);
          for(int j=0;j<currSet.size;j++){
            int id = currSet.eqSetId[j];
            printf("id : %d\n", id);
            printCandidates(query, data, u, candSpace[u].candidates[j], j);
          }
          puts("eqSetIds");
          for(int j=0;j<currSet.size;j++)
          {
            printf("[%d]", currSet.eqSetId[j]);
          }
          puts("");
        }
      }
      puts("");
    }

    void makeEqSet
    (
      const Graph& query, const Graph & data
    )
    {
      #ifdef DEBUG
      printf("a0. makeEqSet is called\n");
      #endif

      for(int i=0;i<uSequenceSize;i++){
        //[0]
        int u = uSequence[i];

        CandidateSpace& currSet = candSpace[u];
        for(int j = 0; j < currSet.size; ++j) {
          candArray[j] = j;
        }
        //[1] initialization
        nextPivot[0] = currSet.size;

        //[2.1] make pivot, by childs set
        for(int idxChQuVer = 0; idxChQuVer < DAG_child_query_size[u]; idxChQuVer++){
          int chQuVer = DAG_child_query[u][idxChQuVer];
          if(query.NECMap[chQuVer] == -1 || query.NECMap[chQuVer] == chQuVer){
            //[2.1] child candidates examine
            CandidateSpace& childSet = candSpace[chQuVer];
            
            int uPos = DAG_child_query_parent_index[u][idxChQuVer];
            
            fill(posArray, posArray + currSet.size, 0);

            for(int idxChQuCand = 0; idxChQuCand < childSet.size; idxChQuCand++){
              int chQuCand = childSet.candidates[idxChQuCand];
              int idxStartCand = 0;

              while(idxStartCand < currSet.size){
                int head = idxStartCand;
                int tail = nextPivot[idxStartCand];

                while(head < tail){
                  bool exist = false;
                  while(posArray[candArray[head]] < childSet.nAdjacent[uPos][candArray[head]]){
                    if(childSet.adjacent[uPos][candArray[head]][posArray[candArray[head]]] < idxChQuCand){
                      posArray[candArray[head]]++;
                    }
                    else{
                      break;
                    }
                  }
                  if(
                    posArray[candArray[head]] < childSet.nAdjacent[uPos][candArray[head]] &&
                    childSet.adjacent[uPos][candArray[head]][posArray[candArray[head]]] == idxChQuCand)
                  {
                    ++head;
                  }
                  else{
                    --tail;
                    swap(candArray[head], candArray[tail]);
                  }
                }

                int idxNextStartCand = nextPivot[idxStartCand];

                if(head != nextPivot[idxStartCand] && tail != idxStartCand){
                  nextPivot[idxStartCand] = head;
                  nextPivot[head] = idxNextStartCand;
                }

                idxStartCand = idxNextStartCand;
              }
            }
          }
        }
        /*
        //[2.2] make pivot, by childs set
        for
        (
          int idx_parQuVer = 0; idx_parQuVer < DAG_parent_query_size[u]; idx_parQuVer++
        )
        {
          int parQuVer = DAG_parent_query[u][idx_parQuVer];
          CandidateSpace & parSet = candSpace[parQuVer];

          for
          (
            int idx_parCand = 0; idx_parCand < parSet.size; idx_parCand++
          )
          {
            int nVisited = 0;

            #ifdef DEBUG
            printf("u : %d uPar : %d parCand : %d\n", u, parQuVer, parSet.candidates[idx_parCand]);
            #endif

            for
            (
              int j = 0; j < currSet.nAdjacent[idx_parQuVer][idx_parCand]; j++  
            )
            {
              int pos = currSet.adjacent[idx_parQuVer][idx_parCand][j];
              arrayToClean[nVisited++] = pos;
              candToPos[pos] = j;
            }

            int idxStartCand = 0;
            while(idxStartCand < currSet.size)
            {
              int head = idxStartCand;
              int tail = nextPivot[idxStartCand];
              
              while(head < tail)
              {
                if(candToPos[candArray[head]]  < 0)
                {
                  ++head;
                }
                else
                {
                  --tail;
                  swap(candArray[head], candArray[tail]);
                }
              }

              int idxNextStartCand = nextPivot[idxStartCand];
              if(head != nextPivot[idxStartCand] && tail != idxStartCand){
                nextPivot[idxStartCand] = head;
                nextPivot[head] = idxNextStartCand;
              }
              idxStartCand = idxNextStartCand;
            }
            
            while(nVisited!=0)
            {
              candToPos[arrayToClean[--nVisited]] = -1;
            }
          }
        }
        */

        //[3] after pivotMaking, partition the cell
        int idxStartCand = 0;
        currSet.szEqSetId = 0;

        int idx_components_cell = 0;

        while(idxStartCand < currSet.size){
          currSet.szEqSetId++;
          const int rep = currSet.candidates[candArray[idxStartCand]];
          sort(candArray + idxStartCand, candArray + nextPivot[idxStartCand]);
          for(int idxCand = idxStartCand; idxCand < nextPivot[idxStartCand]; idxCand++){
            currSet.eqSetId[candArray[idxCand]] = currSet.szEqSetId;
            currSet.components_cell[idx_components_cell++] = candArray[idxCand];
          }
          currSet.intervals_cell[currSet.szEqSetId - 1] = make_pair(idxStartCand, nextPivot[idxStartCand]);
          idxStartCand = nextPivot[idxStartCand];
        }
        currSet = candSpace[u];
      }
    }
  }
  namespace Queue{
    inline void reinsertToQueue(int u, weight_type weight, int position)
    {
      extendable[posExtendable] = extendable[position];
      weightArray[posExtendable] = weightArray[position];  
      posExtendable++;
      extendable[position] = u;
      weightArray[position] = weight;
    }

    inline void insertToQueue(int u, weight_type weight)
    {
      extendable[posExtendable] = u;
      weightArray[posExtendable] = weight;
      posExtendable++;
    }

    inline void popFromQueue(int& current, weight_type& currW, int& currP)
    {
      int optPos = -1;
      weight_type optWeight = WEIGHT_MAX;
      for(int i = 0; i < posExtendable; i++){
        if(weightArray[i] < optWeight)
        {
          optPos = i;
          optWeight = weightArray[i];
        }
        else if(weightArray[i] == optWeight){
          if(extendable[optPos] > extendable[i]){
            optPos = i;
          }
        }
      }
      current = extendable[optPos];
      currW = weightArray[optPos];
      currP = optPos;
      extendable[optPos] = extendable[posExtendable - 1];
      weightArray[optPos] = weightArray[posExtendable - 1];
      posExtendable--;
    }

    inline void removeFromQueue(int n)
    {
      posExtendable -= n;
    }

    inline void clearQueue()
    {
      posExtendable = 0;
    }

    inline bool isEmpty()
    {
      return posExtendable == 0;
    }
  }

  namespace Extendable{
    long long intersectTwoExtendable(
      const int* A, const int sizeA,
      const int* B, const int sizeB,
      int* C, int &sizeC
    ){
      sizeC = 0;
      int idxA = 0, idxB = 0;
      while(idxA < sizeA and idxB < sizeB)
      {
        const int candIdA = A[idxA];
        const int candIdB = B[idxB];
        if(candIdA == candIdB)
        {
          C[sizeC++] = candIdA;
          ++idxA;
          ++idxB;
        }
        else if(candIdA < candIdB)
          ++idxA;
        else
          ++idxB;
      }
      return sizeC;
    }

    long long getExtendableCandidates
    (
      int u, int uIdx,
      int up, int candParentIdx,
      int numParentsMapped
    )
    {
      CandidateSpace & unit = candSpace[u];
      int iecPos = numParentsMapped - 1;
      int uParIdx = DAG_child_query_parent_index[up][uIdx];
      if(numParentsMapped == 1)
      {
        iecSize[u][0] = 0;
        for(int i=0;i<unit.nAdjacent[uParIdx][candParentIdx];i++){
          int cand = unit.adjacent[uParIdx][candParentIdx][i];
          iec[u][0][iecSize[u][0]++] = cand;
        }
      }
      else
      {
        int iecPos = numParentsMapped - 1; 
        iecSize[u][iecPos] = 0;
        int indexSize = iecSize[u][iecPos-1];
        int newIndexSize = 
        unit.nAdjacent[uParIdx][candParentIdx];
        intersectTwoExtendable(
          iec[u][iecPos-1], indexSize,
          unit.adjacent[uParIdx][candParentIdx], 
          newIndexSize,
          iec[u][iecPos], iecSize[u][iecPos]
        );
      }
      return (long long)iecSize[u][iecPos];
    }
  }

  inline void initFs(uint64_t* fs)
  {
    fill(fs, fs + FAILING_SET_SIZE, (uint64_t)(0));
  }

  bool checkRedundant(const uint64_t* f1, int u)
  {
    int uQuo = u >> 6;
    return (f1[uQuo] & (1ULL << (u & 0x3f))) == (uint64_t)(0);
  }

  void copyFs(const uint64_t* srcFs, uint64_t* destFs)
  {
    for(int x=0;x<FAILING_SET_SIZE;x++)
    {
      destFs[x] = srcFs[x];
    }
  }

  void unionFs(const uint64_t* srcFs, uint64_t* destFs)
  {
    for(int x=0;x<FAILING_SET_SIZE;x++)
    {
      destFs[x] |= srcFs[x];
    }
  }

  void addBitFs(uint64_t* destFs, int u)
  {
    int uQuo = u >> 6;
    destFs[uQuo] |= (1ULL << (u & 0x3f));
  }

  void makeFull(uint64_t* destFs)
  {
    for(int x=0;x<FAILING_SET_SIZE;x++)
    {
      destFs[x] = 0xffffffffffffffff;
    }
  }

  void printFs(const uint64_t *fs)
  {
    int u = 0;

    for(int x=0;x<FAILING_SET_SIZE;x++)
    {
      for(int y=0;y<64;y++)
      {
        uint64_t k = (uint64_t)(1) << y;
        if(fs[x] & k){
          printf("[%d]", u);
        }
        u++;
      }
    }
    puts("");
  }

  namespace BipartiteMatching
  { 
    struct Edge
    {
      int head, tail, inverse_index;
      int capacity, flow;
      int candidateIdx;
      Edge(int h = 0, int t = 0,int c = 0,int ii = 0):
        head(h),tail(t),capacity(c),flow(0),
        inverse_index(ii){}
      int residual() const { return capacity - flow; }
      void update(int h, int t, int c, int ii){
        this->head = h;
        this->tail = t;
        this->capacity = c;
        this->flow = 0;
        this->inverse_index = ii;
      }
      void flip(){
        if(this->capacity == 0){
          this->flow = (-1) - (this->flow);
        }
        else{
          this->flow = 1 - (this->flow);
        }
      }
      bool backWardEdge() const{
        return this->capacity == 0;
      }
      void print(){
        printf("from %d to %d cap : %d flow : %d ii : %d\n", 
          this->head, this->tail, this->capacity, this->flow,
          this->inverse_index
        );
      }
    };

    int numLeftVer = 0;

    int *numLeftVerAdj;
    Edge **leftVerAdj;

    int* repVerLeft;

    int* matchedRightFromLeft;
    Edge **matchedEdgeFromLeft;

    bool* isLeafLeft;
    int* currEqSetId;
    uint64_t* fsBiMat;

    int* visRightVerAdj;
    int* numRightVerAdj;
    Edge **rightVerAdj;
    
    int* matchedLeftFromRight;
    Edge** matchedEdgeFromRight;
    
    int idx_right_adj = 0;

    int pointer = 0;

    long long *visLeftVerSearch;
    long long *visRightVerSearch;
    long long idxSearch = 0;

    const int DATAV = 4300000;
    int cntCandFreq[DATAV];
    int visCantFreq[DATAV] = { 0, };
    int idxCandFreq = 0;

    int ** components_cell_visited_per_query_ver;
    int * size_components_cell_visited_per_query_ver;

    void allocate(
      const CandidateSpace candSpace[MAX_NUM_VERTEX],
      const int uSequenceSize,
      const int maxNumDataVertex,
      const Graph& query
    )
    {
      #ifdef DEBUG
      printf("b0. allocate is called\n");
      #endif
      ++idxCandFreq;

      numLeftVer = 0;
      numLeftVerAdj = new int[uSequenceSize]();

      leftVerAdj = new Edge*[uSequenceSize](); 
      for(int i=0;i<uSequenceSize;i++)
      {
        int rn = i;
        if(
          query.NECMap[i] != -1 &&
          query.NECMap[i] != i
        ){
          rn = query.NECMap[i];
        }
        for(int c=0;c<candSpace[rn].size;c++)
        {
          int cand = candSpace[rn].candidates[c];
          if(visCantFreq[cand] < idxCandFreq){
            visCantFreq[cand] = idxCandFreq;
            cntCandFreq[cand] = 1;
          }
          else{
            cntCandFreq[cand]++;
          }
        }
        leftVerAdj[i] = new Edge[2 * candSpace[rn].size]();
      }

      repVerLeft = new int[uSequenceSize];
      matchedRightFromLeft = new int[uSequenceSize]();
      matchedEdgeFromLeft = new Edge*[uSequenceSize]();

      isLeafLeft = new bool[uSequenceSize]();
      currEqSetId = new int[uSequenceSize]();
      fsBiMat = new uint64_t[uSequenceSize]();

      visRightVerAdj = new int[maxNumDataVertex]();
      numRightVerAdj = new int[maxNumDataVertex]();
      
      rightVerAdj = new Edge*[maxNumDataVertex]();
      for(int r=0;r<maxNumDataVertex;r++)
      {
        if(cntCandFreq[r] > 0 && visCantFreq[r] == idxCandFreq)
          rightVerAdj[r] = new Edge[2 * cntCandFreq[r]]();
      }
      
      matchedLeftFromRight = new int[maxNumDataVertex]();
      matchedEdgeFromRight = new Edge*[maxNumDataVertex]();
      
      idx_right_adj = 0;
      pointer = 0;

      visLeftVerSearch = new long long[uSequenceSize]();
      visRightVerSearch = new long long[maxNumDataVertex]();
      idxSearch = 0;

      components_cell_visited_per_query_ver = new int * [uSequenceSize]();
      for(int i=0;i<uSequenceSize;i++)
      {
         components_cell_visited_per_query_ver[i] = new int[maxNumCandidate]();  
      }
      size_components_cell_visited_per_query_ver = new int[uSequenceSize]();
    }
    
    void deallocate(int uSequenceSize)
    {
      #ifdef DEBUG
      printf("b1. deallocate is called\n");
      #endif
      delete[] numLeftVerAdj;
      for(int i=0;i<uSequenceSize;i++)
      {
        delete[] leftVerAdj[i];
      }
      delete[] leftVerAdj;
      delete[] repVerLeft;
      delete[] matchedRightFromLeft;
      delete[] matchedEdgeFromLeft;

      delete[] isLeafLeft;
      delete[] currEqSetId;
      delete[] fsBiMat;
      delete[] visRightVerAdj;
      delete[] numRightVerAdj;
      for(int r=0;r<maxNumDataVertex;r++)
      {
        if(cntCandFreq[r] > 0 && visCantFreq[r] == idxCandFreq)
          delete[] rightVerAdj[r];
      }
      delete[] rightVerAdj;
      delete[] matchedLeftFromRight;
      delete[] matchedEdgeFromRight;
      delete[] visLeftVerSearch;
      delete[] visRightVerSearch;
      
      for(int i=0;i<uSequenceSize;i++)
      {
        delete[] components_cell_visited_per_query_ver[i];
      }
      delete[] components_cell_visited_per_query_ver;
      delete[] size_components_cell_visited_per_query_ver;
    }

    bool findAugmentedPath(int currVer, bool isLeft, int length = 0,
      bool printFlag = false
    )
    {
      if(isLeft){
        visLeftVerSearch[currVer] = idxSearch;
        unionFs(DAG_ancestor[currVer], fsBiMat);
        if(
          nMappedParent[currVer] == 0 ||
          (nMappedParent[currVer] == 1 && iecSize[currVer][0] == -1)
        ){
          for(int i=0;i<numLeftVerAdj[currVer];i++){
            Edge* currEdge = &leftVerAdj[currVer][i];
            if(
              isLeafLeft[currVer] == false &&
              currEqSetId[currVer] != 0 &&
              candSpace[currVer].eqSetId[i] != currEqSetId[currVer]
            ){
              continue;
            }
            if(currEdge->residual() > 0
              && visRightVerSearch[currEdge->tail] != idxSearch
            ){
              if(findAugmentedPath(currEdge->tail, !isLeft, length + 1)){
                Edge* inverseEdge = 
                  &rightVerAdj[currEdge->tail][currEdge->inverse_index];
                
                currEdge->flip();
                inverseEdge->flip();
                
                matchedLeftFromRight[currEdge->tail] = currEdge->head;   
                matchedEdgeFromRight[currEdge->tail] = inverseEdge;

                matchedRightFromLeft[currEdge->head] = currEdge->tail;
                matchedEdgeFromLeft[currEdge->head] = currEdge;       
                return true;
              }
            }
          }
        }
        else{
          int repCurrVer;
          if(nMappedParent[currVer] == DAG_parent_query_size[currVer] + 1)
          {
            repCurrVer = currVer;
          }
          else
          {
            repCurrVer = repVerLeft[currVer];
          }
          int iecPos = nMappedParent[currVer] - 1;
          if(iecSize[currVer][iecPos] == -1){
            --iecPos;
          }
          for(int j = 0;j < iecSize[repCurrVer][iecPos]; j++){
            int i = iec[repCurrVer][iecPos][j];
            Edge* currEdge = &leftVerAdj[currVer][i];
            if(
              isLeafLeft[currVer] == false &&
              currEqSetId[currVer] != 0 &&
              candSpace[repCurrVer].eqSetId[i] != currEqSetId[currVer]
            ){ 
              continue;
            }
            if(currEdge->residual() > 0
              && visRightVerSearch[currEdge->tail] != idxSearch
            ){
              if(findAugmentedPath(currEdge->tail, !isLeft, length + 1)){
                Edge* inverseEdge = 
                  &rightVerAdj[currEdge->tail][currEdge->inverse_index];
                
                currEdge->flip();
                inverseEdge->flip();
                
                matchedLeftFromRight[currEdge->tail] = currEdge->head;   
                matchedEdgeFromRight[currEdge->tail] = inverseEdge;

                matchedRightFromLeft[currEdge->head] = currEdge->tail;
                matchedEdgeFromLeft[currEdge->head] = currEdge;
                return true;
              }
            }
          }
        }
      }
      else{//right
        visRightVerSearch[currVer] = idxSearch;
        if(matchedLeftFromRight[currVer] == -1){
          return true;
        }
        for(int i=0;i<numRightVerAdj[currVer];i++){
          Edge* currEdge = &rightVerAdj[currVer][i];
          if(currEdge->residual() > 0
            && visLeftVerSearch[currEdge->tail] != idxSearch
          ){
            if(findAugmentedPath(currEdge->tail, !isLeft, length + 1)){
              Edge* inverseEdge = 
                &leftVerAdj[currEdge->tail][currEdge->inverse_index];
              
              currEdge->flip();
              inverseEdge->flip();
              return true;
            }
          }
        }
      }
      return false;
    }

    bool canPushed(int startLeftVer, int numPusedFlow)
    {
      initFs(fsBiMat);

      ++idxSearch;

      if(idxSearch >= (long long int)(1e18)){
        idxSearch = 1;
        for(int i=0;i<uSequenceSize;i++){
          visLeftVerSearch[i] = 0;
        }
        for(int i=0;i<maxNumDataVertex;i++){
          visRightVerSearch[i] = 0;
        }
      }

      if(!findAugmentedPath(startLeftVer, true)){
        return false;
      }
      else{
        return true;
      }
    }

    void insertEdge(
      int leftVer, int rightVer
    )
    {
      numLeftVerAdj[leftVer]++;
      
      matchedRightFromLeft[leftVer] = -1;
      if(visRightVerAdj[rightVer] != idx_right_adj)
      {
        visRightVerAdj[rightVer] = idx_right_adj;
        numRightVerAdj[rightVer] = 1;
        matchedLeftFromRight[rightVer] = -1;
      }
      else
      {
        numRightVerAdj[rightVer]++;
      }
      leftVerAdj[leftVer][numLeftVerAdj[leftVer] - 1].update(
        leftVer, rightVer, 1, numRightVerAdj[rightVer] - 1
      );
      rightVerAdj[rightVer][numRightVerAdj[rightVer] - 1].update(
        rightVer, leftVer, 0, numLeftVerAdj[leftVer] - 1
      );
    }

    bool insertCandidates(const Graph &query)
    {
      #ifdef DEBUG
      pT(3); printf("k1. insertCandidates is called\n");
      #endif
      fill(currEqSetId, currEqSetId + uSequenceSize, 0);
      for(int i=0;i<uSequenceSize;i++){
        addBitFs(DAG_ancestor[i], i);
        ++numLeftVer;
        const int leftVer = numLeftVer - 1;

        numLeftVerAdj[numLeftVer - 1] = 0;
        int u = i;
        if(query.NECMap[u] != -1){
          isLeafLeft[u] = true;
        }
        else{
          isLeafLeft[u] = false;
        }
        if(query.NECMap[u] != -1 &&
          query.NECMap[u] != u
        ){
          u = query.NECMap[u];
        }
        repVerLeft[leftVer] = u;
        const CandidateSpace &unit = candSpace[u];
        
        for(int j=0;j<unit.size;j++){
          int v = unit.candidates[j];
          int rightVer = v;
          insertEdge(leftVer, rightVer);
        }
        if(!canPushed(leftVer, 1)){
          return false;
        }
      }
      return true;
    }

    bool init(const Graph &query)
    {
      numLeftVer = 0;
      ++idx_right_adj;
      if(!insertCandidates(query)){
        return false;
      }
      assert(numLeftVer == query.nVertex);
      #ifdef DEBUG
      pT(3); printf("g0. init BipartiteMatching\n");
      pT(3); printf("g0. idx : %d\n", idx_right_adj);
      #endif
      return true;
    }

    void releaseLeft(int leftVer)
    {
      int rightVer = matchedRightFromLeft[leftVer];
      Edge *currEdge = matchedEdgeFromLeft[leftVer];
      Edge *inverseEdge = 
        &rightVerAdj[currEdge->tail][currEdge->inverse_index];

      currEdge->flip();
      inverseEdge->flip();

      matchedLeftFromRight[rightVer] = -1;
      matchedRightFromLeft[leftVer] = -1;
    }

    void reconstructLeft(int leftVer)
    {
      #ifdef DEBUG
      pT(4); printf("r2. reconstructLeft is called\n");
      #endif
      Edge *currEdge = matchedEdgeFromLeft[leftVer];
      Edge *inverseEdge = 
        &rightVerAdj[currEdge->tail][currEdge->inverse_index];

      currEdge->flip();
      inverseEdge->flip();

      int rightVer = currEdge->tail;
      matchedLeftFromRight[rightVer] = leftVer;
      matchedRightFromLeft[leftVer] = rightVer;
    }

    bool canReMatchAfterMapping(
      int uCurr, int candMapped, int canMappedIdx
    )
    {
      #ifdef DEBUG
      pT(3); printf("g3. canReMatchAfterMapping is called\n");
      pT(3); printf("g3. uCurr : %d / candMapped : %d\n", uCurr, candMapped);
      pT(3); printf("g3. cond1 : %d\n", (visRightVerAdj[candMapped] != idx_right_adj));
      #endif
      int mappedEqSetId = candSpace[uCurr].eqSetId[canMappedIdx];
      currEqSetId[uCurr] = mappedEqSetId;
      releaseLeft(uCurr);
      
      if(canPushed(uCurr, 1) == true){
        return true;
      }
      else{
        currEqSetId[uCurr] = 0;
        reconstructLeft(uCurr);
        return false;
      }
    }

    bool canReMatchAfterIntersect(int uC)
    {
      int leftVer = uC;

      #ifdef DEBUG
      pT(3); printf("b1. canReMatchAftetIntersect is called\n");
      pT(3); printf("b1. uC : %d\n", uC);
      pT(3); printf("b1. leftVer : %d\n", leftVer);
      #endif

      releaseLeft(leftVer);

      bool isPushed = canPushed(leftVer, 1);
      
      if(isPushed){
        return true;
      }
      else{
        reconstructLeft(leftVer);
        return false;
      }
    }

    void make_components_cell_visited
    (
      const Graph &query
    )
    {
      for(int u=0;u<uSequenceSize;u++)
      {
        CandidateSpace& currSet = candSpace[u];

        if(u == root)
        {
          size_components_cell_visited_per_query_ver[u] = 0;
          for
          (
            int i = currSet.intervals_cell[currEqSetId[u] - 1].first;
            i < currSet.intervals_cell[currEqSetId[u] - 1].second;
            i++
          )
          {
            components_cell_visited_per_query_ver[u][size_components_cell_visited_per_query_ver[u]++] = 
              currSet.components_cell[i];
          }
        }
        else if(query.NECMap[u] == -1)
        {
          int numParentsMapped = nMappedParent[u];
          int iecPos = numParentsMapped - 1;

          Extendable::intersectTwoExtendable
          (
            currSet.components_cell + currSet.intervals_cell[currEqSetId[u] - 1].first,
            currSet.intervals_cell[currEqSetId[u] - 1].second - currSet.intervals_cell[currEqSetId[u] - 1].first,
            iec[u][iecPos],
            iecSize[u][iecPos],
            components_cell_visited_per_query_ver[u],
            size_components_cell_visited_per_query_ver[u]
          );
        }
        else
        {
          size_components_cell_visited_per_query_ver[u] = 0;
          int rep_u = query.NECMap[u];
          int numParentsMapped = nMappedParent[rep_u];
          int iecPos = numParentsMapped - 1;
          for
          (
            int i = 0;
            i < iecSize[rep_u][iecPos];
            i++
          )
          {
            components_cell_visited_per_query_ver[u][size_components_cell_visited_per_query_ver[u]++] = 
              iec[rep_u][iecPos][i];
          }
        }
        bool flag_print = false;
        if(flag_print)
        {
          printf("u : ");
          printf("%d\n", u);
          printf("n m p : %d ", nMappedParent[u]);
          printf("cell(%d):", currEqSetId[u]);
          for
          (
            int i = currSet.intervals_cell[currEqSetId[u] - 1].first;
            i < currSet.intervals_cell[currEqSetId[u] - 1].second;
            i++
          )
          {
            printf("[%d]", currSet.components_cell[i]);
          }
          if(u != root)
          {
            int rep_u = query.NECMap[u];
            printf("visited :");
            int numParentsMapped = nMappedParent[u];
            int iecPos = numParentsMapped - 1;
            printf("sz:%d", iecSize[u][iecPos]);
            for
            (
              int i = 0;
              i < iecSize[u][iecPos];
              i++
            )
            {
              printf("[%d]", iec[u][iecPos][i]);
            }
          }

          puts("");
          printf("\tvisited cell:");
          for
          (
            int i = 0;
            i < size_components_cell_visited_per_query_ver[u];
            i++
          )
          {
            printf
            (
              "[%d, %d]", 
              components_cell_visited_per_query_ver[u][i],
              candSpace[u].candidates[components_cell_visited_per_query_ver[u][i]]
            );
          }
          printf("\tmatched vertex:");
          printf("[%d]", matchedRightFromLeft[u]); 
          printf("\tNECMap : %d\n", query.NECMap[u]);
          puts("");
        }
      }
    }
  
    int cnt_matching = 0;

    bool find_extension_matching
    (
      int u, int num_matching
    )
    {
      bool success_to_find_all_mats = false;
      
      if(u == uSequenceSize)
      {
        ++cnt_matching;
        if(cnt_matching >= num_matching)
        {
          success_to_find_all_mats = true;
        }
      }
      else
      {
        ++nMappedParent[u];
        int iecPos = nMappedParent[u] - 1;
        iecSize[u][iecPos] = 1;
        for
        (
          int idx = 0;
          idx < size_components_cell_visited_per_query_ver[u];
          idx++
        )
        {
          int component_cell = components_cell_visited_per_query_ver[u][idx];
          
          releaseLeft(u);

          iec[u][iecPos][0] = component_cell;

          if(canPushed(u, 1))
          {
            if(find_extension_matching(u + 1, num_matching))
            {
              success_to_find_all_mats = true;
              break;
            }
          }
          else
          {
            reconstructLeft(u);
          }
        }
        iecSize[u][iecPos] = -1;
        --nMappedParent[u];
      }

      return success_to_find_all_mats;
    }

    int get_num_matching
    (
      const Graph &query,
      int num_matching
    )
    {
      
      #ifdef DEBUG
      pT(2); printf("g3. get_num_matching is called\n");
      #endif
      make_components_cell_visited(query);
      cnt_matching = 0;
      find_extension_matching(0, num_matching);
      return cnt_matching;
    }
  };

  void releaseMapping(int depth)
  {
    #ifdef DEBUG
    pT(2); printf("f3. releaseMapping is called\n");
    #endif
    using namespace Queue;
    removeFromQueue(nInserted[depth]); nInserted[depth] = 0;
  }

  void releaseMatching(int depth, const Graph &query)
  {
    #ifdef DEBUG
    pT(2); printf("f4. releaseMatching is called\n");
    #endif
    for(int j = 0; j < DAG_child_query_size[element[depth].vertex]; j++){
      int uC = DAG_child_query[element[depth].vertex][j];
      iecSize[uC][nMappedParent[uC] - 1] = -1;
    }
    BipartiteMatching::currEqSetId[element[depth].vertex] = 0;
  }

  void releaseQueryVertex(int depth)
  {
    #ifdef DEBUG
    pT(2); printf("z.1. releaseQuerVer is called\n");
    #endif
    using namespace Queue;
    int uCurr = element[depth].vertex;
    reinsertToQueue(uCurr, currWeight[depth], position[depth]); 
    element[depth].address = NULL;
    #ifdef DEBUG
    pT(2); printf("z.2. release childern mapped num\n");
    #endif
    for(int i = 0; i < DAG_child_query_size[uCurr]; ++i){
      --nMappedParent[DAG_child_query[uCurr][i]];
    }
  }

  void propagateExtendableQueryVer(int depth)
  {
    const int uCurr = element[depth].vertex;
    for(int i = 0; i < DAG_child_query_size[uCurr]; ++i){
      int uC = DAG_child_query[uCurr][i];
      MappedParent[uC][nMappedParent[uC]++] = uCurr;
      iecSize[uC][nMappedParent[uC] - 1] = -1;
    }
    initFs(element[depth].failingSet);
    element[depth].visEqSetIdCnt++;
  }

  uint64_t fsCurr[MAX_NUM_VERTEX];

  bool tryMapping(
    const Graph& query,
    int depth,
    int extendableCand, int extendableCandIdx
  )
  {
    using namespace Extendable;
    using namespace Queue;
    
    const int uCurr = element[depth].vertex;
    const CandidateSpace &unit = candSpace[uCurr];

    if(
      BipartiteMatching::canReMatchAfterMapping(
        element[depth].vertex, extendableCand, extendableCandIdx
      ) == false
    )
    {
      copyFs(BipartiteMatching::fsBiMat, fsCurr);
      return false;
    }    
    else
    {
      for(int j = 0; j < DAG_child_query_size[uCurr]; ++j)
      {
        int uC = DAG_child_query[uCurr][j];
        if(query.NECMap[uC] == -1 ||
          query.NECMap[uC] == uC
        ){
          long long weight = 
            getExtendableCandidates(
              uC, j,
              uCurr, extendableCandIdx,
              nMappedParent[uC]
          );
          if(
            BipartiteMatching::canReMatchAfterIntersect(uC)
            == false
          )
          {
            copyFs(BipartiteMatching::fsBiMat, fsCurr);
            releaseMapping(depth);
            releaseMatching(depth, query);
            return false;
          }
          if(
            query.NECMap[uC] == -1 && 
            nMappedParent[uC] == DAG_parent_query_size[uC]
          )
          {
            insertToQueue(uC, weight); ++nInserted[depth];
          }
        }
      }
      for(int j = 0; j < DAG_child_query_size[uCurr]; ++j)
      {
        int uC = DAG_child_query[uCurr][j];
        if(query.NECMap[uC] != -1 &&
          query.NECMap[uC] != uC
        ){
          iecSize[uC][0] = iecSize[query.NECMap[uC]][0]; 
          if(
            BipartiteMatching::canReMatchAfterIntersect(uC) == false
          ){
            copyFs(BipartiteMatching::fsBiMat, fsCurr);
            releaseMapping(depth);
            releaseMatching(depth, query);
            return false;
          }
        }
      }
      return true;
    }
  }

  inline void upAddressPos(int depth)
  {
    #ifdef DEBUG 
    pT(2); printf("t1. upAddressPos is called / depth : %d\n",
    depth);
    #endif
    element[depth].addressPos++;
  }

  long long backtrack
  (
    const Graph& query, const Graph& data, int dataGraphID,
    int num_matching = 1
  ) 
  {   
    using namespace Extendable;
    using namespace Queue;

    int cnt_matching = 0;

    long long cntRecur = 0;
    {
      fill(nInserted, nInserted + nQueryVertex, 0);
      fill(nMappedParent, nMappedParent + nQueryVertex, 0);
    
      clearQueue();
      if(BipartiteMatching::init(query)){
        insertToQueue(root, candSpace[root].size);
      }
      else{
        return cntRecur;
      }
    }

    int depth =  0;
    while(depth >= 0)
    {
      if(answer[dataGraphID])
      {//[1]
        #ifdef DEBUG
        pT(1); puts("1.0 already found answer");
        #endif
        int uCurr = element[depth].vertex;
        element[depth].address = NULL; 
        --depth;
      }
      else if(depth == query.nNonLeaf)
      {
        #ifdef DEBUG
        pT(1); printf("2.0. try to get answer. depth : %d\n", depth);
        #endif
        cnt_matching += BipartiteMatching::get_num_matching(query, num_matching - cnt_matching);
        if(cnt_matching >= num_matching)
        {
          answer[dataGraphID] = true;
        }
        #ifdef DEBUG
        pT(1);
        printf("2.1. cnt_matching : %d, num_matching : %d\n", 
        cnt_matching, num_matching);
        #endif 
        --depth;
        if(depth >= 0)
        {
          releaseMapping(depth);
          //element[depth].addressPos++;
          upAddressPos(depth);
          releaseMatching(depth, query);
          makeFull(element[depth].failingSet);
        }
      }
      else if (element[depth].address == NULL)
      {//[3]        
        ++cntRecur;

        int uCurr;
        popFromQueue(uCurr, currWeight[depth], position[depth]);

        #ifdef DEBUG
        pT(1); puts("3.0 extendable new query vertex");
        pT(1); printf("3.1. uCurr : %d / depth : %d\n", uCurr, depth);
        #endif

        if(depth == 0)
        {
          element[depth].update(
            candSpace[uCurr].candidates,
            candSpace[uCurr].size,
            0,
            uCurr
          );
        }
        else
        {
          element[depth].update(
            iec[uCurr][nMappedParent[uCurr]-1],
            iecSize[uCurr][nMappedParent[uCurr]-1],
            0,
            uCurr
          );
        }
        propagateExtendableQueryVer(depth);

        #ifdef DEBUG
        pT(1); printf("3.2. depth : %d\n", depth);
        pT(1); printf("3.3. number of extendable : %d\n",
          element[depth].addressSize
        );
        #endif
      }
      else if(element[depth].addressPos == element[depth].addressSize)
      {//[6]
        #ifdef DEBUG
        pT(1); puts("4.0 already search all extendable candidates");
        pT(1); printf("4.1. depth : %d\n", depth);
        #endif

        releaseQueryVertex(depth);
        const int uCurr = element[depth].vertex;

        if(depth > 0) 
        {
          #ifdef DEBUG
          pT(1); printf("4.2 back to depth from whole searched query vertex: %d\n", depth);
          #endif
          releaseMapping(depth - 1);
          releaseMatching(depth - 1, query);
          if(
            checkRedundant(
              element[depth].failingSet, element[depth - 1].vertex
            ) == true
          ){
            #ifdef DEBUG
            pT(2); printf("4.3.A. redundant: %d\n", depth);
            #endif
            copyFs(element[depth].failingSet, element[depth - 1].failingSet);
            element[depth - 1].addressPos = element[depth - 1].addressSize;
          }
          else{
            #ifdef DEBUG
            pT(2); printf("4.3.A. nonredundant: %d\n", depth);
            #endif
            unionFs(element[depth].failingSet, element[depth - 1].failingSet);
            upAddressPos(depth - 1);
          }
        }
        --depth;
      }
      else
      {//[7]
        #ifdef DEBUG
        pT(1); puts("5.0 choose new extendable candidate");
        pT(1); printf("5.1. depth : %d\n", depth);
        #endif

        const int uCurr = element[depth].vertex;

        int extendableCandIdx;
        int extendableCand;
        if(depth == 0){
          extendableCandIdx = element[depth].addressPos;
        }
        else{
          extendableCandIdx = element[depth].address[element[depth].addressPos];
        }
        extendableCand = candSpace[element[depth].vertex].candidates[extendableCandIdx];
      
        int extEqSetId = candSpace[uCurr].eqSetId[extendableCandIdx];
        
        #ifdef DEBUG
        pT(1);
        printf("5.3. addressPos : %d, extendableCandIdx : %d\n",
        element[depth].addressPos, extendableCandIdx);
        pT(1); printf("5.2. u : %d / v : %d / pos : %d / depth : %d\n",
        element[depth].vertex, extendableCand, 
        element[depth].addressPos, depth);
        pT(1); printf("5.3. eqSetId : %d\n", extEqSetId);
        pT(1); 
        printf("5.4. visEqSetId : %lld visEqSetIdCnt : %lld\n", 
        element[depth].visEqSetId[extEqSetId], 
        element[depth].visEqSetIdCnt);
        #endif

        if(
          element[depth].visEqSetId[extEqSetId]
          == element[depth].visEqSetIdCnt
        ){
          #ifdef DEBUG
          pT(1); printf("5.5.0. already visited cell\n");
          #endif
          upAddressPos(depth);
        }
        else{
          element[depth].visEqSetId[extEqSetId] = element[depth].visEqSetIdCnt;
          if(tryMapping(query, depth, extendableCand, extendableCandIdx)){
            #ifdef DEBUG
            pT(1); printf("5.5.1. mapping success\n");
            #endif
            ++depth;
          }
          else{
            if(
              checkRedundant(fsCurr, element[depth].vertex) == true
            ){
              #ifdef DEBUG
              pT(1); printf("5.5.2. \n");
              #endif
              copyFs(fsCurr, element[depth].failingSet);
              element[depth].addressPos = element[depth].addressSize;
            }
            else{
              #ifdef DEBUG
              pT(1); printf("5.5.3. \n");
              #endif
              unionFs(fsCurr, element[depth].failingSet);
              upAddressPos(depth);
            }
          }
        }
      }
    }
    nMatch += cnt_matching;
    return cntRecur;
  }
}
#endif