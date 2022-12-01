#include <algorithm>
#include <fstream>
#include "randomwalk.h"

using namespace std;

double fRand(double fMin, double fMax);
void PrintQuery(int, RandomWalk* , vector<int>, vector<set<int> >, ofstream&);

//usage ./a.out datagraph_filename #querygraphs #queryvertex sparse_query_filename non-sparse_query_file_name
int main(int argc, char* argv[])
{
	srand(time(NULL));

	RandomWalk walker(argv[1]);
    //cerr << "walker instance ready" << '\n';
	int queryCount = atoi(argv[2]);				//the number of query graphss to generate
	int queryVertexCount = atoi(argv[3]);		//the number of vertices in a query graph
    ofstream sparseQueryFile(argv[4]);          //sparse query file name
    ofstream nonsparseQueryFile(argv[5]);       //non sparse query file name

	for(int i = 0; i < queryCount; i++)
	{	
		int currNumEdge = 0;

		vector<int> vertices;
		vector<set<int> > edges;
		int prevVertex = walker.Initialize();
        //cerr << "walker.Initialize() done" << '\n';
		vertices.push_back(prevVertex);
		set<int> edgeSet;
		edges.push_back(edgeSet);

		while(vertices.size() < queryVertexCount)
		{
			int currVertex = walker.Walk();
			if(currVertex == -1)
				break;

            //if currVertex is a new vertex
			if(std::find(vertices.begin(), vertices.end(), currVertex) == vertices.end()) 
			{
				vertices.push_back(currVertex);
				set<int> s;
				edges.push_back(s);
			}
			
			//...add edge (prev, curr)
			int prevIndex = -1;
			int currIndex = -1;
			for(int j = 0; j < vertices.size(); j++)
			{
				if(vertices[j] == prevVertex)	prevIndex = j;
				if(vertices[j] == currVertex)	currIndex = j;
			}

			//we assume there is no loop
			if(prevIndex != currIndex)
			{
				if(find( edges[prevIndex].begin(), edges[prevIndex].end(), currIndex) == edges[prevIndex].end())
				{
					edges[prevIndex].insert(currIndex);
					edges[currIndex].insert(prevIndex);

					currNumEdge += 2;
				}
			}
			prevVertex = currVertex;

			if(walker.WalkCount() == -1 || walker.WalkCount() > 1000000)
			{
                //cerr << "walker.WalkCount():" <<  walker.WalkCount() << '\n';
				break;
			}
		}
		if(walker.WalkCount() == -1 || walker.WalkCount() > 1000000)
		{
			//cout << "too many walks or walkCount = -1 " << endl << flush;
			i--;
			continue;
		}
		double avgDegreeMin = currNumEdge / (double)queryVertexCount; 
		//PrintQuery(i, &walker, vertices, edges, minQueryFile); //print min

        //cerr << "copying backup_vertices...";
		vector<int> backup_vertices;
		backup_vertices.assign( vertices.begin(), vertices.end() );
        //cerr << "  Done!" << '\n';

        //cerr << "copying backup_edges...";
		vector<set<int> > backup_edges;
		backup_edges.assign( edges.begin(), edges.end() );
        //cerr << "  Done!" << '\n';

		int backupCurrNumEdge = currNumEdge;
		//<end min query>

		//<start calculate edge pool>
		vector<pair<int, int> > edge_pool;
		
		for(int j = 0; j < vertices.size(); j++)
		{
			set<int> orin_edge = walker.AdjacentSet(vertices[j]);
			set<int>::iterator orin_it = orin_edge.begin();
			while(orin_it != orin_edge.end())
			{
				int index = -1;
				for(int k = 0; k < vertices.size(); k++)
				{
					if(vertices[k] == *orin_it)
					{
						index = k;
						break;
					}
				}
				if(index != -1)
				{
					if(std::find(edges[j].begin(), edges[j].end(), index) == edges[j].end())
					{
						pair<int, int> p(j, index);
						edge_pool.push_back(p);
					}
				}
				orin_it++;
			}
		}
		//<end calculate edge pool>
		
		//<start extract induced graph>
        //cerr << "copying backup_edge_pool...";
		vector<pair<int, int> > backup_edge_pool;
		backup_edge_pool.assign( edge_pool.begin(), edge_pool.end() );
        //cerr << "  Done!" << '\n';

		while(edge_pool.size() > 0)
		{
			int x = edge_pool[0].first;
			int y = edge_pool[0].second;
			if(x != y and std::find(edges[x].begin(), edges[x].end(), y) == edges[x].end())
			{
				edges[x].insert(y);
				edges[y].insert(x);
				currNumEdge+=2;
			}
			edge_pool.erase(edge_pool.begin());
		}

		double avgDegreeMax = currNumEdge / (double)queryVertexCount;
		//PrintQuery(g_index, &walker, vertices, edges, greaterQueryFile); //print max

        //cerr << "copying backup_edge_pool...";
		edge_pool.clear();
		edge_pool.assign( backup_edge_pool.begin(), backup_edge_pool.end() );
        //cerr << "  Done!" << '\n';
		//<end extract induced graph>
		//cerr << "avgDegreeMax: "<<avgDegreeMax<<endl;
		
        //<version 2.0>
        //sparse: min~3, non-sparse: (3,max]
        if(avgDegreeMax <= 3 || avgDegreeMin > 3)
        {
            //cerr << "avgDegreeMax: " << avgDegreeMax << ", avgDegreeMin: " << avgDegreeMin << '\n';
            --i;
            continue;
        }
        
        double sparseNumEdge = (int)(fRand(avgDegreeMin, 3) * queryVertexCount);
        double nonsparseNumEdge = (int)(fRand(3, avgDegreeMax) * queryVertexCount);
		
		//<start extarct avg degree graph>
		//pick edge until avgDegree is satisfied
		//stop if there is no more edge

        //cerr << "copying backup_vertices...";
		vertices.clear();
		vertices.assign( backup_vertices.begin(), backup_vertices.end() );
        //cerr << "  Done!" << '\n';

        //cerr << "copying backup_edges...";
		edges.clear();
		edges.assign( backup_edges.begin(), backup_edges.end() );
        //cerr << "  Done!" << '\n';

		currNumEdge = backupCurrNumEdge;
        //sparse
		while(edge_pool.size() > 0 && sparseNumEdge > currNumEdge)
		{
			int rand_index = walker.Random(0, edge_pool.size() - 1);
			int x = edge_pool[rand_index].first;
			int y = edge_pool[rand_index].second;
			if(x != y and std::find(edges[x].begin(), edges[x].end(), y) == edges[x].end())
			{
				edges[x].insert(y);
				edges[y].insert(x);
				currNumEdge+=2;
			}
			edge_pool.erase(edge_pool.begin() + rand_index);
		}
		PrintQuery(i, &walker, vertices, edges, sparseQueryFile); //print sparse query file
        //non sparse
		while(edge_pool.size() > 0 and nonsparseNumEdge >= currNumEdge)
		{
			int rand_index = walker.Random(0, edge_pool.size() - 1);
			int x = edge_pool[rand_index].first;
			int y = edge_pool[rand_index].second;
			if(x != y and std::find(edges[x].begin(), edges[x].end(), y) == edges[x].end())
			{
				edges[x].insert(y);
				edges[y].insert(x);
				currNumEdge+=2;
			}
			edge_pool.erase(edge_pool.begin() + rand_index);
		}
		PrintQuery(i, &walker, vertices, edges, nonsparseQueryFile); //print non-sparse query file
		
        //<end extarct avg degree graph>
	}//for numQuery	

    sparseQueryFile.close();
    nonsparseQueryFile.close();
	return 0;
}

void PrintQuery(int a_graph_index, RandomWalk* a_walker, vector<int> a_vertices, vector<set<int> > a_edges, ofstream& a_file)
{
	//graph id
	a_file << "t # " << a_graph_index << endl;
	//vertex
	for(int j = 0; j < a_vertices.size(); j++)
	{
		a_file << "v " << j << " " << a_walker->Label(a_vertices[j]) << endl;
	}
	//edge
	for(int k = 0; k < a_edges.size(); k++)
	{
		for(set<int>::iterator j = a_edges[k].begin(); j != a_edges[k].end(); j++)
		{
			if(*j > k) //caution: we assume there is no loop (i.e. *j==k)
				a_file << "e " << k << " " << *j << " 0" << endl;
		}
	}
}

//https://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f*(fMax-fMin);
}
