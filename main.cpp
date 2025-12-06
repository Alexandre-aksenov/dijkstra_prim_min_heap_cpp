//Assignment for week 4 of the course "C++ for C programmers". Prim's algorithm. The classes of the assignment of Week 3 have been used and modified. 

//Three nontrivial classes are defined and their methods are provided.
//Graph: a symmetric graph with labeled edges (labels: unsigned int). Supports initialization from file or by copy from another instance of Graph.
//Priority Queue: the min-heap structure. Supports addition of a new element with given priority, extraction of minimal priority, decreasing priority of an element. The set of priorities is a template type.
//VertexLabeledGraph: a Graph together with vertex labels and a Priority Queue, containing a subset of vertices. The mehods of this class keep the labels and the Queue consistent. Supports initialisation from a Graph, Dijkstra's and Prim's algorithms.

//An attempt has been made in creating a "universal" Priority Queue class, supporting both Integer and Double types (adapted to both assignments of Week 3 and Week 4), and in creating two different types for Graphs (Dijkstra's and Prim's algorithms do not have access to modifying the edges of the graph). 

// The extension towards MST with different colors has not not being adressed for the moment. The author begs for reviewer's understanding regarding this point.


#include <iostream>
#include <vector>
#include <list>
#include <tuple>

//For reading from file:
#include<fstream>
#include<iterator>


using namespace std;


const unsigned int FlagNonAdj=0;

/*
Class 'Graph': a graph represented by two adjacency matrices: 

    1. AdjValues: symmetric matrix of booleans. 

    2. Costs: symmetric matrix of unsigned int, 0 on diagonal.
        The values are interpreted as:
        value>=0 : distance btw 2 vertices
        0 if vertices are not adjacent.

Edge weights may be zero or non-negative. 
*/


class Graph{
public:
    // Constructor: totally disconnected graph
    Graph(unsigned int N=20){ 
        size=N;
        
        AdjValues = new bool* [N];
        for (unsigned int i=0; i<size; ++i)
           AdjValues[i]= new bool[N];
        
        for (unsigned int i=0; i<size; ++i){
            for (unsigned int j=0; j<size; ++j)
                AdjValues[i][j]=false;
        }
        
        Costs = new unsigned int* [N];
        for (unsigned int i=0; i<size; ++i)
           Costs[i]= new unsigned int[N];
        
        for (unsigned int i=0; i<size; ++i){
            for (unsigned int j=0; j<size; ++j)
                Costs[i][j]=FlagNonAdj;
        }
    }
    
    // Constructor from file.
    Graph(istream_iterator<unsigned int> start, istream_iterator<unsigned int> end){ 
        vector<unsigned int> raw_nbs(start, end);
    
        //Check the residue modulo 3.
        if ((raw_nbs.size() % 3) !=1)
        {
            cout << "Wrong file format.\n" ;
            throw ("Wrong file format.");
        }

        //Reading the first value.
        auto it= raw_nbs.begin();
        //Read the elements of 'raw_nbs' as formatted arrays.
    
        size=  *it; //raw_nbs[1]-> unsigned int : V(G)
        ++it;
            
        //interpret other entries as: extremity1 (unsigned int) _ extremity2 (unsigned int) _ Cost (unsigned int)

        //The following lines replicate the base constructor Graph(size)
        AdjValues = new bool* [size];
        for (unsigned int i=0; i<size; ++i)
           AdjValues[i]= new bool[size];
        
        for (unsigned int i=0; i<size; ++i){
            for (unsigned int j=0; j<size; ++j)
                AdjValues[i][j]=false;
        }
        
        Costs = new unsigned int* [size];
        for (unsigned int i=0; i<size; ++i)
           Costs[i]= new unsigned int[size];
        
        for (unsigned int i=0; i<size; ++i){
            for (unsigned int j=0; j<size; ++j)
                Costs[i][j]=FlagNonAdj;
        }
        
        //Add edges
        while (it != raw_nbs.end())
        {
            unsigned int extr1 = *it; 
            ++it;
        
            unsigned int extr2 = *it; 
            ++it;
        
            unsigned int cost = *it; 
            ++it;

            addEdge( extr1, extr2, cost); // This "adds" every edge twice
        }        
    }

    // Constructor by deep copy. Used in constructing a VertexLabeledGraph.
    Graph(const Graph& src){
        size= src.size;
        
        AdjValues = new bool* [size];
        for (unsigned int i=0; i<size; ++i)
           AdjValues[i]= new bool[size];
        
        for (unsigned int i=0; i<size; ++i){
            for (unsigned int j=0; j<size; ++j)
                AdjValues[i][j]=src.adjacent(i,j);
        }
        
        Costs = new unsigned int* [size];
        for (unsigned int i=0; i<size; ++i)
           Costs[i]= new unsigned int[size];
        
        for (unsigned int i=0; i<size; ++i){
            for (unsigned int j=0; j<size; ++j)
                Costs[i][j]=src.get_edge_value(i,j);
        }
    }
    
    ~Graph(){
        for (unsigned int i=0; i<size; ++i)
            delete AdjValues[i];
        
        delete AdjValues; 
    
        for (unsigned int i=0; i<size; ++i)
            delete Costs[i];
    
        delete Costs; 
    }
    
    unsigned int V() const {
        return size;
    }
    
    unsigned int E() const {
        unsigned int count=0;
        for (unsigned int i=0; i<size; ++i){
            for (unsigned int j=i+1; j<size; ++j)
                if (AdjValues[i][j])
                    count++;
        }
    return count;
    }
    
    
    bool adjacent(unsigned int x, unsigned int y) const {
        if ((x<size)&&(y<size))
            {return (AdjValues[x][y]);} 
        else
            throw invalid_argument( "index out of bounds" );
    }
    
    //lists all nodes y such that there is an edge from x to y.
    vector<unsigned int> neighbors(unsigned int x) const { 
        vector<unsigned int> result;
        for (unsigned int i=0; i<size; ++i){
            if (adjacent(i,x))
                result.push_back(i);
        }
        return result; 
    }
    
    void addEdge(unsigned int x, unsigned int y, unsigned int value=1){
        
        if ((x<size)&&(y<size))
            {
            AdjValues[x][y]=true;
            AdjValues[y][x]=true;
            Costs[x][y]= value;
            Costs[y][x]= value;
            }
        else
            throw invalid_argument( "index out of bounds" );
    }
    
    void deleteEdge(unsigned int x, unsigned int y){
        if ((x<size)&&(y<size))
            {
            AdjValues[x][y]=false;
            AdjValues[y][x]=false;
            Costs[x][y]=FlagNonAdj;
            Costs[y][x]=FlagNonAdj;
            }
        else
            throw invalid_argument( "index out of bounds" );
    }
    
    unsigned int get_edge_value(unsigned int x, unsigned int y) const {
        if ((x>=size)||(y>=size))
            throw invalid_argument( "index out of bounds" );
        else
           return Costs[x][y];
    }
    
    void set_edge_value (unsigned int x, unsigned int y, unsigned int v){
        if ((x>=size)||(y>=size))
            throw invalid_argument( "index out of bounds" );
        else
            if (adjacent(x,y))
                {Costs[x][y]=v;
                 Costs[y][x]=v;}
    }


    bool is_symmetric() const { //checking correction of a Graph, only useful for debugging
        
        //check no vertex is a loop
        for (unsigned int i=0; i<size; ++i){
            if (AdjValues[i][i]) return false;
            if (Costs[i][i]>0) return false;
        }
        
        //Check AdjValue is symmetric
        for (unsigned int i=0; i<size; ++i){
            for (unsigned int j=i+1; j<size; ++j){
                if (AdjValues[i][j] != AdjValues[j][i]) return false;
            }
        };

        
        for (unsigned int i=0; i<size; ++i){
            for (unsigned int j=i+1; j<size; ++j){
                if (Costs[i][j] != Costs[j][i]) return false;
            }
        };
        
        //Check AdjValue corresponds to Costs
        for (unsigned int i=0; i<size; ++i){
            for (unsigned int j=0; j<size; ++j){
                if ( (Costs[i][j]>0) != AdjValues[i][j] ) return false;
            }
        };
        
        return true;
    }

private:
    unsigned int size; //enforce const correctness?
    bool ** AdjValues;
    unsigned int ** Costs;
    
};


//Class Priority Queue.
//labels: template (unsigned int).
template <class T>
void PrintVector(vector<T> vec) //For checking correctness
{
    for(unsigned int i=0; i<vec.size(); ++i)
    cout << vec[i] << ' ';
    cout << endl;
}

template <class LabelType> //Must allow comparison, and have a fixed value '0'. Can be specialized to 'double' or 'unsigned int'.
class heapElem
{public:
    heapElem(unsigned int ind=0, LabelType lbl=0):VertexIndex(ind),label(lbl){};
    unsigned int VertexIndex;
    LabelType label;
};


template <class LabelType> //Must allow comparison, and have a fixed value '0'. Can be specialized to 'double' or 'unsigned int'.
class PriorityQueue{
public:
    PriorityQueue(unsigned int MaxNbIndexes=50): Arr(vector<heapElem<LabelType>>(1)), MaxLbl(0), VertexLocations(vector<unsigned int>(MaxNbIndexes)) {}
        //vector of 1 element, representing empty queue

    ~PriorityQueue(){
        Arr.erase(Arr.begin(),Arr.end()); //can this be improved?
        VertexLocations.erase(VertexLocations.begin(),VertexLocations.end());
    }
    
    void ClearPQ(){
    //Erases the contents of a PQ, leaving it empty.
    //see: https://stackoverflow.com/questions/8848575/fastest-way-to-reset-every-value-of-stdvectorint-to-0
        Arr.erase((Arr.begin()+1),Arr.end());
        MaxLbl=0;
        std::fill(VertexLocations.begin(),VertexLocations.end(),0);
    }
    
    
    
    vector<LabelType> ArrayLabels(){ //return the tree of labels in form of an array. Returns the 'meaningless' element nb 0 as the first element of output
        unsigned int N=Arr.size(); //The size of heap (nb of meaningful elements) is N-1.
        vector<LabelType> res(N);
        for (unsigned int i=0; i<N; ++i)
            res[i]=Arr[i].label;
        return res;
    }
    
    vector<unsigned int> ArrayIndices(){ //return the tree of indices in form of an array. Returns the 'meaningless' element nb 0 as the first element of output
        unsigned int N=Arr.size();
        vector<unsigned int> res(N);
        for (unsigned int i=0; i<N; ++i)
            res[i]=Arr[i].VertexIndex;
        return res;
    }
    
    
    vector<unsigned int> CopyVertexLocs(){ //return the "locations" of the elements, which are already in the heap 
        unsigned int L=VertexLocations.size();
        vector<unsigned int> res(L);
        for (unsigned int i=0; i<L; ++i)
            res[i]=VertexLocations[i];
        return res;
    }
    
    void heapify(unsigned int i){ //i : index in the array Arr, which is (possibly) violating the min-heap property
    // If it does, the element number 'i' descends in the heap.
    // Called after extracting the minimal element
        unsigned int N=Arr.size();
        if ((i==0)||(i>=N))
            {cout << "index out of bounds" << endl;}
        else
            {
            unsigned int left  = 2*i;
            unsigned int right = 2*i+1;
            unsigned int smallest;
        
            if((left <= N) && (Arr[left].label < Arr[i].label) )
                smallest = left;
            else
                smallest = i;
        
            if((right <= N) && (Arr[right].label < Arr[smallest].label) )
                smallest = right;

            if(smallest != i)
                {
                SwapPQ (i,smallest);
                heapify(smallest);
                }
            }
    }
    
    void chgPriority(unsigned int ind, LabelType priority){
    // called by "insert" and "chgPriorityVertex"
    // 'ind' is the index in the array. 'priority' is the new priority for the element 'ind', smaller than the ancient.
        unsigned int N=Arr.size();
        
        if ((ind==0)||(ind>=N))
            {cout << "index out of bounds" << endl;}
        else
            {
            if(priority > Arr[ ind ].label)
                {
                cout<< "New value is less than current value, can’t be inserted" <<endl;
                return;
                }
            Arr[ ind ].label = priority;
            unsigned int i=ind; //'i' points to the element, whose priority has been changed, while this element goes up in the heap.
            while( (i > 1) && (Arr[ i/2 ].label > Arr[ i ].label))
                {
                SwapPQ(i/2 , i );
                i = i/2;
                }
            }
    }
    
    void chgPriorityVertex(unsigned int vert, LabelType priority){
        //called in Dijkstra's algo, line 20.
        if ((vert>= VertexLocations.size())||(VertexLocations.at(vert) ==0))
        {
            cout << "Vertex does not exist in the heap" << endl;
        }
        else
        {
            chgPriority(VertexLocations.at(vert),priority);
        }
    
    }

    heapElem<LabelType> minPriority(){
        // REMOVES the top element of the heap, and returns heapElem. Called in Dijkstra's algo, line 15.
        unsigned int N=Arr.size();
    
        if (N<=1)
        {
            cout << "Trying to remove an element from empty heap." << endl;
            throw  invalid_argument( "empty heap" );
        }
        else
        {
            //Save the top vertex for return
            heapElem<LabelType> res=Arr[1];
    
            //"Remove" the top vertex from heap
            VertexLocations[res.VertexIndex]=0;
        
            if (N>2)
            {
                //Move the "bottom" vertex to top
                Arr[1]=Arr.at(N-1);
                VertexLocations[Arr[1].VertexIndex]=1;
    
                Arr.pop_back();
    
                //Heapify the vertex, which is now on top
                heapify(1);
            }
            else
            {
                Arr.pop_back(); //No vertex to move
            }
        
            return res;
        }
    
    }
    
    bool contains(unsigned int queue_element){
        // Check whether 'queue_element' is a vertex index. Called in Dijkstra's algo, line 17.
        if (queue_element>=VertexLocations.size())
        {
            cout << "Index of the vertex queried is too large" <<endl;
            return false;
        }
        else
        {
            return (VertexLocations.at(queue_element)>0);
        }
        
    }
    
    void Insert(unsigned int vindex, LabelType lbl){ //insert a new entry  with label 'lbl' and vertex index 'vindex'
        if ((vindex>= VertexLocations.size() ) )
        {
            cout << "the vertex is out of bounds" << endl;
            throw  invalid_argument( "index out of bounds" );
        }
        else
            if (VertexLocations.at(vindex)>0) //can be replaced by call of "contains"
            {
                cout << "the vertex already exists" << endl;
            }
            else
            {
        
    //insert with a bigger label, than all existing elements.
                unsigned int N=Arr.size();
                LabelType biggerLbl;
        
                if (lbl< MaxLbl)
                    biggerLbl = MaxLbl;
                else
                {
                    biggerLbl = lbl;
                    MaxLbl = lbl;
                };
        
                Arr.push_back(heapElem<LabelType>(vindex,biggerLbl));
                VertexLocations[vindex]=N;
        
    //decrease the label to 'lbl', which is smaller
                chgPriority(N, lbl);
            }
    
    }
    
    
    unsigned int top(){
        //Returns: the vertex index of the top element, without removing. Called in Dijkstra's algo, line 13.
        unsigned int N=Arr.size();
    
        if (N<=1)
        {
            cout << "Trying to remove an element from empty heap." << endl;
            throw  invalid_argument( "empty heap" );
        }
        else
        {
            return(Arr.at(1).VertexIndex);
        }
    }
    
    unsigned int size(){
        // Called in Dijkstra's algo, line 12 in the restricted form "Q is not empty".
        unsigned int N=Arr.size();
    
        if (N==0)
        {
            cout << "The heap is invalid, showing Arr==0." << endl; //This should never happen
            throw  invalid_argument( "Invalid heap" );
        }
        else
        {
            return (N-1);
        }
    }
    

    //Checks the validity of heap: Arr[0] should exist;
    //each label should be bigger than its parent's;
    //each VertexIndex should correspond to VertexLocations.
    //This checking method should always return true. Its complexity is O(N) 
    bool is_valid() 
    {
        unsigned int N=Arr.size();
        if (N==0)
            return false;
        
        
        //if VertexLocations too small
        unsigned int L=VertexLocations.size();
        if (L<N-1)
            return false;
        
        for (unsigned int i=2; i<N; ++i)
        {
            if (Arr.at(i).label< Arr.at(i/2).label)
                return false;
        }
        
        for (unsigned int i=1; i<N; ++i)
        {
            if (L<=(Arr.at(i)).VertexIndex)
                return false;
            if (VertexLocations.at((Arr.at(i)).VertexIndex) != i)
                return false;
        }
        
        for (unsigned int j=0; j<L; ++j)
        {
            if (VertexLocations[j]>=N)
                return false;
            
            if ((VertexLocations[j]!=0)&& ( Arr.at(VertexLocations[j]).VertexIndex != j) )
                return false;
        }
        
        return true;

    }

private:
    void SwapPQ(unsigned int i, unsigned int j){ //Swaps the entries i and j in the array, and the corresponding entries in the array of Vertex Locations.
        if ((i==0)||(j==0)||(i>=Arr.size())||(j>=Arr.size())) //this should never happen
            {cout << "index out of bounds" << endl;}
        else
            {
            std::swap(VertexLocations[Arr[i].VertexIndex],VertexLocations[Arr[j].VertexIndex]);
            
            std::swap(Arr[i],Arr[j]);
            /*
            unsigned int tempVertexIndex=Arr[j].VertexIndex;
            double tempLabel=Arr[j].label;
            
            Arr[j].VertexIndex=Arr[i].VertexIndex;
            Arr[j].label=Arr[i].label;
            
            Arr[i].VertexIndex=tempVertexIndex;
            Arr[i].label=tempLabel;*/
            }
    }
    
    vector<heapElem<LabelType>> Arr;
    //MaxLbl>=0 is the hignest label, which has been observed in the history of the heap. It is set at zero at the heap's creation, updated by method "Insert".
    LabelType MaxLbl; 

    //Locations of all elements of the initial graph (bijection to VertexIndex). This array is modified at every swapping. 
    //Used in Dijkstra's algorithm at the step of finding dist[v] (line 18)
    vector<unsigned int> VertexLocations;
};


void PrintPair(pair<unsigned int, unsigned int> p){
    cout << "(" << p.first << " , " << p.second << ")";
}



void PrintPairVector(vector<pair<unsigned int, unsigned int>>  PVect){
    // Used for testing the Prim implementation
    for (auto p : PVect)
    {
        PrintPair(p);
        cout << endl;
    }
}


//Class VertexLabeledGraph contains fields:
//Graph G
//VertexLabel * Labels
//PriorityQueue Q

//Methods:
//constructor of empty graph
//constructor from Graph (simple)
//destructor
//set label (5 methods)
//
//PrimMST

//Properties of valid-ness: 
//The underlying graph and the priority queue should be valid;
//The priority queue should contain exactly the vertices with status 1 (Labels[x].status == 1)  and 
//should satisfy Q.Arr[VertexLocations[x]].label == Labels[x].label wherever the vertex has status 1.

//Printing a list,
//see https://stackoverflow.com/questions/16229729/printing-out-contents-of-a-list-from-the-c-list-library
void PrintList(list<unsigned int> lst)
{for (list<unsigned int>::iterator i = lst.begin(); i != lst.end(); ++i)
    cout << *i << " ";
cout << endl;
}


class VertexLbl
{public:
    VertexLbl(unsigned int prevind=0, unsigned int lbl=0):prev(prevind),label(lbl),status(0){};
    
    unsigned int prev;
    unsigned int label;
    int status;
};

class VertexLabeledGraph{
public:
    VertexLabeledGraph(unsigned int N=20): size(N), G(Graph(N)), Labels(vector<VertexLbl> (N)), Q(PriorityQueue<unsigned int> (N)){};
    
    //A vertex-labeled copy of a graph.
    VertexLabeledGraph(const Graph Ginput): size(Ginput.V()), G(Ginput), Labels(vector<VertexLbl> (Ginput.V())), Q(PriorityQueue<unsigned int> (Ginput.V())){}; 
    
    
    ~VertexLabeledGraph(){//the destructors of G and of Q will be called automatically
        Labels.erase(Labels.begin() , Labels.end());
    }
    
    //Output V(G), E(G) to check correctness
    unsigned int V(){
        return G.V();
    }
    
    
    unsigned int E(){
        return G.E();
    }
    
    void printQindexes(){ //for testing
        PrintVector(Q.ArrayIndices()); 
        cout << "\n";
    }
    
    //2 functions for reading labels
    
    int get_node_status(unsigned int x){
        if (x<size)
            {return Labels[x].status;}
        else
            throw invalid_argument( "index out of bounds" );
    }
    
    VertexLbl get_node_value (unsigned int x){
        if (x<size)
            {return Labels[x];}
        else
            throw invalid_argument( "index out of bounds" );
    }
    
    //3 functoins, which modify labels and the Queue accordingly.
    heapElem<unsigned int> set_status_pop_open(){ //Removes the element of lowest priority from heap and Returns it. status: 1->2. Used in Dijksta, l.13,15
        if (Q.size() == 0)
        {
            cout << "Trying to extract a vertex from empty Open Set\n";
            throw invalid_argument("Empty Open Set");
        }
        else
        {
            heapElem<unsigned int> vheap=Q.minPriority();
            if (get_node_status(vheap.VertexIndex) != 1)
            {
                cout << "The graph claims that the min-label vertex is not in the open set.\n";
                throw invalid_argument("Inconsistent Open Set btw Graph and Queue");
            }
            else
            {
                set_node_status(vheap.VertexIndex,2);
                return vheap;
            }
        }
    }
    
    void set_status_insert(unsigned int x, unsigned int prevind, unsigned int  a){ //Changes the labels of 'x' to 'a', adds to Q. status: 0->1. 
    //Called in Dijkstra' algo , lines 20,21 
        if (x>=size)
            throw invalid_argument( "index out of bounds" );
        else
            {
                if ((Labels[x].status!=0)|| (Q.contains(x)))
                    throw invalid_argument( "Vertex has been visited already." );
                else
                {
                    set_node_previous(x,prevind);
                    set_node_label(x,a);
                    set_node_status(x,1);
                    Q.Insert(x,a);
                }
            }
    }

    void set_status_decreaselbl(unsigned int x, unsigned int prevind, unsigned int a){ //Changes the labels of 'x' to 'a', adds to Q. status remains 1
    //Called in Dijkstra' algo , lines 20,21  
        if (x>=size)
            throw invalid_argument( "index out of bounds" );
        else
            {
                if ((Labels[x].status!=1)|| (!Q.contains(x)))
                    throw invalid_argument( "Vertex is not in open set." );
                else
                {
                    set_node_previous(x,prevind);
                    set_node_label(x,a);
                    Q.chgPriorityVertex(x,a);
                }
            }
    }
    
    //Dijkstra and Prim algorithms.
    //Utility functions    
    void ClearLabels(){ 
        for (unsigned int i=0; i<size; ++i)
        {
            Labels[i].prev=0;
            Labels[i].label=0.0;
            Labels[i].status=0;
        }
        Q.ClearPQ();
    }

    std::list<unsigned int> TracePathToVertex(unsigned int init, unsigned int w){ //trace the path using the information in the label "previous".
        list<unsigned int> res;
        unsigned int current=w;
        double currentlbl=get_node_value(w).label;
        res.push_front(w);
        
        while (current!=init)
        {
            unsigned int previous=get_node_value(current).prev;
            if (get_node_value(previous).label < currentlbl)
            {
                res.push_front(previous);
                current = previous;
                currentlbl=get_node_value(previous).label;
            }
            else
            {
                cout << "The labels are not decreasing" << endl;
                break;
            }
        }
        
        return res;
    }
    
    //Outputs a Minimal Spanning Tree as a vector of edges. Take the vertices with status 2 (the Connected component of the initial vertex) into account. Appends the pair (init, sum of costs) at the end.
    vector<pair<unsigned int, unsigned int>> MakeTree(unsigned int init=0){
        if (init>=size)
            throw invalid_argument( "index out of bounds" );
        else
        {
            vector<pair<unsigned int, unsigned int>> res;
            unsigned int TotalCost=0;
            for (unsigned int i=0;  i<size; ++i)
            {
                if ((get_node_status(i) == 2)&&(i != init))
                {
                    res.push_back( make_pair(i,get_node_value(i).prev));
                    TotalCost += get_node_value(i).label;
                }
            }
            
            res.push_back( make_pair(init, TotalCost) );
            return res;
        }
        
    }
    
    // Dijkstra algorithm
    // see https://en.wikipedia.org/wiki/Dijkstra's_algorithm
    std::list<unsigned int> path(unsigned int u, unsigned int w, int PrintInterm=0){
        vector<unsigned int> NeighborsV(V());
        
        ClearLabels(); //lines 3,5,6,7 
        
        set_status_insert(u,u,0); //line 10
        if(PrintInterm) { cout << "Node " << u << " to open set with label 0\n";}
        
        while (Q.size()>0) //Main loop l.12
        {

            heapElem<unsigned int> vheap= set_status_pop_open(); //lines 13,15
            unsigned int v=vheap.VertexIndex  ;
            if(PrintInterm) { cout << "Node " << v << " to closed set with label " << vheap.label <<endl;}
            if (v!=w)
            {
            
                NeighborsV=G.neighbors(vheap.VertexIndex);
                
                for(unsigned int i=0; i<NeighborsV.size(); ++i)
                {if (get_node_status(NeighborsV[i]) <2) //Neighbor not in the closed set: l.17 
                {
                    unsigned int alt= vheap.label +  G.get_edge_value(v,NeighborsV[i]); // line 18
                    if (get_node_status(NeighborsV[i]) ==0) //Add to Open set, line 19 
                    {
                        set_status_insert(NeighborsV[i],v,alt); //lines 20,21
                        if(PrintInterm) { cout << "Node " << (NeighborsV[i]) << " to open set with label " << alt <<endl;}
                    }
                    else 
                    {
                        if (alt<get_node_value(NeighborsV[i]).label) // line 18
                        {
                            set_status_decreaselbl(NeighborsV[i],v,alt); //lines 20,21
                            if(PrintInterm) { cout << "Node " << (NeighborsV[i]) << " changes its label to " << alt <<endl;}
                        }
                    }
                }
            
                }
            
            }
            else break; //if v==w , the algo found the path required, break the while loop 
        }

        std::list<unsigned int> res;        
        if (get_node_value(w).status ==2)
        {
            res=TracePathToVertex(u,w);
        }
        
        return res;
    }
    
    //PrimMST
    //see: https://en.wikipedia.org/wiki/Prim%27s_algorithm
    //The steps are mentioned in comments to the present implementation 
    
    // Simplifying assumption : the graph is connected"
    // Output: the cost and the tree (the edges)
    //see: https://www.quora.com/How-can-I-encode-a-tree-into-a-string-format-such-that-the-tree-can-be-constructed-back-from-the-string-encoding-Also-encoding-should-be-efficient-both-time-and-space-wise

    //If the graph is disconnected, output the Spanning Tree of the connected component of the initial vertex 'u'.
    //Append the initial vertex and the cost to the output vector.
    vector<pair<unsigned int, unsigned int>> PrimMST(unsigned int u, int PrintInterm=0){
        vector<unsigned int> NeighborsV(V());
        
        ClearLabels(); //instructions 1,2
        
        set_status_insert(u,u,0); //instructions 1,2
        if(PrintInterm) { cout << "Node " << u << " to open set with label 0\n";}
        
        while (Q.size()>0) //Main loop, instruction 3
        {

            heapElem<unsigned int> vheap= set_status_pop_open(); //instructions 3a, 3b
            unsigned int v=vheap.VertexIndex  ; 
            if(PrintInterm) { cout << "Node " << v << " to closed set with label " << vheap.label <<endl;}
            
            NeighborsV=G.neighbors(vheap.VertexIndex); 
                
            for(unsigned int i=0; i<NeighborsV.size(); ++i) //instruction 3c
            {if (get_node_status(NeighborsV[i]) <2) //Neighbor not in the closed set: l.17 
            {
                unsigned int alt= G.get_edge_value(v,NeighborsV[i]); // line 18
                if (get_node_status(NeighborsV[i]) ==0) //Add to Open set, line 19 
                {
                    set_status_insert(NeighborsV[i],v,alt); //lines 20,21
                    if(PrintInterm) { cout << "Node " << (NeighborsV[i]) << " to open set with label " << alt <<endl;}
                }
                else 
                {
                    if (alt<get_node_value(NeighborsV[i]).label) // line 18
                    {
                        set_status_decreaselbl(NeighborsV[i],v,alt); //lines 20,21
                        if(PrintInterm) { cout << "Node " << (NeighborsV[i]) << " changes its label to " << alt <<endl;}
                    }
                }
            }
            }
            
        }
        
        return MakeTree(u); //instruction 4
    }
    
    bool is_valid()
    {
        if ((G.V() != size)|| (Q.CopyVertexLocs().size() != size) || (Labels.size() != size))
        {
            cout << "The sizes of the labeled graph do not match";
            return false;
        }
        
        //The underlying graph and the priority queue should be valid;
        if (!(G.is_symmetric() ))
        {
            cout << "The underlying graph in not valid";
            return false;
        }
        
        if (!(Q.is_valid() ))
        {
            cout << "The underlying queue in not valid";
            return false;
        }

        
        //The priority queue should contain exactly the vertices with status 1 (VertexLbl[x].status == 1)  
        for (unsigned int i=0; i<size; ++i)
        {
            if ( (Labels[i].status == 1) != Q.contains(i)  )
            {
                cout << "The graph labels and the queue do not agree about the status of the vertex" << i <<endl;
                return false;
            }
        }
        
        //The labels shoul be consistent between the Queue and the array Labels.
        vector<unsigned int> VertQ = Q.ArrayIndices();
        vector<unsigned int> LabelsQ = Q.ArrayLabels();
        for (unsigned int x=0; x<VertQ.size(); ++x) 
        {    
            if (LabelsQ.at(x) != Labels.at(VertQ.at(x)).label)
            {
                cout << "The graph labels and the queue do not agree about the current label of a vertex" << x <<endl;
                return false;
            }
        }
        
        return true;
    }
    
private:
    const unsigned int size;
    const Graph G;
    vector<VertexLbl> Labels;
    PriorityQueue<unsigned int> Q;
    
    
    //3 methods, which modify one label
    
        void set_node_previous(unsigned int x, unsigned int prevind){
        if (x<size)
            {Labels[x].prev=prevind;}
        else
            throw invalid_argument( "index out of bounds" );
    }
    
    void set_node_label(unsigned int x, double a){
        if (x<size)
            {Labels[x].label=a;}
        else
            throw invalid_argument( "index out of bounds" );
    }
    
    
    
    void set_node_status(unsigned int x, int stat){
        if (x<size)
            {Labels[x].status=stat;}
        else
            throw invalid_argument( "index out of bounds" );
    }
};


// Five tests.
//Test: create a very small graph (0-(1)-1-(2)-2), then destroy
void test_create_graph()
{
    cout << "-- Test a simple graph. --" << endl;

    //Create graph 
    Graph Gtest(3);

    Gtest.addEdge(0,1);
    Gtest.addEdge(1,2,2.0);

    cout << "Is I2 symmetric? " << Gtest.is_symmetric() << endl;

    //Print neighbors of vertex 2
    vector<unsigned int> neighbors1=Gtest.neighbors(1);
    cout<< "The neighbors of 1 are: ";

    for(unsigned int i=0; i<neighbors1.size(); ++i)
        cout << neighbors1[i] << ' ';

    cout << endl;

    //Print "lengths" of all 3 edges
    cout << "The 'Lengths' of edges of I2 are: 0-(" << Gtest.get_edge_value(0,1) << ")-1-(" << Gtest.get_edge_value(1,2) << ")-2-(" << Gtest.get_edge_value(2,0) << ")\n";

    //Expand to K3
    Gtest.addEdge(0,2,3.0);
    cout << "Is K3 symmetric? " << Gtest.is_symmetric() << endl;
    cout << endl;
}


// Test: Create a very small graph ( 0-(1)-1-(2)-2 ), load it into a VertexLabeledGraph, check correctness.
void test_VertexLabeledGraph()
{ 
    cout << "-- Test a vertex-labeled graph. --" << endl;

    //Create graph
    Graph Gtest(3);

    Gtest.addEdge(0,1);
    Gtest.addEdge(1,2,2.0);

    cout << "Is I2 symmetric? " << Gtest.is_symmetric() << endl;


    //Print neighbors of vertex 2
    vector<unsigned int> neighbors1=Gtest.neighbors(1);
    cout<< "The neighbors of 1 are: ";

    for(unsigned int i=0; i<neighbors1.size(); ++i)
        cout << neighbors1[i] << ' ';


    cout << endl;

    //Print "lengths" of all 3 edges
    cout << "The 'Lengths' of edges of I2 are: 0-(" << Gtest.get_edge_value(0,1) << ")-1-(" << Gtest.get_edge_value(1,2) << ")-2-(" << Gtest.get_edge_value(2,0) << ")\n";

    //Expand to K3
    Gtest.addEdge(0,2,3.0);
    cout << "Is K3 symmetric? " << Gtest.is_symmetric() << endl;

    //Expand to VertexLabeledGraph
    VertexLabeledGraph GTestLabeled(Gtest);
    cout << "Is labeled K3 valid? " << GTestLabeled.is_valid() << endl;
    
    cout << "It has " << GTestLabeled.V() << " vertices," << endl;
    cout << "and " <<    GTestLabeled.E() << " edges," << endl;

    Gtest.deleteEdge(1,2); //modify the original
    cout << "The original has been modified.\n";

    cout << "Now, the original graph graph has " << Gtest.V() << " vertices,\n and " << Gtest.E() << " edges.\n";
    cout << "Is it symmetric? " << Gtest.is_symmetric() << endl;
    

    cout << "The labeled test graph has " << GTestLabeled.V() << " vertices,\n and " << GTestLabeled.E() << " edges.\n";
    cout << "Is labeled K3 still valid? " << GTestLabeled.is_valid() << endl;
    
    //Testing the methods, which modify the queue.
    GTestLabeled.set_status_insert(0,0,0);
    cout << "Vertex 0 has been added to the open set.\n";
    cout << "The queue is now: ";
    GTestLabeled.printQindexes();
    
    heapElem<unsigned int> zeroLabeled = GTestLabeled.set_status_pop_open();
    cout << "Vertex 0 has been removed from the open set.\n";
    cout << "The queue is now: "; 
    GTestLabeled.printQindexes();
    cout << "The weight of vertex 0 is " << zeroLabeled.label << endl;
    
    GTestLabeled.set_status_insert(1,0,1);
    GTestLabeled.set_status_insert(2,0,2);
    cout << "Vertices 1,2 have been added to the open set.\n";
    cout << "The queue is now: "; 
    GTestLabeled.printQindexes();
    
    GTestLabeled.set_status_decreaselbl(2,1,1);
    cout << "Vertex 2 has its label decreased.\n";
    cout << "The queue is now: "; 
    GTestLabeled.printQindexes();
}


// Test: Dijkstra algorithm.
void test_D()
{
    //Create graph
    Graph Gtest(9);
    
    //Gtest.addEdge(0,1,4);
    Gtest.addEdge(0,1,3);
    Gtest.addEdge(0,2,3);
    Gtest.addEdge(1,3,1);
    Gtest.addEdge(0,4,7);
    Gtest.addEdge(2,4,4);
    Gtest.addEdge(3,4,3);
    Gtest.addEdge(3,5,1);
    Gtest.addEdge(4,5,1);
    Gtest.addEdge(4,6,5);
    Gtest.addEdge(5,7,2);
    Gtest.addEdge(4,8,2);
    Gtest.addEdge(5,8,4);
    Gtest.addEdge(6,8,5);
    Gtest.addEdge(7,8,3);

    VertexLabeledGraph GTestLabeled(Gtest);

    unsigned int StartVtx, EndVtx;
    StartVtx = 4;
    EndVtx =0;
    
    cout << "-- Test Dijkstra's algorith: finding a path from vertex " << StartVtx << " to vertex " << EndVtx << endl; //To make analogous change in Main 6,7 .
    list<unsigned int> Path1=GTestLabeled.path(StartVtx,EndVtx,1);
    
    PrintList(Path1);

    cout << endl;
}



//Test: Execute Prim on a small graph.
void test_Prim_small()
{
    //Create graph
    Graph Gtest(9);
    
    Gtest.addEdge(0,1,4);
    Gtest.addEdge(0,2,3);
    Gtest.addEdge(1,3,1);
    Gtest.addEdge(0,4,7);
    Gtest.addEdge(2,4,4);
    Gtest.addEdge(3,4,3);
    Gtest.addEdge(3,5,1);
    Gtest.addEdge(4,5,1);
    Gtest.addEdge(4,6,5);
    Gtest.addEdge(5,7,2);
    Gtest.addEdge(4,8,2);
    Gtest.addEdge(5,8,4);
    Gtest.addEdge(6,8,5);
    Gtest.addEdge(7,8,3);

    VertexLabeledGraph GTestLabeled(Gtest);

    unsigned int StartVtx = 4;
    
    cout << "-- Test Prim's algorithm on small graph from the vertex " << StartVtx << " --" << endl;
    vector<pair<unsigned int, unsigned int>> Tree4=GTestLabeled.PrimMST(StartVtx,1);
    unsigned int CostTree= Tree4[Tree4.size()-1].second ;

    Tree4.pop_back();
    PrintPairVector(Tree4);
    cout << "Cost: " << CostTree << endl;
    
    cout << endl;
}


// Test 5: Execute Prim on the graph "SampleTestData".
// Print out a Minimum Spanning Tree and its cost (equal 30).
void test_Prim()
{
    ifstream infile("SampleTestData.txt");
    istream_iterator<unsigned int> start(infile);
    istream_iterator<unsigned int> end;
    
    Graph Gfile(start, end) ;
    
    VertexLabeledGraph GfileLabeled(Gfile);

    unsigned int StartVtx = 1;
    cout << "-- Testing Prim's algorithm from the vertex " << StartVtx << " --" << endl;

    vector<pair<unsigned int, unsigned int>> Tree0=GfileLabeled.PrimMST(StartVtx,1);
    //The 2nd argument of "PrimMST" controls the printout of the intermediate steps.
    //Change to 1 in order to print them,  change to 0 in order to hide.
    
    PrintPairVector(Tree0);
    cout << "Cost: " << Tree0[Tree0.size()-1].second << endl;
}


int main()
{
    test_create_graph();
    test_VertexLabeledGraph();
    test_D();
    test_Prim_small();
    test_Prim();
    return 0;
}