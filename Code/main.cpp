#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <climits>
#include <cstring>
#include <vector>
#include <time.h>
#define num_core 4
#define V 10000
//#define debug 1
#define showEdge

#define dropRate 0.7
using namespace std;

/*
 * Indexed min priority queue
 * Supports insertion in O(log N), deletion of any key (regardless of whether
 * the key is the minimum key or not) in O(log N) and changes to key values
 * in O(log N), where N is the number of
 * elements currently in the PQ
 */
class MinIndexedPQ {
    int NMAX, N, *heap, *index, *keys;
    
    void swap(int i, int j) {
        int t = heap[i]; heap[i] = heap[j]; heap[j] = t;
        index[heap[i]] = i; index[heap[j]] = j;
    }
    
    void bubbleUp(int k)    {
        while(k > 1 && keys[heap[k/2]] > keys[heap[k]])   {
            swap(k, k/2);
            k = k/2;
        }
    }
    
    void bubbleDown(int k)  {
        int j;
        while(2*k <= N) {
            j = 2*k;
            if(j < N && keys[heap[j]] > keys[heap[j+1]])
                j++;
            if(keys[heap[k]] <= keys[heap[j]])
                break;
            swap(k, j);
            k = j;
        }
    }
    
public:
    // Create an empty MinIndexedPQ which can contain atmost NMAX elements
    MinIndexedPQ(int NMAX)  {
        this->NMAX = NMAX;
        N = 0;
        keys = new int[NMAX + 1];
        heap = new int[NMAX + 1];
        index = new int[NMAX + 1];
        for(int i = 0; i <= NMAX; i++)
            index[i] = -1;
    }
    MinIndexedPQ()  {
        this->NMAX = 0;
        N = 0;
        keys = NULL;
        heap = NULL;
        index = NULL;
    }
    MinIndexedPQ(MinIndexedPQ &obj)  {
        this->NMAX = obj.NMAX;
        this->N = obj.N;
        this->keys = new int[NMAX + 1];
        this->heap = new int[NMAX + 1];
        this->index = new int[NMAX + 1];
        std::memcpy(this->keys, obj.keys, sizeof(int)*(NMAX + 1));
        std::memcpy(this->heap, obj.heap, sizeof(int)*(NMAX + 1));
        std::memcpy(this->index, obj.index, sizeof(int)*(NMAX + 1));
    }
    
    ~MinIndexedPQ() {
        delete [] keys;
        delete [] heap;
        delete [] index;
    }
    
    void clear(){
        this->N = 0;
        for(int i = 0; i <= NMAX; i++)
            index[i] = -1;
    }
    
    void CopyPQ(MinIndexedPQ &obj)  {
        this->NMAX = obj.NMAX;
        this->N = obj.N;
        this->keys = new int[NMAX + 1];
        this->heap = new int[NMAX + 1];
        this->index = new int[NMAX + 1];
        std::memcpy(this->keys, obj.keys, sizeof(int)*(NMAX + 1));
        std::memcpy(this->heap, obj.heap, sizeof(int)*(NMAX + 1));
        std::memcpy(this->index, obj.index, sizeof(int)*(NMAX + 1));
    }
    
    // check if the PQ is empty
    bool isEmpty()  {
        return N == 0;
    }
    
    // check if i is an index on the PQ
    bool contains(int i)    {
        return (index[i] != -1) && (index[i]!=-2);
    }
    
    bool has_deleted(int i){
        return index[i] == -2;
    }
    
    void check_validity(int i){
        if(i<0 || i >= NMAX){
            printf("-------%d-------", i);
        }
        if(contains(i)){
            printf("contain_error");
        }
        if(N>NMAX){
            printf("over");
        }
    }
    // associate key with index i; 0 < i < NMAX
    void insert(int i, int key) {
        //check_validity(i);
        N++;
        index[i] = N;
        heap[N] = i;
        keys[i] = key;
        bubbleUp(N);
    }
    
    
    // delete the minimal key and return its associated index
    // Warning: Don't try to read from this index after calling this function
    int deleteMin() {
        int min = heap[1];
        swap(1, N--);
        bubbleDown(1);
        index[min] = -2;
        heap[N+1] = -1;
        return min;
    }
    int popMin(){
        return heap[1];
    }
    int popKey(){
        int min = heap[1];
        return keys[min];
    }
    
    // decrease the key associated with index i to the specified value
    void decreaseKey(int i, int key)    {
        keys[i] = key;
        bubbleUp(index[i]);
    }
    
    // merge the PQ with another PQ
    void merge(MinIndexedPQ pq2){
        int n2 = pq2.N;
        int i, edg;
        for(i=0;i<n2;i++){
            edg = pq2.deleteMin();
            if(this->contains(edg)){
                if(pq2.keys[edg] < this->keys[edg]){
                    this->keys[edg]=pq2.keys[edg];
                    bubbleUp(this->index[edg]);
                }
            }
            else if (this->has_deleted(edg)){
                continue;
            }
            else{
                insert(edg, pq2.keys[edg]);
            }
        }
    }
};


// representation of directed edge to vertex 'to' having weight 'weight'
struct Edge {
    int to, weight;
};
struct linkage {
    int n1, n2, dist;
};

typedef vector <Edge*> adjList;



void run_test()  {
    int i, u, v, j, cost, num_e, num_visit;
    double start_time, end_time;
    
    Edge *tmp;
    adjList *G;
    omp_set_dynamic(0);
    omp_set_num_threads( num_core );
    // Assuming vertices are labeled 0...V-1
    // freopen("input1.txt", "r", stdin);
    num_e=0;
    FILE *fp = fopen ("1wc.txt","r");
    if (fp == 0) {
        printf("fopen failed. ");
        return;
    }
    G = new adjList[V];
       // Enter E undirected edges (u, v, weight)
    while(fscanf(fp, "%d %d %d", &u, &v, &cost)!=EOF){
        //printf("%d %d %d\n", u, v, cost);
        tmp = new Edge;
        num_e++;
        tmp->to = v; tmp->weight = cost;
        G[u].push_back(tmp);    // add edge to adjacency list of u
        tmp = new Edge;
        tmp->to = u; tmp->weight = cost;
        G[v].push_back(tmp);    // add edge to adjacency list of v
    }
    fclose(fp);
    start_time=omp_get_wtime();
    MinIndexedPQ PQ_array[num_core];
    int dist_array[num_core][V];
    bool marked_array[num_core][V];
    int marked_by[V];
    for(i=0;i<V;i++){
        marked_by[i]=-1;
    }
    int edge_to_array[num_core][V];
    int merge_flag[num_core];
    for(i=0;i<num_core;i++){
        merge_flag[i]=0;
    }
    int applicant[num_core];
    int final_result;
    linkage links[V];
    bool end_thread;
    num_visit=0;
    int *dist;
    int *edgeTo;
    bool *marked;
    
    
#pragma omp parallel shared(marked_array, merge_flag, G, PQ_array, dist_array, edge_to_array, applicant, final_result, num_e, num_visit, marked_by, links, start_time, end_time) private( u,v, i, end_thread, j, dist, edgeTo, marked)
    {
        int thread_id = omp_get_thread_num();
        int src = (V/num_core)*thread_id; // Taking a vertex as source
        end_thread=false;
        #pragma omp critical(mark_visit)
        {
            while(marked_by[src]!=-1 && V-num_visit-1>dropRate*V){
                src=rand() % V;
            }
            marked_by[src]=thread_id;
        }
        dist=new int[V];
        edgeTo=new int[V];
        marked=new bool[V];
        edgeTo[src]=-1;
        for(i=0; i<V; i++){
            dist[i] = INT_MAX;
            marked[i]=false;
        }
        dist[src] = 0;
        //memset(marked, false, V*sizeof(bool));
        
        MinIndexedPQ pq(V);
        pq.insert(src, 0);
        
        while(!pq.isEmpty() && num_visit<V-1)    {
            //first check if anyone wants to merge: 0 means no request; 1 means data uploading; 2 means upload finished
            if(merge_flag[thread_id]==2){
                //download data to the local data structure
                int uploader=applicant[thread_id];
                for(i=0; i<V; i++){
                    if(PQ_array[uploader].contains(i) && dist[i]>dist_array[uploader][i]){
                        dist[i]=dist_array[uploader][i];
                        edgeTo[i]=edge_to_array[uploader][i];
                    }
                    marked[i] = marked[i] || marked_array[uploader][i];
                }
                pq.merge(PQ_array[uploader]);
                // mark the thread as ready to accept new merges
                #pragma omp atomic
                merge_flag[thread_id]*=0;
                #ifdef debug
                printf("thread %d receives thread %d\n", thread_id, uploader);
                #endif
            }
            else if (merge_flag[thread_id]!=0){
                #ifdef debug
                printf("Thread %d: The merge flag is %d when the add edge loop start, wait for next iteration\n", thread_id, merge_flag[thread_id]);
                #endif
            }
            
            
            u = pq.popMin();
            bool surrender = true;
            int aim;
            
            //if the node visited by the thread but unvisited by merged in thread, skip the node
            if((marked[u]==true || marked_by[u]==thread_id)&&u!=src){
                u=pq.deleteMin();
                continue;
            }
            //if the new vertice has been visited by another thread
            //start merging: copy all data to chared memory, set flag, and clear all data to start at a new point
            if(marked_by[u]!=-1 && u!=src && merge_flag[thread_id]==0){
                #pragma omp critical(mark_visit)
                {
                    aim = marked_by[u];
                    if(merge_flag[thread_id]==0 && merge_flag[aim]==0 && aim!=thread_id){
                        merge_flag[aim]=1;
                        applicant[aim]=thread_id;
                        u=pq.deleteMin();
                        for(i=0;i<V;i++){
                            if(marked_by[i]==thread_id){
                                marked_by[i]=aim;
                            }
                        }
                        surrender=false;
                    }
                }
                if(surrender){
                    continue;
                }
                //prepare to merge into another thread
                #ifdef debug
                printf("thread %d merging thread %d on %d\n", thread_id, aim, u);
                #endif
                
                std::memcpy(&dist_array[thread_id], dist, sizeof(int)*V);
                std::memcpy(&marked_array[thread_id], marked, sizeof(bool)*V);
                std::memcpy(&edge_to_array[thread_id], edgeTo, sizeof(int)*V);
                PQ_array[thread_id].CopyPQ(pq);
                #pragma omp atomic
                merge_flag[marked_by[u]]++;
                
                
                #pragma omp atomic capture
                {
                    j=num_visit;
                    num_visit++;
                }
                links[j].n1=u;
                links[j].n2=edgeTo[u];
                links[j].dist=dist[u];
                #ifdef debug
                printf("add edge %d %d %d by thread %d #%d\n", u, edgeTo[u], dist[u], thread_id, j);
                #endif
                
                //if too close to finish, end the thread
                if(V - num_visit - 1 < dropRate*V){
                    #ifdef debug
                    printf("thread %d ends\n", thread_id);
                    #endif
                    end_thread=true;
                    break;
                }
                //reset all data structure
                
                src = (V/num_core)*thread_id; // Taking a vertex as source
                int ct=0;
                while(marked_by[src]!=-1){
                    src=rand() % V;
                    ct++;
                    if(ct>99){
                        if(V - num_visit - 1 < dropRate*V){
                            end_thread=true;
                            break;
                        }else{
                            ct=0;
                        }
                    }
                }
                if(end_thread){
                    break;
                }
                #pragma omp critical(mark_visit)
                {
                    while(marked_by[src]!=-1 && V - num_visit - 1 < dropRate*V){
                        src=rand() % V;
                    }
                    marked_by[src]=thread_id;
                }
                edgeTo[src]=-1;
                for(i=1; i<V; i++){
                    dist[i] = INT_MAX;
                    marked[i]=false;
                }
                dist[src] = 0;
                pq.clear();
                pq.insert(src, 0);
                continue;
            }
            
            //if the new vertice has not been visited, then do the loop for adding edges
            else{
                surrender=true;
                #pragma omp critical(mark_visit)
                {
                    if(marked_by[u]==-1){
                        marked_by[u]=thread_id;
                        surrender=false;
                    }
                }
                if(surrender && u!=src){
                    continue;
                }
                u=pq.deleteMin();
                #ifdef debug
                printf("node %d visited by thread %d\n", u, thread_id);
                #endif
                marked[u] = true;
                if(src!=u){
                    #pragma omp atomic capture
                    {
                        j=num_visit;
                        num_visit++;
                    }
                    links[j].n1=u;
                    links[j].n2=edgeTo[u];
                    links[j].dist=dist[u];
                    #ifdef debug
                    printf("normal add edge %d %d %d by thread %d #%d\n", u, edgeTo[u], dist[u], thread_id, j);
                    #endif
                }
        
                //add all edges that links to the visited vertice to the indexable PQ
                for(adjList::iterator it = G[u].begin(); it != G[u].end(); it++)    {
                    v = (*it)->to;
                    if(marked[v]) continue;
                    if((*it)->weight < dist[v])    {
                        dist[v] = (*it)->weight;
                        edgeTo[v] = u;
                        if(pq.contains(v)) {
                            pq.decreaseKey(v, dist[v]);
                        }
                        else {
                            pq.insert(v, dist[v]);
                        }
                    }
                }
            }
        }
        
        
        if(!end_thread){
            end_time=omp_get_wtime();
            #ifdef showEdge
            printf("Edges in MST (by thread %d):\n", thread_id);
            cost=0;
            for(i=0, cost=0; i<V-1; i++)  {
                cost += links[i].dist;
                printf("%d %d %d  #%d\n", links[i].n1, links[i].n2, links[i].dist, i);
            }
            #endif
            printf("Total cost of MST : %d\n", cost);
            printf("time=%f dropRate=%f\n", end_time - start_time, dropRate);
        }
        delete [] marked;
        delete [] dist;
        delete [] edgeTo;
    }
    
    // print all the edges in the MST
    
    
    return;
}


// function to

int main(){
    run_test();
    return 0;
}
