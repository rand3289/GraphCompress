// GCP (Graph ComPress) rearranges bytes in a high entrophy file so that it can be further compressed
// by any off-the-shelf compression utility like zip, 7zip, bzip, gzip, rar etc...
// After the file is compressed and extracted back it can be transformed into the original using GUC (Graph UnCompress)
//
// author:  toandrey(at)yahoo(dot)com
// licence: GPLv3
//
// Algorithm:
// A file is read in 16bit chunks to form unsigned short int.
// Each number represents a node in a graph with 0xFF00 nodes.
// The sequence of numbers read from a file represents a PATH through this graph.
// Nodes are arranged in 255 partitions with 256 nodes each.
// Optimization is done to rearrange the nodes among partitions to minimize the number of edges between partitions.
// When path is saved, edges local to one partition are stored as an index of that node in that partition
// This transforms the 16 bit short node ID to an 8 bit node index in the cluster (clusters have 256 nodes max)
// Even though indexes are stored as 16 bit, it introduces lots of zeros (from high byte=0) into the output file.
// Cluster index to node ID mappings are saved and the original path through the graph is saved.

#include <limits>
#include <algorithm>
#include <functional> // bind()
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
#include <chrono>
using namespace std::chrono;

//template <typename T>
//ostream& operator<<(ostream& os, vector<T> v){
//  for(T& t: v){ os << t << " "; }
//  return os;
//}


typedef unsigned short Node;
static const size_t NodeMax = std::numeric_limits<Node>::max() + 1; // max cardinal value Node can represent


class Graph{
    vector<vector<Node> > graph; // These are paths.  For each node being visited what is the next node?
    vector<size_t> counts; // counts stores the number of times each node is visited in the path
    vector<int> remapped; // first 256 nodes are reserved and need to be remapped to unused higher index nodes
    int mapped = 0;

    void countNodes(){ // fill vector counts with number of times every node is visited
        for(auto& v: graph){
            for(Node n: v){
                ++counts[n];
            }
        }
    }

    void replaceConnection(const Node old1, const Node new1){ // while remapping nodes to new IDs, replace edges pointing to them
        for(auto& g: graph){
            for(Node& c: g){
                if(c == old1){
                    c = new1;
                }
            }
        }
    }

    void remap(){ // if any of the lower 0xFF nodes have connections, we need to remap them to higher IDs
        int j = 256;
        for(int i = 0; i<= 255; ++i){ // for all nodes with id 0 through 255
            if( -1 == remapped[i] ){ continue; }  // it is not in use
            while( -1 != remapped[j] ){ ++j; } // find next unused node ID with ID > 255
//            cout << i << "->" << j << " ";
	    graph[i].swap(graph[j]); // vectors of connections
            remapped[i] = j;
            remapped[j] = i;
            replaceConnection(i,j);
        }
        int cnt = 256 - std::count(remapped.begin(), remapped.end(), -1);
        cout << cnt << " nodes were remapped in the path(graph)." << endl;
//        cout << "first -1 node: " << &*find(remapped.begin(), remapped.end(), -1) - &remapped[0];
//        cout << ", last -1 node: " << &*find(remapped.rbegin(), remapped.rend(), -1) - &remapped[0] << endl;
    }

public:
    Graph(): graph(NodeMax), counts(NodeMax,0), remapped(NodeMax,-1) {}
    const vector<Node>& getConnections(Node node) const { return graph[node]; }
    bool isPopular(Node node) const { return counts[node] > 2; }
    void finalize(){ remap(); countNodes(); } // finalize Graph's internal structures

    bool addPath(Node from, Node to) {
        if( -1 == remapped[from] ){
            remapped[from] = from;
            ++mapped; // count how many nodes are in use
        }
        if( -1 == remapped[to] ){ // if node does not have any paths through it yet
            remapped[to] = to;
            ++mapped; // count how many nodes are in use
        }
        graph[from].push_back(to);
        return mapped < NodeMax-256; // reserve lower 256 (0 through 255) nodes
    }

    bool followPath(Node& dataInOut){ // TODO: it's a bit tricky how it returns data... change it.
        auto& path = graph[dataInOut];
        if(path.empty()) {
            return false;
        }
        dataInOut = path.front();
        path.erase(path.begin()); // inefficient - use dequeue???
        return true;
    }

    ostream& serialize(ostream& os) const {
        for(int i=0; i<256; ++i){ // first 256 nodes were remapped to reserve ids [0-255]
            Node node = (-1 == remapped[i]) ? 0 : remapped[i];
            os.write( (char*)&node, sizeof(node) );
        }
        return os;
    }
};


class Partition;
class Node2part { // provides a lookup of Partition* by Node
    vector<Partition*> lookup;
public:
    Node2part():lookup(NodeMax, nullptr) {}
    void setPart(Partition* part, Node node){ lookup[node] = part; }
    Partition* getPart(Node node){ return lookup[node]; }
    void swapParts(Node node1, Node node2){ std::swap(lookup[node1], lookup[node2]); }
    Node useOrReplace(Node from, Node to) const;
    long long int countAllHyperEdges(const Graph& graph) const;
};


class Partition{
    vector<Node> local;

    Node& pickRandom(const Graph& graph){ // pick a popular random node from this partition
        const int idx = rand() % local.size();
        Node* pnode = &local[ idx ];                               // pick a random vertex in this partition
        for(int i=1; !graph.isPopular(*pnode) && i<local.size(); ++i){ // try to avoid vertices that occur less than 3 times in the path
            pnode = &local[ (idx+i)%local.size() ];
        }
        return *pnode;
    }

public:

    bool addNode(Node node){  // add node to partition
        if(local.size() >= 256){ return false; }
        local.push_back(node);
        return true;
    }

    int indexOf(Node node) const { // returns -1 if nodeID is not found
        return &*find(local.rbegin(), local.rend(), node) - &local[0]; // in reverse so that &*rend() is just before &local[0]
    }

    ostream& serialize(ostream& os, int partIdx) const {
        for(int j=0; j<local.size(); ++j){
            Node n = local[j];
            if(n == (partIdx+1)*256+j){ // if node has not been swapped out of it's original place when partitions were initialized
                n = 0;                  // write zero to file because we can recreate it in decompressor
            }
            os.write( (char*) &n, sizeof(n) );
        }
        return os;
    }

    // There are nodes in partitions that were placed there by default then swapped out and later placed back at a different spot.
    // Put them back into their original places so they can be stored as a 0 in the file during serialize().
    int order(const int partIdx){
        int reordered = 0;
        for(int i=0; i < local.size(); ++i){
            Node& current = local[i];
            if(partIdx == current/256 &&  i != current%256 ){ // if it belongs to this partition but not in it's place
                std::swap(current, local[current%256]);
                ++reordered;
            }
        }
        return reordered;
    }

    void optimize(vector<Partition>& parts, Node2part& node2part, const Graph& graph){
        Partition* part2 = this;
        while(part2==this){
            part2 = &parts[ rand() % parts.size() ];
        }

        Node& n1 = pickRandom(graph);
        Node& n2 = part2->pickRandom(graph); // this is a reference to a node in part2 !!!

        long long int edgeCount = node2part.countAllHyperEdges(graph);
        node2part.swapParts(n1, n2);
        long long int newEdgeCount = node2part.countAllHyperEdges(graph);

        if(edgeCount <= newEdgeCount){ // new count is worse or same.  undo changes (stable order)
            node2part.swapParts(n1, n2); // swap them back!
        } else {
            std::swap(n1, n2);
            cout << '.';
            cout.flush();
        }
    }
}; // Partition


#include <cassert>
Node Node2part::useOrReplace(Node from, Node to) const {
    if( lookup[from] == lookup[to] ){      // are nodes "from" and "to" in the same partition???
        int idx = lookup[to]->indexOf(to); // since partitions are less than 256 entries long, this returns a "short ID"
        assert(idx!=-1);                   // lookup is out of date with partition!!!
        return idx;
    }
    return to;
}

long long int Node2part::countAllHyperEdges(const Graph& graph) const {
    long long int count = 0;
    for(int i=0; i < NodeMax; ++i){
        const Partition* partition = lookup[i]; // partition for node i
        for(Node conn: graph.getConnections(i) ){      // for each node i is connected to
            if( partition != lookup[conn] ){    // if partitions are different
                ++count;
            }
        }
    }
    return count;
}


int kbhit(){ // from https://www.flipcode.com/archives/_kbhit_for_Linux.shtml
    const int STDIN = 0;
    fd_set rdset;
    FD_ZERO(&rdset);
    FD_SET(STDIN, &rdset);

    timeval timeout;
    timeout.tv_sec  = 0;
    timeout.tv_usec = 0;

    return select(STDIN + 1, &rdset, NULL, NULL, &timeout);
}


int main(int argc, char* argv[]){
    if(argc != 2){
        cout << "Graph Compress - file sparcing utility (entrophy shifter)" << endl << "usage: " << argv[0] << " filename" << endl;
        exit(1);
    }
    string filename = argv[1];
    ifstream inFile (filename, ios::in|ios::binary);
    if(!inFile.is_open()){
        cerr << "ERROR: failed to open " << filename << endl;
        exit(2);
    }

    string outFName = filename + ".comp";
    ofstream outFile (outFName, ios::out|ios::binary);  // open output file now so that error does not happen 10 hours later
    if(!outFile.is_open()){
        cerr << "ERROR: failed to open output file " << outFName << endl;
        exit(3);
    }

/******************************** done with files ******************************/
    Node head = 0;                   // treat file as a graph path and
    inFile.read( (char*) &head, sizeof(head)); // read the first node in the path
    size_t count = inFile.gcount();
    Node node = head;
    Node next = 0;

    Graph graph;
    for(; inFile.read( (char*) &next, sizeof(next)); count+= inFile.gcount() ){
        if( !graph.addPath(node,next) ){
            cout << "filled the path with " << count/2-1 << " edges." << endl; // two bytes per node - first node
            inFile.unget(); inFile.unget(); // we could not use this data, put it back into stream
            break;
        }
        node = next;
        next = 0; // in case only 1 byte will be read into next
    }

    cout << "according to tellg(), read " << inFile.tellg() << " bytes from " << filename << endl;
    inFile.close();
    graph.finalize();
//TODO: what happens after remap() if head is one of the first 256 nodes...  FIX THIS here or in guc???
//TODO: this just uses the data read so far and throws the rest away... FIX THIS
//TODO: I am reading 2 bytes at a time... what if there is only 1 byte left.. FIX THIS
//TODO: No need to use count += inFile.gcount() ... use count = tellg() instead???

    Node2part node2part;     // parts map partition->nodes.  Node2part maps node->partition.
    vector<Partition> parts(255); // at first nodes are distributed in-order among partitions
    Partition* partPtr = &parts[0];
    for(int i = 256; i < NodeMax; ++i){ // for all nodes (lower 256 [0-255] nodes are reserved)
        while( ! partPtr->addNode(i) ){  // if we can't add it to the last partition, add it to a new one
            ++partPtr;
        }
        node2part.setPart(partPtr,i);
    }

    const long long int hyper = node2part.countAllHyperEdges(graph);
    cout << count/2-1 << " edges formed (" << count << " bytes read) from " << filename << endl;
    cout << "allocated nodes into " << parts.size() << " partitions." << endl;
    cout << hyper << " hyper edges found." << endl;
    cout << "Compressing. Press enter to stop and save." << endl;

    auto start = system_clock::now();
    long long optCount = 0;
    while(!kbhit()){
        for(auto& p: parts){ // let each partition to optimize itself
            p.optimize(parts, node2part, graph);
        }
        optCount += parts.size();
    }

    int runtime = duration_cast<seconds>(system_clock::now() - start).count();
    double perSec = ((double) optCount) / (double)runtime;
    cout << "Performed "<< optCount << " optimizations in " << runtime << " seconds.  " << perSec << " per second." << endl;
    const long long int hyper2 = node2part.countAllHyperEdges(graph);
    cout << hyper2 << " hyper edges remaining. " << (hyper-hyper2) << " less." << endl;

/*************************************** save output  ***********************************/
    cout << "Saving graph data." << endl;
    graph.serialize(outFile); // write out first 256 remapped nodes

    cout << "Reordering nodes." << endl;
    int reordered = 0;
    for(int i=0; i<parts.size(); ++i){
        reordered += parts[i].order(i);
    }
    cout << reordered << " nodes reordered." << endl;

    cout << "Saving partition data." << endl;
    for(int i=0; i< parts.size(); ++i){
        parts[i].serialize(outFile, i);
    }

    cout << "Saving path data." << endl;
    outFile.write( (char*) &count, sizeof(count)); // save the length of the path
    outFile.write( (char*) &head,  sizeof(head));  // fist node in the path (start)

    Node prev = head;
    Node data = head;
    while ( graph.followPath(data) ){
        Node n = node2part.useOrReplace(prev, data); // Node n = data; // testing ONLY !!!
        outFile.write( (char*) &n, sizeof(n) );
        prev = data;
    }
    outFile.close();
}
