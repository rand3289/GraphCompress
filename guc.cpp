#include <iostream>
#include <fstream>
#include <vector>
#include <string.h> // memset()
using namespace std;

typedef unsigned short int USI;

int main(int argc, char* argv[]) {
    if(argc != 2){
        cout << "usage: " << argv[0] << " filename" << endl;
        exit(1);
    }
    string filename = argv[1];

    ifstream inFile (filename, ios::in|ios::binary);
    if(!inFile.is_open()){
        cerr << "ERROR: failed to open " << filename << endl;
        exit(2);
    }

    ofstream outFile (filename+".out", ios::out|ios::binary);
    if(!outFile.is_open()){
        cerr << "ERROR: failed to open output file " << filename << endl;
        exit(3);
    }

    /********************************************************************/
    USI table[256][256];
    inFile.read( (char*) &table, sizeof(table) );

    size_t psize; // TODO: change it to unsigned int???
    inFile.read( (char*) &psize, sizeof(psize) );
    // TODO: use the psize while reading data !!!

    USI next;
    vector<USI> data;
    while( inFile.read( (char*) &next, sizeof(next)) ){
        data.push_back(next);
    }
    cout << "read " << 2*data.size() << " bytes from " << filename << endl;
    /*********************************************************************/

    USI* flat = reinterpret_cast<USI*>(table);
    for(int i=0; i<=0xFFFF; ++i){ // renormalize entries stored as 0
        flat[i] = flat[i] ? flat[i] : i;
    }

    USI inverted[256*256];
    memset(inverted, 0, sizeof(inverted));
    for(int i=0; i<255; ++i){ // i = cluster
        for(int j=0; j<256; ++j){
            inverted[ table[i+1][j] ] = i;
        }
    }

    int cluster = 0;
    for(USI& i: data){
        if(i>255){ // all nodes above 255 are absolute
            cluster = inverted[i];
        } else {
            i = table[cluster+1][i]; // substitute i with the real id of the node
        }
    }

    // undo reservations for the first 256 nodes
    for(int i=0; i<256; ++i){
        for(int j=1; j<data.size(); ++j){ // skip the first (head)
            if(table[0][i] == data[j]){
                data[j] = i;
            }
        }
    }

    for(USI j: data){ // not by reference!!!
        outFile.write( (char*) &j, sizeof(j) );
    }
    outFile.close();
}
