#ifndef _MAPS_H_
#define _MAPS_H_
#include <vector>

using namespace std;

//TASK:          COMPUTES THE POINT TO ELEMENT MAP (HASH TABLE)
// alorithm: loop over all elements and in each element loop over all nodes
// and then add that element number to the hash_table "p_to_e" for that node.
int create_p_to_e(int ntri, int tri[][3], vector<int> *p_to_e);

// TASK: COMPUTES THE ELEMENT CONTAINING THAT EDGE (EXCEPT THE ORIGINAL ELEMENT ITSELF)
// on success return the element number
// on failure returns -1
int elm_contain_edge(int pt1, int pt2, int itself, vector<int> *p_to_e);

//TASK: makes element surrounding element map or simply neighbors
void make_nbrs(int nn, int ntri, int tri[][3], int nbrs[][3], vector<int> *p_to_e);

// creates p2p (point to point map)
// alorithm: loop over triangles 
//              for each vertex add the global node number of other two vertices 
//           end
//           sort
//           remove repeated occurance of node numbers 
//           end
// sample : p2p[i=0..(nn-1)][0..number of connected nodes]
int create_p2p(int ntri, int **tri, int nn, vector<int> *p2p);


#endif
