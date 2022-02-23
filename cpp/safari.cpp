#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h> 

using namespace std;

template <typename S>
ostream& operator<<(ostream &os, const vector<S> &vector){
    // Printing all the elements
    // using <<
    for (auto element : vector) {
        os << element << " ";
    }
    return os;
}
template<typename Ra, typename Rb>
vector<vector<int_fast32_t> > nodes_parallelization(Ra &N, Rb &chunks){

    vector<vector<int_fast32_t> > content(2, vector<int_fast32_t>(chunks, 0));

    int_fast32_t sep = (int) N/chunks;

    content[0][0] = 0;
    content[1][0] = sep;

    for (int i=0; i < chunks - 1; i++){
        content[0][i+1] += sep + content[0][i];
        content[1][i+1] += sep + content[1][i];
        if (i == chunks - 2)
            content[1][i+1] = N-1;
    }

    return content;
}

int main(){

    int N = 4;
    double leaves = 5;

    vector<vector<int_fast32_t> > range = nodes_parallelization(leaves, N);

    for (int i=0; i < N; i++)
        printf("%i %i\n", range[0][i], range[1][i]);

    return 0;
}