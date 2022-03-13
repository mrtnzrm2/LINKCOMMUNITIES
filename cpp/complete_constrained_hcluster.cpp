#include <vector>
#include <math.h>
#include <Rcpp.h>
#include <numeric>

using namespace std;

// [[Rcpp::plugins(cpp14)]]

vector<int> find_member(vector<int> &id, int &N, int k){
    vector<int> members_k;
    for (int i=0; i < N; i++){
        if (id[i] == k)
            members_k.push_back(i);
    }
    return members_k;
}

double complete_linkage(vector<vector<double> > &distance, vector<int> &members_k, int node){
    double max_distance = -1, dist1, dist2;
    int mem_size = members_k.size();

    for (int i = 0; i <= mem_size/2; i++){
        if (node > members_k[i])
            dist1 = distance[node][members_k[i]];
        else
            dist1 = distance[members_k[i]][node];
        if (node > members_k[mem_size-i-1])
            dist2 = distance[node][members_k[mem_size-i-1]];
        else
            dist2 = distance[members_k[mem_size-i-1]][node];
        if (dist1 > dist2){
            if (dist1 > max_distance)
                max_distance = dist1;
        } else{
             if (dist2 > max_distance)
                max_distance = dist2;
        }
            
    }
    return max_distance;
}

// [[Rcpp::export]]
vector<vector<int> > comctr_hclust(vector<vector<double> > &distance, vector<int> &id, int &N_train, int &N_test, int &K){

    vector<vector<int> > communities_train;
    for (int i=0; i < K; i++){
        communities_train.push_back(find_member(id, N_train, i+1));
    }

    vector<vector<int> > communities_test(K, vector<int>(N_test,0));
    vector<int> test_ids(N_test);
    iota(test_ids.begin(), test_ids.end(), 0);
    double min_dist, com_dist;
    int com, node, count=0;
    
    while (test_ids.size() >= 1){
        if (int(count*100/N_test)%10 == 0)
            Rcpp::Rcout << floor(count*100/N_test) << "%" << " ";
        min_dist = INT32_MAX;
        for (int j=0; j < K; j++){
            for (int k=0; k < test_ids.size(); k++){
                com_dist = complete_linkage(distance, communities_train[j], test_ids[k]+N_train);
                if (com_dist < min_dist){
                    com = j;
                    node = k;
                    min_dist = com_dist;
                }
            }
        }
        communities_train[com].push_back(test_ids[node]+N_train);
        communities_test[com][test_ids[node]]++;
        test_ids.erase(test_ids.begin()+node);
        count++;
    }
    Rcpp::Rcout << endl;
    return communities_test;
}
