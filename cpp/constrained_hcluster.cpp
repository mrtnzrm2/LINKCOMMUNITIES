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
        dist1 = distance[node][members_k[i]];
        dist2 = distance[node][members_k[mem_size-i-1]];
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

double single_linkage(vector<vector<double> > &distance, vector<int> &members_k, int node){
    double min_distance = INT32_MAX, dist1, dist2;
    int mem_size = members_k.size();

    for (int i = 0; i <= mem_size/2; i++){
        dist1 = distance[node][members_k[i]];
        dist2 = distance[node][members_k[mem_size-i-1]];
        if (dist1 > dist2){
            if (dist2 < min_distance)
                min_distance = dist2;
        } else{
             if (dist1 < min_distance)
                min_distance = dist1;
        }
            
    }
    return min_distance;
}

double average_linkage(vector<vector<double> > &distance, vector<int> &members_k, int node){
    double ave_distance = 0;
    int mem_size = members_k.size();

    for (int i = 0; i < mem_size; i++)
            ave_distance += distance[node][members_k[i]];
    if (mem_size > 0)
        ave_distance /= mem_size;
    else
        ave_distance = 0;
    return ave_distance;
}

int print_progress(int count, int N, int flag){
    if (int(count*100/N)%10 == 0 && flag != int(count*100/N)){
             Rcpp::Rcout << floor(count*100/N) << "%" << " ";
             flag = int(count*100/N);
    }
    return flag;
}

// [[Rcpp::export]]
vector<vector<int> > constrained_hclust(vector<vector<double> > &distance, vector<int> &id, int &N_train, int &N_test, int K){
    vector<vector<int> > communities_train;
    for (int i=0; i < K; i++){
        communities_train.push_back(find_member(id, N_train, i+1));
    }

    vector<vector<int> > communities_test(K, vector<int>(N_test,0));
    double min_dist, dist;
    int com, pre=-1;
    
    for (int i=0; i < N_test; i++){
        pre = print_progress(i, N_test, pre);
        min_dist = INT32_MAX;
        for (int k=0; k < K; k++){
            dist = complete_linkage(distance, communities_train[k], i);
            if (dist < min_dist){
                com = k;
                min_dist = dist;
            }
        }
        communities_test[com][i]++;
    }
    Rcpp::Rcout << endl;
    return communities_test;
}

// [[Rcpp::export]]
vector<vector<int> > centroid_hclust(vector<vector<double> > &distance, vector<int> &id, int &N_train, int &N_test, int &K){

    vector<vector<int> > communities_test(K, vector<int>(N_test,0));
    double min_dist, dist;
    int com, pre=-1;
    
    for (int i=0; i < N_test; i++){
        pre = print_progress(i, N_test, pre);
        min_dist = INT32_MAX;
        for (int k=0; k < K; k++){
            dist = distance[i][k];
            if (dist < min_dist){
                com = k;
                min_dist = dist;
            }
        }
        communities_test[com][i]++;
    }
    Rcpp::Rcout << endl;
    return communities_test;
}
