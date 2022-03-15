#include <vector>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <Rcpp.h>

using namespace std;

// [[Rcpp::plugins(cpp14)]]

double distance3d(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2){
    return sqrt(pow(z2-z1,2)+pow(y2-y1,2)+pow(x2-x2,2));
}

// [[Rcpp::export]]
vector<vector<double> > disma3d(vector<double> &w, vector<double> &dist, vector<double> &sim, int &N){

    vector<vector<double> > distance(N, vector<double>(N,0));
    for (int i=0; i < N; i++){
        for (int j=0; j < i; j++){
            distance[j][i] = distance3d(w[i], dist[i], sim[i], w[j], dist[j], sim[j]);
        }
    }
    return distance;
}

double distance2d(double x1, double y1, double x2, double y2){
    return sqrt(pow(y2-y1,2)+pow(x2-x1,2));
}

// [[Rcpp::export]]
vector<vector<double> > disma2d_truncated(vector<double> &dist_train,
                                        vector<double> &sim_train,
                                        vector<double> &dist_test,
                                        vector<double> &sim_test,
                                        int &N_train,
                                        int &N_test){

    vector<vector<double> > distance_train_test(N_test, vector<double>(N_train,0));
    for (int i=0; i < N_test; i++){
        for (int j=0; j < N_train; j++){
            distance_train_test[i][j] = distance2d(dist_train[j], sim_train[j], dist_test[i], sim_test[i]);
        }
    }

    return distance_train_test;
}

vector<vector<double> > K_centroid(vector<double> &dist, vector<double> &sim, vector<int> &id, int N, int K){
    vector<vector<double> > K_coords(K, vector<double>(2,0));
    int count_k;
    for (int k=0; k < K; k++){
        count_k=0;
        for (int i=0; i < N; i++){
            if (id[i] == k+1){
                K_coords[k][0] += dist[i];
                K_coords[k][1] += sim[i];
                count_k++;
            }
        }
        K_coords[k][0] /= count_k;
        K_coords[k][1] /= count_k;
    }
    return K_coords;
}

// [[Rcpp::export]]
vector<vector<double> > disma2d_centroid(vector<double> &dist_train,
                                vector<double> &sim_train,
                                vector<double> &dist_test,
                                vector<double> &sim_test,
                                vector<int> &id,
                                int &N_train,
                                int &N_test,
                                int &K){

    vector<vector<double> > K_coords = K_centroid(dist_train, sim_train, id, N_train, K);

    vector<vector<double> > distance_train_test(N_test, vector<double>(K,0));
    for (int i=0; i < N_test; i++){
        for (int k=0; k < K; k++){
            distance_train_test[i][k] = distance2d(K_coords[k][0], K_coords[k][1], dist_test[i], sim_test[i]);

        }
    }

    return distance_train_test;
}

// [[Rcpp::export]]
Rcpp::List distma2centroid(vector<double> &dist_train,
                                vector<double> &sim_train,
                                vector<double> &dist_test,
                                vector<double> &sim_test,
                                vector<int> &id_train,
                                vector<int> &id_test,
                                int &N_train,
                                int &N_test,
                                int &K){

    vector<vector<double> > K_coords_train = K_centroid(dist_train, sim_train, id_train, N_train, K);
    vector<vector<double> > K_coords_test = K_centroid(dist_test, sim_test, id_test, N_test, K);
    vector<vector<double> > K_coords(K, vector<double>(2,0));
    
    for (int k=0; k < K; k++){
        for (int j=0; j < 2; j++){
            K_coords[k][j] = (K_coords_train[k][j] + K_coords_test[k][j])/2;
        }
    }
    vector<vector<double> > leave2centroid(2*K, vector<double>(max(N_test, N_train),0));

    for (int k=0; k < K; k++){
        for (int i=0; i < N_train; i++){
            leave2centroid[k][i] = distance2d(K_coords[k][0], K_coords[k][1], dist_train[i], sim_train[i]);
        }
    }
    for (int k=0; k < K; k++){
        for (int i=0; i < N_test; i++){
            leave2centroid[k+K][i] = distance2d(K_coords[k][0], K_coords[k][1], dist_test[i], sim_test[i]);
        }
    }
    return Rcpp::List::create(leave2centroid, K_coords);
}

// [[Rcpp::export]]
Rcpp::List distma2zero(vector<double> &dist_train,
                                vector<double> &sim_train,
                                vector<double> &dist_test,
                                vector<double> &sim_test,
                                vector<int> &id_train,
                                int &N_train,
                                int &N_test,
                                int &K){

    vector<vector<double> > K_coords = K_centroid(dist_train, sim_train, id_train, N_train, K);
    vector<vector<double> > leave2centroid(2*K, vector<double>(max(N_test, N_train),0));

    for (int k=0; k < K; k++){
        for (int i=0; i < N_train; i++){
            leave2centroid[k][i] = distance2d(K_coords[k][0], K_coords[k][1], dist_train[i], sim_train[i]);
        }
    }
    for (int k=0; k < K; k++){
        for (int i=0; i < N_test; i++){
            leave2centroid[k+K][i] = distance2d(K_coords[k][0], K_coords[k][1], dist_test[i], sim_test[i]);
        }
    }
    return Rcpp::List::create(leave2centroid, K_coords);
}