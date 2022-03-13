#include <vector>
#include <math.h>
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
    return sqrt(pow(y2-y1,2)+pow(x2-x2,2));
}

// [[Rcpp::export]]
vector<vector<double> > disma2d(vector<double> &dist_train,
                                vector<double> &sim_train,
                                vector<double> &dist_test,
                                vector<double> &sim_test,
                                int &N_train,
                                int &N_test){

    vector<vector<double> > distance_train_train(N_train, vector<double>(N_train,0));
    for (int i=0; i < N_train; i++){
        for (int j=0; j < i; j++){
            distance_train_train[i][j] = distance2d(dist_train[i], sim_train[i], dist_train[j], sim_train[j]);
        }
    }

    vector<vector<double> > distance_test_test(N_test, vector<double>(N_test,0));
    for (int i=0; i < N_test; i++){
        for (int j=0; j < i; j++){
            distance_test_test[i][j] = distance2d(dist_test[i], sim_test[i], dist_test[j], sim_test[j]);
        }
    }
    vector<vector<double> > distance_train_test(N_train, vector<double>(N_test,0));
    for (int i=0; i < N_train; i++){
        for (int j=0; j < N_test; j++){
            distance_train_test[i][j] = distance2d(dist_train[i], sim_train[i], dist_test[j], sim_test[j]);
        }
    }

    vector<vector<double> > distance(N_train+N_test, vector<double>(N_train+N_test,0));

    for (int i=0; i < N_train; i++){
        for (int j=0; j < i; j++){
            distance[i][j] = distance_train_train[i][j];
        }
    }
    for (int i=0; i < +N_test; i++){
        for (int j=0; j < i; j++){
            distance[i+N_train][j+N_train] = distance_test_test[i][j];
        }
    }
    for (int i=0; i < N_train; i++){
        for (int j=0; j < N_test; j++){
            distance[i][j+N_test] = distance_train_test[i][j];
        }
    }
    return distance;
}