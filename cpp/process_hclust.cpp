#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "hclust-cpp/fastcluster.h"

// [[Rcpp::plugins(cpp14)]]

void vec_times_k(std::vector<int>& v, int k) {
    std::transform(v.begin(), v.end(), v.begin(), [k](int &c){ return c*k; });
}

double get_pss(std::vector<int>& lc, int& size, int& M, int& n, double& th) {
  vec_times_k(lc, -1);
  std::sort(lc.begin(), lc.end());
  vec_times_k(lc, -1);
  double pss = 0;
  double D;
  if (size >= n) {
    for (int i=0; i < n - 1; i++) {
      for (int j=(i + 1); j < n; j++) {
        D = static_cast<double>(lc[j]) / static_cast<double>(lc[i]);
        if (D > th)
          pss += D * (static_cast<double>(lc[i] + lc[j]) / static_cast<double>(M));
        else
          pss -= D * (static_cast<double>(lc[i] + lc[j]) / static_cast<double>(M));
      }
    }
    pss /= 0.5 * n * (n - 1);
  } else {
    for (int i=0; i < size - 1; i++) {
      for (int j=(i + 1); j < size; j++) {
        D = static_cast<double>(lc[j]) / static_cast<double>(lc[i]);
        if (D > th)
          pss += D * (static_cast<double>(lc[i] + lc[j]) / static_cast<double>(M));
        else
          pss -= D * (static_cast<double>(lc[i] + lc[j]) / static_cast<double>(M));
      }
    }
    pss /= 0.5 * size * (size - 1);
  }
  return pss;
}

double approx_nodes_by_edges(int& M) {
  double x;
  double m = static_cast<double>(M);
  double x1 = (1 + sqrt(1 + 4 * m)) / 2;
  x1 = floor(x1);
  double x2 = (1 - sqrt(1 + 4 * m)) / 2;
  x2 = floor(x2);
  if (x2 <= 0)
    x = x1;
  else
    x = x2;
  return x;
}

void get_sizes(
  int* labels,
  std::vector<int>& lcs,
  std::vector<int>& nds,
  std::vector<int>& ulabels,
  std::vector<int>& source,
  std::vector<int>& target,
  int& n
) {

  std::vector<std::vector <int>> node_buffer(ulabels.size());
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < ulabels.size(); j++) {
      if (labels[i] == ulabels[j]) {
        lcs[j]++;
        node_buffer[j].push_back(source[i]);
        node_buffer[j].push_back(target[i]);
      }
    }
  }
  for (int j = 0; j < ulabels.size(); j++) {
    std::vector<int>::iterator ip;
    std::sort(node_buffer[j].begin(), node_buffer[j].end());
    ip = std::unique(
      node_buffer[j].begin(),
      node_buffer[j].begin() + node_buffer[j].size()
    );
    node_buffer[j].resize(
      std::distance(node_buffer[j].begin(), ip)
    );
    nds[j] = node_buffer[j].size();
  }
}

std::vector<double> simplify_height_to_k_end(
  int &n,
  double* height,
  std::vector<double>& sim_height,
  int* size
) {
  double h = height[0];
  std::vector<double> sim_k;
  for (int i = 0; i < n - 1; i++) {
    if (i < n - 2) {
      if (height[i + 1] != h) {
        sim_k.push_back(n - i);
        sim_height.push_back(h);
        h = height[i + 1];
        ++(*size);
      }
    } else {
      if (height[i] != height[i - 1]) {
        sim_k.push_back(n - i);
        h = height[i];
        sim_height.push_back(h);
        ++(*size);
      }
    }
    
  }
  return sim_k;
}

std::vector<double> simplify_height_to_k_start(
  int &n,
  double* height,
  std::vector<double>& sim_height,
  int* size
) {
  double h = height[0];
  std::vector<double> sim_k;
  for (int i = 0; i < n - 1; i++) {
    if (i == 0) {
      sim_k.push_back(n);
      sim_height.push_back(h);
      ++(*size);
    }
    if (height[i] != h && i != 0) {
      h = height[i];
      sim_k.push_back(n - i);
      sim_height.push_back(h);
      ++(*size);
    }
  }
  return sim_k;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix process_hclust_fast(
  int& n,
  std::vector<std::vector <double>>& distmat,
  std::vector<int>& source,
  std::vector<int>& target,
  int& nodes,
  int& nss,
  double& th
)
{
  // From distance matrix to upper triangular array
  double* tri_distmat = new double[(n*(n-1))/2];
  for (int i = 0, k = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      tri_distmat[k] = distmat[i][j];
      k++;
    }
  }
  // hclust
  int* merge = new int[2*(n-1)];
  double* height = new double[n-1];
  hclust_fast(
    n,
    tri_distmat,
    2, // Average method
    merge,
    height
  );
  // simplify height to k
  std::vector<double> k;
  std::vector<double> sim_height;
  int* K = new int;
  *K = 0;
  k = simplify_height_to_k_start(n, height, sim_height, K);
  std::cout << "Number of simplified heights: " << *K << std::endl;
  Rcpp::NumericMatrix features_main(4, *K);
  for (int i = 0; i < *K; i++)
    features_main(3, i) = sim_height[i];
  // 0: K
  // 1: Dc
  // 2: SS
  // 3: Height
  // THE GAME STARTS
  for (int i=0; i < *K; i++) {
    int* labels = new int[n];
    cutree_k(
      n,
      merge,
      k[i],
      labels
    );
    std::vector<int> ulabels(labels, labels + n);
    std::vector<int>::iterator ip;
    std::sort(ulabels.begin(), ulabels.end());
    ip = std::unique(
      ulabels.begin(),
      ulabels.begin() + n
    );
    ulabels.resize(
      std::distance(ulabels.begin(), ip)
    );
    int number_lcs = ulabels.size();
    // K
    features_main(0, i) = k[i];
    // Get LCs sizes in order
    std::vector<int> lc_size(number_lcs, 0);
    std::vector<int> node_size(number_lcs, 0);
    get_sizes(
      labels, lc_size,
      node_size, ulabels,
      source, target,
      n
    );
    double Dc;
    for (int j=0; j < number_lcs; j++) {
      if (lc_size[j] > 1 && node_size[j] > 2) {
        Dc = static_cast<double>(lc_size[j] -  node_size[j] + 1) /
            static_cast<double>((node_size[j] - 1) * (node_size[j] - 1));
        features_main(1, i) +=  (static_cast<double>(lc_size[j]) / static_cast<double>(n)) * Dc;
      }
    }
    features_main(2, i) = get_pss(
      lc_size, number_lcs,
      n, nss, // ss average range
      th
    );
    delete[] labels;
  }
  // Delete phase
  delete K;
  delete[] merge;
  delete[] height;
  delete[] tri_distmat;
  // Return
  return features_main;
}

