#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "hclust-cpp/fastcluster.h"

// [[Rcpp::plugins(cpp14)]]

std::vector<int> link_memberships_area(
  int& area,
  int& leaves,
  int* label,
  std::vector<int>& source,
  std::vector<int>& target
) {
  std::vector<int> lms;
  for (int i = 0; i < leaves; i++) {
    if (source[i] == area || target[i] == area)
      lms.push_back(label[i]);
  }
  return lms;
}

std::vector<int> hist_memberships_area(
  std::vector<int>& lms
) {
  int lms_size = lms.size();
  std::vector<int> c_lms = lms;
  // Find unique LC memberships
  std::vector<int>::iterator ip;
  std::sort(lms.begin(), lms.end());
  ip = std::unique(
    lms.begin(),
    lms.begin() + lms_size
  );
  lms.resize(
    std::distance(lms.begin(), ip)
  );
  std::vector<int> hist(lms.size(), 0);
  for (int i = 0; i < lms_size; i++) {
    for (int j = 0; j < lms.size(); j++) {
      if (c_lms[i] == lms[j])
        hist[j]++;
    }
  }
  return hist;
}

bool if_area_in_chain(
  Rcpp::NumericMatrix& chain,
  int& area,
  int& nodes
) {
  bool is_area = false;
  for (int i = 0; i < nodes; i++) {
    if (chain(i, 0) == area)
      is_area = true;
  }
  return is_area;
}

bool area_check(
  int& area,
  int& K,
  std::vector<int>& hist,
  std::vector<int>& indeg,
  std::vector<int>& outdeg,
  double p
) {
  bool state = false;
  int hist_size = hist.size();
  for (int i = 0; i < hist_size; i++) {
    if (hist[i] >= p * indeg[area] || hist[i] >= p * outdeg[area])
      state = true;
  }
  return state;
}

void aggregate_node(
  int& area,
  int& area2,
  int& inode,
  int& K,
  int& nodes,
  Rcpp::NumericMatrix& chain,
  std::vector<int>& alinks,
  std::vector<int>& ajlinks,
  std::vector<int>& ahist,
  std::vector<int>& ajhist,
  std::vector<int>& in_deg,
  std::vector<int>& out_deg,
  std::vector<std::string>& labels
) {
  for (int i = 0; i < alinks.size(); i++) {
    for (int j = 0; j < ajlinks.size(); j++) {
      if (alinks[i] == ajlinks[j]) {
        if (
          (ahist[i]/in_deg[area - 1] > 0.4 &&
          ajhist[j]/in_deg[area2 - 1] > 0.4) ||
          (ahist[i]/out_deg[area - 1] > 0.4 &&
          ajhist[j]/out_deg[area2 - 1] > 0.4) 
        ) {
          if (!if_area_in_chain(chain, area2, nodes)) {
            std::cout << "Step: " << K << " ladder: " << inode << " area name: " << labels[area2 - 1] << std::endl;
            chain(inode, 0) = area2;
            chain(inode, 1) = K;
            inode++;
          }
        }
      }
    }
  }
}

template<typename T>
std::vector<T> sum_v_a_b(std::vector<T>& a, std::vector<T>& b, int& nodes) {
  std::vector<T> sum_v(nodes, 0);
  for (int i = 0; i < nodes; i++) {
    sum_v[i] = a[i] + b[i];
  }
  return sum_v;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix area_hierarchy(
  int& area,
  int& leaves,
  int& nodes,
  std::vector<int>& source,
  std::vector<int>& target,
  std::vector<int>& indeg,
  std::vector<int>& outdeg,
  std::vector<std::string>& labs,
  std::vector<std::vector <double>>& distmat
  ) {
  // Sum up degrees
  std::vector<int> sum_deg = sum_v_a_b(indeg, outdeg, nodes);
  // From distance matrix to upper triangular array
  double* tri_distmat = new double[(leaves * (leaves - 1)) / 2];
  for (int i = 0, k = 0; i < leaves; i++) {
    for (int j = i+1; j < leaves; j++) {
      tri_distmat[k] = distmat[i][j];
      k++;
    }
  }
  // hclust
  int* merge = new int[2 * (leaves - 1)];
  double* height = new double[leaves - 1];
  hclust_fast(
    leaves,
    tri_distmat,
    2, // Average method
    merge,
    height
  );
  // Define area hierarchy
  Rcpp::NumericMatrix area_h(nodes, 2);
  for (int i = 0; i < nodes; i++) {
    area_h(i, 0) = -1;
  }
  int K;
  // Start building up the hierarchy
  for (int i = 1, i_nodes = 0; i < leaves - 1; i++) {
    K = leaves - i;
    int* labels = new int[leaves];
    // Cut tree
    cutree_k(
      leaves,
      merge,
      K,
      labels
    );
    // Assign link membership to the area's links
    std::vector<int> area_links = link_memberships_area(
      area,
      leaves,
      labels,
      source,
      target
    );
    // Compute a histogram of the area's links memberships
    std::vector<int> area_hist = hist_memberships_area(area_links);
    // Check if most of the area's LC belong to a single LC
    if (
      area_check(
        area,
        K,
        area_hist,
        indeg,
        outdeg,
        0.6
      )
    ) {
      if (!if_area_in_chain(area_h, area, nodes)) {
        std::cout << "Step: " << K << " ladder: " << i_nodes << " area name: " << labs[area - 1] << std::endl;
        area_h(i_nodes, 0) = area;
        area_h(i_nodes, 1) = leaves - i;
        i_nodes++;
      }
      for (int j = 0; j < nodes; j++) {
        if (area == j)
          continue;
        std::vector<int> area_j_links = link_memberships_area(
          j,
          leaves,
          labels,
          source,
          target
        );
        std::vector<int> area_j_hist = hist_memberships_area(area_j_links);
        aggregate_node(
          area,
          j,
          i_nodes,
          K,
          nodes,
          area_h,
          area_links,
          area_j_links,
          area_hist,
          area_j_hist,
          indeg,
          outdeg,
          labs
        );
      }
    }
    else {
      
    }
    delete[] labels; 
  }
  return area_h;
}