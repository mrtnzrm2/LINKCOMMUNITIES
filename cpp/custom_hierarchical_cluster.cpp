#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <chrono>
#include <math.h>
#include <stdio.h>
#include <limits.h>
using namespace std::chrono;
using namespace std;

vector<vector<long double> > read_csv(const string fname){
	vector<vector<long double> > content;
	vector<long double> row;
	string line, word;
 
	fstream file (fname, ios::in);
	if(file.is_open())
	{
		while(getline(file, line))
		{
			row.clear();
 
            stringstream str(line);

            while(getline(str, word, ','))
                row.push_back(stod(word));
            content.push_back(row);
		}
	}
	else
		cout<<"Could not open the file\n";

	return content;
}

template <typename T, typename I>
void custom_dis_matrix(vector<vector<T> > &matrix, const I &N){
    
    T max_dis = -1;

    for (int i=0; i < N-1; i++){
        for (int j=0; j < N-i-1; j++){
            if (matrix[i][j] > 0){
                matrix[i][j] = -log10(matrix[i][j]);
                if (max_dis < matrix[i][j])
                    max_dis = matrix[i][j];
            }
        }
    }

    for (int i=0; i < N-1; i++){
        for (int j=0; j < N-i-1; j++){
            if (matrix[i][j] <= 0)
                matrix[i][j] = max_dis + 1;
        }
    }
}

template <typename T, typename I>
vector<vector<T> > up_triangular_rows(vector<vector<T> > &matrix, I &N){
    vector<vector<T> > dist_h;
    vector<T> row;

    for (int i=0; i < N - 1; i++){
        row.clear();
        for (int j=0; j < N-i-1; j++){
            row.push_back(matrix[i][j+i+1]);
        }
        dist_h.push_back(row);
    }

    return dist_h;
}

template <typename T, typename I>
struct hclust{
    vector<vector<I> > merge;
    vector<T> height; 
};

template <typename I>
struct nodes{
    I root;
    vector<I> members;
};

template <typename I>
struct saver{
    I u;
    I v;
};

template <typename T, typename I>
T single_linkage(vector<I> &m1, vector<I> &m2, vector<vector<T> >  &dist){

    T min_dist=INT_MAX;
    T dist1;
    T dist2;
    I N2 = m2.size();

    for (int i=0; i < m1.size(); i++){
        for (int j=0; j <= N2/2; j++){

            if (m1[i] > m2[j] && m1[i] > m2[N2-j-1]){
                dist1 = dist[m2[j]][m1[i]-m2[j]-1];
                dist2 = dist[m2[N2-j-1]][m1[i]-m2[N2-j-1]-1];
                if (dist1 > dist2){
                    if (min_dist > dist2)
                        min_dist = dist2;
                } else{
                    if (min_dist > dist1)
                        min_dist = dist1;
                }
            } else if (m1[i] < m2[j] && m1[i] < m2[N2-j-1]){
                dist1 = dist[m1[i]][m2[j]-m1[i]-1];
                dist2 = dist[m1[i]][m2[N2-j-1]-m1[i]-1];
                if (dist1 > dist2){
                    if (min_dist > dist2)
                        min_dist = dist2;
                } else{
                    if (min_dist > dist1)
                        min_dist = dist1;
                }
            } else if (m1[i] < m2[j] && m1[i] > m2[N2-j-1]){
                dist1 = dist[m1[i]][m2[j]-m1[i]-1];
                dist2 = dist[m2[N2-j-1]][m1[i]-m2[N2-j-1]-1];
                if (dist1 > dist2){
                    if (min_dist > dist2)
                        min_dist = dist2;
                } else{
                    if (min_dist > dist1)
                        min_dist = dist1;
                }
            } else if (m1[i] > m2[j] && m1[i] < m2[N2-j-1]){
                dist1 = dist[m2[j]][m1[i]-m2[j]-1];
                dist2 = dist[m1[i]][m2[N2-j-1]-m1[i]-1];
                if (dist1 > dist2){
                    if (min_dist > dist2)
                        min_dist = dist2;
                } else{
                    if (min_dist > dist1)
                        min_dist = dist1;
                }
            } else
                cout << "The same node is in two different clusters\n";
        }
    }
    return min_dist;
}

long double average_linkage(vector<int_fast32_t> &m1, vector<int_fast32_t> &m2, vector<vector<long double> >  &dist){

    long double dist1 = 0;
    for (int i=0; i < m1.size(); i++){
        for (int j=0; j < m2.size(); j++){

            if (m1[i] == m2[j]){
                cout << "The same node is in two different clusters\n";
            } else if (m1[i] > m2[j]){
                dist1 += dist[m2[j]][m1[i]-m2[j]-1];
            } else{
                dist1 += dist[m1[i]][m2[j]-m1[i]-1];
            }
        }
    }

    return dist1/(m1.size()*m2.size());
}

template<typename I>
vector<vector<int_fast32_t> > nodes_parallelization(I &N, const I &chunks){

    vector<vector<int_fast32_t> > content(2, vector<I>(chunks, 0));

    I sep = (int) N/chunks;

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

template <typename T, typename I>
hclust<T,I> hierarchical_clustering_parallel(vector<vector<T> >  &dist, const I &steps, const I &th, const I &chunks){
    hclust<T,I> H;
    H.merge = vector<vector< int_fast32_t> > (steps-1,vector<int_fast32_t>(2,0));
    H.height = vector<long double> (steps-1,0);
    vector<long double> dist_vector;
    vector<saver<I> > id_saver(chunks);

    vector<int_fast32_t> inter;
    vector<vector<int_fast32_t> > range;

    long double min_dist;
    long double dist_g1;
    long double dist_g2;

    int_fast32_t N;
    int_fast32_t u;
    int_fast32_t v;

    int pi, x, y;

    vector<nodes<I> > hierarchy(steps);
    for (int i=0; i < steps; i++){
        hierarchy[i].root = -(i+1);
        hierarchy[i].members.push_back(i);
    }

    for (int e=0; e < steps-1; e++){

        // if ((int) 100*e/(steps-1)%10 == 0 && (int) 100*e/(steps-1) >= 1)
        //     cout << (int) 100*e/(steps-1)<< '%' << endl;

        min_dist = INT_MAX;
        N = hierarchy.size();
        if (N <= 1)
            continue;
        if (N > chunks){
            dist_vector.assign(chunks, INT_MAX);
            range = nodes_parallelization(N, chunks);
            #pragma omp parallel for default(none)  private(pi, x, y) shared(hierarchy, dist, chunks, th, e, dist_g1, dist_g2, N, dist_vector, id_saver, range)
            for (pi=0; pi < chunks; pi++){
                // cout << omp_get_num_threads() << endl;
                for (x=range[0][pi]; x < range[1][pi]; x++){
                    for (y=0; y <= (N-x-1)/2; y++){
                        if (e < th){
                            dist_g1 = single_linkage(hierarchy[x].members, hierarchy[y+x+1].members, dist);
                            dist_g2 = single_linkage(hierarchy[x].members, hierarchy[N-y-1].members, dist);
                        } else{
                            dist_g1 = average_linkage(hierarchy[x].members, hierarchy[y+x+1].members, dist);
                            dist_g2 = average_linkage(hierarchy[x].members, hierarchy[N-y-1].members, dist);
                        }
                        if (dist_g1 < dist_g2){
                            if(dist_vector[pi] > dist_g1){
                                dist_vector[pi] = dist_g1;
                                id_saver[pi].u = x;
                                id_saver[pi].v = y+x+1;
                            }
                        } else{
                            if(dist_vector[pi] > dist_g2){
                                dist_vector[pi] = dist_g2;
                                id_saver[pi].u = x;
                                id_saver[pi].v = N-y-1;
                            }
                        }
                    }
                }
            }
            for (int pj=0; pj <= chunks/2; pj++){
                if (dist_vector[pj] > dist_vector[chunks-pj-1]){
                    if (min_dist > dist_vector[chunks-pj-1]){
                         min_dist = dist_vector[chunks-pj-1];
                         u = id_saver[chunks-pj-1].u;
                         v = id_saver[chunks-pj-1].v;
                    }
                } else{
                    if (min_dist > dist_vector[pj]){
                        min_dist = dist_vector[pj];
                        u = id_saver[pj].u;
                        v = id_saver[pj].v;
                    }  
                }
            }
        } else{
            for (int x=0; x < N-1; x++){
                for (int y=0; y <= (N-x-1)/2; y++){
                    if (e < th){
                        dist_g1 = single_linkage(hierarchy[x].members, hierarchy[y+x+1].members, dist);
                        dist_g2 = single_linkage(hierarchy[x].members, hierarchy[N-y-1].members, dist);
                    } else{
                        dist_g1 = average_linkage(hierarchy[x].members, hierarchy[y+x+1].members, dist);
                        dist_g2 = average_linkage(hierarchy[x].members, hierarchy[N-y-1].members, dist);
                    }
                    if (dist_g1 < dist_g2){
                        if(min_dist > dist_g1){
                            min_dist = dist_g1;
                            u = x;
                            v = y+x+1;
                        }
                    } else{
                        if(min_dist > dist_g2){
                            min_dist = dist_g2;
                            u = x;
                            v = N-y-1;
                        }
                    }
                }
            }
        }

        H.height[e] = min_dist/2;
        H.merge[e][0] = hierarchy[u].root;
        H.merge[e][1] = hierarchy[v].root;
        hierarchy[u].root = e+1;
        inter.clear();
        inter.reserve(hierarchy[u].members.size() + hierarchy[v].members.size());
        inter.insert(inter.end(), hierarchy[u].members.begin(), hierarchy[u].members.end());
        inter.insert(inter.end(), hierarchy[v].members.begin(), hierarchy[v].members.end());
        hierarchy[u].members = inter;
        hierarchy.erase(hierarchy.begin()+v);
    }

    return H;
}

template<typename T, typename I>
hclust<T,I> hierarchical_clustering(vector<vector<T> >  &dist, const I &steps, const I &th){
    hclust<T,I> H;
    H.merge = vector<vector<I> > (steps-1,vector<I>(2,0));
    H.height = vector<T> (steps-1,0);
    vector<I> inter;
    T min_dist;
    T dist_g1;
    T dist_g2;
    I N;
    I u;
    I v;

    vector<nodes<I> > hierarchy(steps);
    for (int i=0; i < steps; i++){
        hierarchy[i].root = -(i+1);
        hierarchy[i].members.push_back(i);
    }

    for (int e=0; e < steps-1; e++){

        if ((int) 100*e/(steps-1)%10 == 0 && (int) 100*e/(steps-1) >= 1)
            cout << (int) 100*e/(steps-1)<< '%' << endl;

        min_dist = INT_MAX;
        N = hierarchy.size();
        if (N <= 1)
            continue;
        for (int x=0; x < N-1; x++){
            for (int y=0; y <= (N-x-1)/2; y++){
                if (e < th){
                    dist_g1 = single_linkage(hierarchy[x].members, hierarchy[y+x+1].members, dist);
                    dist_g2 = single_linkage(hierarchy[x].members, hierarchy[N-y-1].members, dist);
                } else{
                    dist_g1 = average_linkage(hierarchy[x].members, hierarchy[y+x+1].members, dist);
                    dist_g2 = average_linkage(hierarchy[x].members, hierarchy[N-y-1].members, dist);
                }
                if (dist_g1 < dist_g2){
                    if(min_dist > dist_g1){
                        min_dist = dist_g1;
                        u = x;
                        v = y+x+1;
                    }
                } else{
                    if(min_dist > dist_g2){
                        min_dist = dist_g2;
                        u = x;
                        v = N-y-1;
                    }
                }
            }
        }

        H.height[e] = min_dist/2;
        H.merge[e][0] = hierarchy[u].root;
        H.merge[e][1] = hierarchy[v].root;
        hierarchy[u].root = e+1;
        inter.clear();
        inter.reserve(hierarchy[u].members.size() + hierarchy[v].members.size());
        inter.insert(inter.end(), hierarchy[u].members.begin(), hierarchy[u].members.end());
        inter.insert(inter.end(), hierarchy[v].members.begin(), hierarchy[v].members.end());
        hierarchy[u].members = inter;
        hierarchy.erase(hierarchy.begin()+v);
    }

    return H;
}

template <typename S>
ostream& operator<<(ostream &os, const vector<S> &vector){
    // Printing all the elements
    // using <<
    for (auto element : vector) {
        os << element << " ";
    }
    return os;
}

template <typename I>
void out_csv_2d(vector<vector<I> > &matrix, const string &out_path){

    I sx = matrix.size();
    I sy = matrix[0].size();
	ofstream out(out_path);
	for (int i=0; i < sx; i++){
		for (int j=0; j < sy; j++){
			if (j < sy -1)
				out << matrix[i][j] << ',';
			else
				out << matrix[i][j];
		}
		out << '\n';
	}
}

template <typename T>
void out_csv_1d(vector<T> &matrix, const string &out_path){
	ofstream out(out_path);
	for (int i=0; i < matrix.size(); i++){
		out << matrix[i] << '\n';
	}
}

int main(){

    string out_path;
    const string csv_path = "../CSV/91x40/sim.csv"; // merged/similarity/tracto2016/zz_model/
	const int_fast32_t leaves = 999; // 6868
    const int_fast32_t threshold = leaves - 181; //
    const int_fast32_t chunks = 6;
    vector<vector<long double> > distance = read_csv(csv_path);

    custom_dis_matrix(distance, leaves);

    auto start = high_resolution_clock::now();
    hclust<long double, int_fast32_t> hyr = hierarchical_clustering(distance, leaves, threshold);
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    cout << duration.count() << endl;
    
    vector<vector<int_fast32_t> > merge = hyr.merge;
    vector<long double> height = hyr.height;

    out_path = "../CSV/91x40/merge_th_181.csv";
    out_csv_2d(merge, out_path);
    out_path = "../CSV/91x40/height_th_181.csv";
    out_csv_1d(height, out_path);

    return 0;
}