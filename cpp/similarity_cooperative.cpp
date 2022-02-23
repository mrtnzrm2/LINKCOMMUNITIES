#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <iomanip>

using namespace std;

vector<vector<long double> > read_csv(string &fname, bool &header){
	vector<vector<long double> > content;
	vector<long double> row;
	string line, word;

    fstream file (fname, ios::in);
    
    if (!header){
        int start = 0;
        if(file.is_open()){
            while(getline(file, line)){
                if (start > 0){
                    row.clear();
    
                    stringstream str(line);
        
                    while(getline(str, word, ',')){
                        row.push_back(stold(word));
                    }  
                    content.push_back(row);
                }
                start++;
            }
        }
        else
            cout<<"Could not open the file\n";       
    } else{
        if(file.is_open()){
            while(getline(file, line)){
                row.clear();
    
                stringstream str(line);

                while(getline(str, word, ',')){
                        row.push_back(stold(word));
                    }  
                content.push_back(row);
            }
        }
        else
            cout<<"Could not open the file\n";
    }

    int flag=0;
    for (int i=0; i < content.size(); i++){
        if (content.size() != content[i].size())
            flag = 1;
    }
    if (flag == 1)
        cout << "Warning: asymmetric csv\n";

	return content;
}

long double mean_nonzero(vector<long double> &v){
	int non_zero = 0;
	long double mean_v = 0;

	for (int i=0; i < v.size(); i++){
		mean_v += v[i];
		if (v[i] != 0){
			non_zero++;
		}
	}
	
	if (non_zero != 0)
		mean_v /= non_zero;
	else{
		cout << "Vector with only zero entries" << endl;
	}

	return mean_v;
}

vector<vector<long double> > create_aik(vector<vector<long double> > &data, const string &average, const int N){

	vector<vector<long double> > aik(N, vector<long double>(N,0));
	vector<long double> vec(N);

	aik = data;

	for (int i=0; i < N; i++){
		if (average.compare("ALPHA")){
			aik[i][i] = mean_nonzero(aik[i]);
		} else if (average.compare("BETA")){
			for (int j=0; j < N; j++)
				vec[j] = aik[j][i];
			aik[i][i] = mean_nonzero(vec);
		} else{
			cout << "Average is not ALPHA or BETA\n";
		}
	}

	return aik;
}

vector<vector<long double> > create_aki(vector<vector<long double> > &data, const string &average, const int N){

	vector<vector<long double> > aki(N, vector<long double>(N,0));

	for (int i=0; i < N; i++){
		for (int j=0; j < N; j++){
			aki[i][j]=data[j][i];
		}
	}

	for (int i=0; i < N; i++){
		if (average.compare("ALPHA")){
			aki[i][i] = mean_nonzero(aki[i]);
		} else if (average.compare("BETA")){
			aki[i][i] = mean_nonzero(data[i]);
		} else{
			cout << "Average is not ALPHA or BETA\n";
		}
	}

	return aki;
}

vector<vector<int_fast32_t> > create_id_matrix(vector<vector<long double> > &matrix, const int &N){
	
	vector<vector<int_fast32_t> > id_matrix(N, vector<int_fast32_t>(N,-1));
	int id = 0;

	for (int i=0; i < N; i++){
		for (int j=0; j < N; j++){
			if (matrix[i][j] > 0){
				id_matrix[i][j] = id;
				id++;
			}
		}
	}

	return id_matrix;
}

long double jaccard_p(vector<long double> &u, vector<long double> &v, const int &N){
	long double JACP = 0;
	long double p;

	for (int i=0; i < N; i++){
		if (u[i] != 0 && v[i] !=0){
			p = 0;
			for (int j=0; j < N; j++){
				p += max(u[j]/u[i], v[j]/v[i]);
			}
			if (p != 0)
				JACP += 1/p;
			else
				cout << "Vectors in jaccardp are both zero";
		}
	}

	return JACP;
}

vector<vector<long double> > up_triangular_rows(vector<vector<long double> > &matrix, const int &N){
    vector<vector<long double> > dist_h;
    vector<long double> row;

    for (int i=0; i < N - 1; i++){
        row.clear();
        for (int j=0; j < N-i-1; j++){
            row.push_back(matrix[i][j+i+1]);
        }
        dist_h.push_back(row);
    }

    return dist_h;
}

vector<vector<long double> > calculate_sim_matrix(vector<vector<long double> > &matrix, vector<vector<long double> > &aik, vector<vector<long double> > &aki, const int &N, const int &leaves){

	vector<vector<int_fast32_t> > id_matrix = create_id_matrix(matrix, N);
	vector<vector<long double> > leave_matrix(leaves, vector<long double>(leaves,-1));
	int_fast32_t col_id;
	int_fast32_t row_id;

	for (int i =0; i < N; i++){
		for (int j=0; j < N; j++){
			row_id = id_matrix[i][j];
			if (row_id == -1)
				continue;
			for (int k=0; k < N; k++){
				col_id = id_matrix[i][k];
				if (k != j && col_id != -1){
					if (row_id < col_id){
						if (leave_matrix[row_id][col_id] == -1)
							leave_matrix[row_id][col_id] = jaccard_p(aki[j], aki[k], N);	
					} else{
						if (leave_matrix[col_id][row_id] == -1)
							leave_matrix[col_id][row_id] = jaccard_p(aki[j], aki[k], N);
					}
				}
				col_id = id_matrix[k][j];
				if (k != i && col_id != -1){
					if (row_id < col_id){
						if (leave_matrix[row_id][col_id] == -1)
							leave_matrix[row_id][col_id] = jaccard_p(aik[i], aik[k], N);
					} else{
						if (leave_matrix[col_id][row_id] == -1)
							leave_matrix[col_id][row_id] = jaccard_p(aik[i], aik[k], N);
					}
				}
			}
		}
	}
	
	return leave_matrix;
}

void out_csv(vector<vector<long double> > &matrix, string &out_path, const int leaves){
	ofstream out(out_path);
	for (int i=0; i < leaves; i++){
		for (int j=0; j < leaves; j++){
			if (j < leaves -1)
				out << fixed << setprecision(10) << matrix[i][j] << ',';
			else
				out << fixed << setprecision(10) << matrix[i][j];
			// if (j < leaves -1)
			// 	out << matrix[i][j] << ',';
			// else
			// 	out << matrix[i][j];
		}
		out << '\n';
	}
}

void out_csv_up_triangular(vector<vector<long double> > &matrix, string &out_path, const int leaves){
	ofstream out(out_path);
	for (int i=0; i < leaves-1; i++){
		for (int j=0; j < leaves-i-1; j++){
			if (j < leaves-i-2)
				out << fixed << setprecision(10) << matrix[i][j] << ',';   // 10
			else
				out << fixed << setprecision(10) << matrix[i][j];  // 10
			// if (j < leaves-i-2)
			// 	out <<  matrix[i][j] << ',';   // 10
			// else
			// 	out << matrix[i][j];  // 10
		}
		out << '\n';
	}
}

template <typename T, typename I>
void nlog10_matrix(vector<vector<T> > &matrix, const I &N){

    for (int i=0; i < N; i++){
        for (int j=0; j < N; j++){
            if (matrix[i][j] > 0)
                matrix[i][j] = -log10(matrix[i][j]);
        }
    }
}


int main(){

	bool header;
	string out_path;
	string csv_path;

	const string average = "ALPHA";
	const int N = 107;
	const int leaves = 6868; // 6880

	header = false;
	csv_path = "../CSV/merged/imputation/tracto2016/zz_model/fln_4_r_6_3.csv"; // ../CSV/merged/imputation/tracto2016/zz_model/NONULL/
	cout << "Loading fln\n";
    vector<vector<long double> > fln_matrix = read_csv(csv_path, header);
	cout << "Finished\n";

	cout << "Formatting data\n";
	nlog10_matrix(fln_matrix, N);
	cout << "Finished\n";

	cout << "Computing AIK and AKI\n";
	vector<vector<long double> > aik = create_aik(fln_matrix, average, N);
	vector<vector<long double> > aki = create_aki(fln_matrix, average, N);
	cout << "Finished\n";

	cout << "Compute similarity matrix\n";
	vector<vector<long double> > sim_matrix = calculate_sim_matrix(fln_matrix, aik, aki, N, leaves);
	cout << "Finished\n";

	cout << "Saving\n";
	out_path = "../CSV/merged/similarity/tracto2016/zz_model/sim_full_l10_4_r_6_3.csv";
	out_csv(sim_matrix, out_path, leaves);
	sim_matrix = up_triangular_rows(sim_matrix, leaves);
	out_path = "../CSV/merged/similarity/tracto2016/zz_model/sim_l10_4_r_6_3.csv";
	out_csv_up_triangular(sim_matrix, out_path, leaves);
	cout << "End\n";

    return 0;
}