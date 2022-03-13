#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <numeric>
#include <cstdio>
#include <utility>
#include <stdio.h>
#include <set>

using namespace std;

vector<vector<long double> > read_csv(char fname[], bool &header){
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

vector<vector<int_fast16_t> > read_icsv(char fname[], bool &header){
	vector<vector<int_fast16_t> > content;
	vector<int_fast16_t> row;
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
                        row.push_back(stoi(word));
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
                        row.push_back(stoi(word));
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

template <typename gfloat, typename gint>
vector<vector<gint> > create_id_matrix(vector<vector<gfloat> > &matrix, gint &N){
	
	vector<vector<gint> > id_matrix(N, vector<gint>(N,-1));
	gint id = 0;

	for (gint i=0; i < N; i++){
		for (gint j=0; j < N; j++){
			if (matrix[i][j] > 0){
				id_matrix[i][j] = id;
				id++;
			}
		}
	}

	return id_matrix;
}

template <typename gfloat>
gfloat average(vector<gfloat> const& v){
    if(v.empty()){
        return 0;
    }

    auto const count = static_cast<gfloat>(v.size());
    return std::reduce(v.begin(), v.end()) / count;
}

template <typename gfloat>
double stdev(vector<gfloat> v, gfloat ave)
{
    gfloat E=0;
    // Quick Question - Can vector::size() return 0?
    gfloat inverse = 1.0 / static_cast<gfloat>(v.size());
    for(unsigned int i=0;i<v.size();i++)
    {
        E += pow(static_cast<gfloat>(v[i]) - ave, 2);
    }
    return sqrt(inverse * E);
}

template <typename gfloat>
struct summary{
    gfloat mean;
    gfloat stdn;
};

template <typename gfloat>
void nlog10_matrix(vector<vector<gfloat> > &matrix){
    
    for (int_fast16_t i=0; i < matrix.size(); i++){
        for (int_fast16_t j=0; j < matrix[i].size(); j++){
            if (matrix[i][j] > 0)
                matrix[i][j] = -log10(matrix[i][j]);
        }
    }
}

template <typename gfloat>
gfloat find_max_special(vector<vector<gfloat> > &matrix){
     gfloat maxe=0;

    for (int_fast16_t i=0; i < matrix.size(); i++){
        for (int_fast16_t j=0; j < matrix[i].size(); j++){
            if (matrix[i][j] > maxe)
                maxe = matrix[i][j];
        }
    }
    return maxe;
}

template <typename gfloat>
void replace_max(vector<vector<gfloat> > &matrix){

    gfloat maxe = find_max_special(matrix);

    for (int_fast16_t i=0; i < matrix.size(); i++){
        for (int_fast16_t j=0; j < matrix[i].size(); j++){
            if (matrix[i][j] == -1)
                matrix[i][j] = maxe + 1;
        }
    }
}

// template <typename gfloat, typename gint>
// vector<gint> RnB_neigborhood(vector<vector<gint> > &id, vector<vector<gfloat> > &red, vector<vector<gfloat> > &blue, gint &N){

//     set<gint> s;
//     vector<gint> neigborhood;
//     gint id_red, id_blue;

//     for (gint i=0; i < N; i++){
//         for (gint j=0; j < N; j++){
//             id_red = id[i][j];
//             if (id_red == -1 || red[i][j] == 0)
//                 continue;
//             for (gint k=0; k < N; k++){
//                 id_blue = id[i][k];
//                 if (blue[i][k] == 1 && id_blue != -1 && id_blue != id_red)
//                     neigborhood.push_back(id_blue);

//                 id_blue = id[k][j];
//                 if (blue[k][j] == 1 && id_blue != -1 && id_blue != id_red)
//                     neigborhood.push_back(id_blue);
//             }
//         }
//     }

//     for (gint i=0; i < neigborhood.size(); i++)
//         s.insert(neigborhood[i]);
//     neigborhood.assign(s.begin(), s.end());

//     return neigborhood;
// }

template <typename gfloat, typename gint>
vector<gint> link_community_ids(vector<vector<gint> > &id, vector<vector<gfloat> > &color, gint &N ){
    set<gint> s;
    vector<gint> ids;

    for (gint i=0; i < N; i++){
        for (gint j=0; j < N; j++){
            if (color[i][j] == 1 && id[i][j] != -1)
                ids.push_back(id[i][j]);
        }
    }

    for (gint i=0; i < ids.size(); i++)
        s.insert(ids[i]);
    ids.assign(s.begin(), s.end());

    return ids;
}

template <typename gfloat, typename gint>
summary<gfloat> RnB_similarity(vector<vector<gfloat> > &matrix, vector<vector<gfloat> > &red, vector<vector<gfloat> > &blue, vector<vector<gfloat> > &sim, gint &N){

    summary<gfloat> values;
    vector<gfloat> sim_sum; 
    vector<vector<gint> > id  = create_id_matrix(matrix, N);
    vector<gint> red_ids = link_community_ids(id, red, N);
    vector<gint> blue_ids = link_community_ids(id, blue, N);
    gfloat maxe = find_max_special(sim) + 1;

    for (auto ired=red_ids.cbegin(); ired != red_ids.cend(); ++ired){
       for (auto iblue=blue_ids.cbegin(); iblue != blue_ids.cend(); ++iblue){
            if (*ired > *iblue){
                sim_sum.push_back(sim[*iblue][*ired - *iblue - 1]);
            }
            else if (*ired < *iblue){
                sim_sum.push_back(sim[*ired][*iblue - *ired - 1]);
            } else{
                sim_sum.push_back(maxe);
            }
       }
    }

    if (sim_sum.size() > 0){
        values.mean = average(sim_sum);
        values.stdn = stdev(sim_sum, values.mean);
        // printf("mean= %Lf, std= %Lf\n", values.mean, values.stdn);
    }
    else{
        cout << "In RnB function. No division by zero\n";
        values.mean = -1;
        values.stdn = -1;
    }
        
    return values;
}

template <typename gfloat>
void print_matrix(vector<vector<gfloat> > &v){
    for (auto const& row : as_const(v)){
        for (auto const& i : as_const(row)){
            cout << i << " ";
        }
        cout << "\n";
    }
}

template <typename gfloat>
void print_vector(vector<gfloat> &v){
    for (auto const& i : as_const(v)){
        cout << fixed << setprecision(10) << i << " ";
    }
    cout << "\n";
}

int main(){

    //// Initializing variables

    // structs
    summary<long double> val;
    // vectors 2D
    vector<vector<long double> > red;
    vector<vector<long double> > blue;
    vector<vector<long double> > fln;
    vector<vector<long double> > sim;
    // vectors 1D
    vector<long double> vmean;
    vector<long double> vstd;
    // chars
    char fln_path[] = "../CSV/merged/imputation/tracto2016/zz_model/fln_4_r_6_3.csv";
    char sim_path[] = "../CSV/merged/similarity/tracto2016/zz_model/sim_l10_4_r_6_3.csv";
    char common[] = "../CSV/merged/RnB/tracto2016/zz_model/";
    char red_name[] = "blue_src_123456_tgt_123456_lcom_special.csv";
    char blue_name[] = "blue_src_123456_tgt_123456_lcom_special.csv";
    char red_path[sizeof(common)/sizeof(common[0])+sizeof(red_name)/sizeof(red_name[0])];
    char blue_path[sizeof(common)/sizeof(common[0])+sizeof(blue_name)/sizeof(blue_name[0])];
    // ints
    int_fast16_t N=107;
    // bools
    bool header;

    //// Start Code

    // load fln
    header = false;
    cout << "Loading fln\n";
    fln = read_csv(fln_path, header);
    cout << "Finished\n";

    // load similarities
    header = true;
    cout << "Loading similarities csv\n";
    sim = read_csv(sim_path, header);
    cout << "Finished\n";

    cout << "Formatting similarity matrix\n";
    nlog10_matrix(sim);
    replace_max(sim);
    cout << "Finished\n";

    // load red
    header = false;

    strcpy(red_path, common);
    strcat(red_path, red_name);
    printf("Loading red csv: %s\n", red_name);
    red = read_csv(red_path, header);
    cout << "Finished\n";

    // load blue
    header = false;
    strcpy(blue_path, common);
    strcat(blue_path, blue_name);
    printf("Loading blue csv: %s\n", blue_name);
    blue = read_csv(blue_path, header);
    cout << "Finished\n";

    // Compute distance
    cout << "Getting average similarity between RnB\n";
    val = RnB_similarity(fln, red, blue, sim, N);
    cout << "Finished\n";

    // Stack means and stds
    vmean.push_back(val.mean);
    vstd.push_back(val.stdn);

    // print maens and stds
    print_vector(vmean);
    print_vector(vstd);
    
    return 0;
}