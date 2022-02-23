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
        
                    while(getline(str, word, ','))
                        row.push_back(stold(word));
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

                while(getline(str, word, ','))
                    row.push_back(stold(word));
                content.push_back(row);
            }
        }
        else
            cout<<"Could not open the file\n";
    }

	return content;
}

template <typename T, typename I>
vector<vector<I> > create_id_matrix(vector<vector<T> > &matrix, I &N){
	
	vector<vector<I> > id_matrix(N, vector<I>(N,-1));
	I id = 0;

	for (I i=0; i < N; i++){
		for (I j=0; j < N; j++){
			if (matrix[i][j] > 0){
				id_matrix[i][j] = id;
				id++;
			}
		}
	}

	return id_matrix;
}

template <typename I>
struct node{
    I i, f;
};

template <typename I>
struct range
{
    node<I> x_range, y_range;
};

template <typename I>
range<I> get_block_ranges(node<I> &u, vector<I> &clust_size)
{   
    I range_x = 0;
    for (I i=0; i < u.i; i++)
    {
        range_x += clust_size[i];
    }
    node<I> range_block_u_x = {.i=range_x, .f=range_x+clust_size[u.i]};

    I range_y = 0;
    for (I i=0; i < u.f; i++)
    {
        range_y += clust_size[i];
    }
    node<I> range_block_u_y = {.i=range_y, .f=range_y+clust_size[u.f]};

    range<I> range_block = {.x_range=range_block_u_x, .y_range=range_block_u_y};

    return(range_block);
}

template <typename F, typename I>
long double block_similarity_uv(node<I> &u, node<I> &v, vector<vector<F> > &matrix, vector<vector<F> > &similarity, vector<I> &clust_size, I &N)
{
    vector<vector<I> > id_matrix = create_id_matrix(matrix, N);
    range<I> range_u = get_block_ranges(u, clust_size);
    range<I> range_v = get_block_ranges(v, clust_size);
    F jacp=0;
    I id_1, id_2, t=0;

    if (u.i == v.i && u.f != v.f)
    {
        for (I u_i=range_u.x_range.i; u_i < range_u.x_range.f; u_i++)
        {
            for (I u_j=range_u.y_range.i; u_j < range_u.y_range.f; u_j++)
            {
                id_1 = id_matrix[u_i][u_j];
                if (id_1 < 0)
                    continue;
                for (I v_j=range_v.y_range.i; v_j < range_v.y_range.f; v_j++)
                {
                    id_2 = id_matrix[u_i][v_j];
                    if (id_2 < 0)
                        continue; 
                    if (id_1 > id_2)
                    {
                        jacp += similarity[id_2][id_1 - id_2 - 1];
                        t++;
                    } else if (id_1 < id_2)
                    {
                        jacp += similarity[id_1][id_2 - id_1 - 1];
                        t++;
                    } else{
                        cout << "There is no similarity between the same edges\n";
                    }
                }
            }
        }
    } else if (u.i != v.i && u.f == v.f)
    {
        for (I u_j=range_u.y_range.i; u_j < range_u.y_range.f; u_j++)
        {
            for (I u_i=range_u.x_range.i; u_i < range_u.x_range.f; u_i++)
            {
                id_1 = id_matrix[u_i][u_j];
                if (id_1 < 0)
                    continue;
                for (I v_i=range_v.x_range.i; v_i < range_v.x_range.f; v_i++)
                {
                    id_2 = id_matrix[v_i][u_j];
                    if (id_2 < 0)
                        continue;
                    if (id_1 > id_2)
                    {   
                        jacp += similarity[id_2][id_1 - id_2 - 1];
                        t++;
                    } else if (id_1 < id_2)
                    {
                        jacp += similarity[id_1][id_2 - id_1 - 1];
                        t++;
                    } else{
                        cout << "There is no similarity between the same edges\n";
                    }
                }
            }
        }
    } else if (u.i == v.i && u.f == v.f)
    {
        // if (u.i != u.f){
        //     printf("%i %i\n",  u.i, u.f);
        //     printf("%i %i %i %i\n", range_u.x_range.i, range_u.x_range.f, range_u.y_range.i, range_u.y_range.f);
        // }

        for (I u_i=range_u.x_range.i; u_i < range_u.x_range.f; u_i++)
        {
            for (I u_j=range_u.y_range.i; u_j < range_u.y_range.f; u_j++)
            {
                id_1 = id_matrix[u_i][u_j];
                if (id_1 < 0)
                    continue;
                for (I v_j=u_j; v_j < range_u.y_range.f; v_j++)
                {
                    id_2 = id_matrix[u_i][v_j];
                    if (id_2 < 0)
                        continue;
                    if (id_1 > id_2)
                    {   
                        // if (u.i == 5 && u.f == 1 && u.i==v.i && u.f==v.f){
                        //     printf("%i %i %Lf\n", id_1, id_2, similarity[id_2][id_1 - id_2 - 1]);
                        // }
                        jacp += similarity[id_2][id_1 - id_2 - 1];
                        t++;
                    } else if (id_1 < id_2)
                    {
                        // if (u.i == 5 && u.f == 1 && u.i==v.i && u.f==v.f){
                        //     printf("%i %i %Lf\n", id_1, id_2, similarity[id_2][id_1 - id_2 - 1]);
                        // }
                        jacp += similarity[id_1][id_2 - id_1 - 1];
                        t++;
                    } 
                }
                for (I v_i=u_i; v_i < range_u.x_range.f; v_i++)
                {
                    id_2 = id_matrix[v_i][u_j];
                    if (id_2 < 0)
                        continue;
                    if (id_1 > id_2)
                    {
                        // if (u.i == 5 && u.f == 1 && u.i==v.i && u.f==v.f){
                        //     printf("%i %i %Lf\n", id_1, id_2, similarity[id_2][id_1 - id_2 - 1]);
                        // }
                        jacp += similarity[id_2][id_1 - id_2 - 1];
                        t++;
                    } else if (id_1 < id_2)
                    {
                        // if (u.i == 5 && u.f == 1 && u.i==v.i && u.f==v.f){
                        //     printf("%i %i %Lf\n", id_1, id_2, similarity[id_2][id_1 - id_2 - 1]);
                        // }
                        jacp += similarity[id_1][id_2 - id_1 - 1];
                        t++;
                    }
                }
            }
        }
    } else
    {
        cout << "Not acceptable node format\n";
    }

    if (t <= 0)
        cout << "No elements found. Division by zero\n";
    else
        jacp /= t;

    return jacp;
}

template <typename F, typename I>
struct combined_list{
    vector<vector<I> > nodes;
    vector<F> similarity;
};

template <typename F, typename I>
combined_list<F, I> block_similarity(vector<vector<F> > &matrix, vector<vector<F> > &similarity, vector<I> &clust_size, I &N)
{
    vector<vector<I> > node_list;
    vector<I> row;
    vector<F> sim_list;
    node<I> u, v;
    F sim;
    I n = clust_size.size();

    for (I i=0; i < n; i++)
    {
        u.i = i;
        for (I j=0; j < n; j++)
        {
            u.f = j;
            for (I k_i=i; k_i < n; k_i++)
            {
                v.i = i;
                v.f = k_i;
                sim = block_similarity_uv(u, v, matrix, similarity, clust_size, N);
                row.clear();
                row.push_back(i);
                row.push_back(j);
                row.push_back(i);
                row.push_back(k_i);
                node_list.push_back(row);
                sim_list.push_back(sim);

                printf("%i %i %i %i %Lf\n", i, j, i, k_i, sim);

            }
            for (I k_j=j; k_j < n; k_j++)
            {
                v.i = k_j;
                v.f = j;
                sim = block_similarity_uv(u, v, matrix, similarity, clust_size, N);
                row.clear();
                row.push_back(i);
                row.push_back(j);
                row.push_back(k_j);
                row.push_back(j);
                node_list.push_back(row);
                sim_list.push_back(sim);

                printf("%i %i %i %i %Lf\n", i, j, k_j, j, sim);
            }
        }
    }

    combined_list<F,I> list_r = {.nodes=node_list, .similarity=sim_list};
    return list_r;
}

template <typename F, typename S>
void out_csv_2d(vector<vector<F> > &matrix, S &out_path){

    int_fast16_t sx = matrix.size();
    int_fast16_t sy = matrix[0].size();
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

template <typename F, typename S>
void out_csv_1d(vector<F> &matrix, S &out_path){
	ofstream out(out_path);
	for (int i=0; i < matrix.size(); i++){
		out << matrix[i] << '\n';
	}
}
 
int main(){
    string out_path;
    int N = 107;

    string csv_path = "../CSV/merged/similarity/tracto2016/zz_model/sim_l10_4_r_6_3.csv"; // merged/similarity/tracto2016/zz_model/
    bool header = true;
    vector<vector<long double> > sim = read_csv(csv_path, header);

    csv_path = "../CSV/merged/imputation/tracto2016/zz_model/fln_4_r_6_3.csv";
    header = false;
    vector<vector<long double> > fln = read_csv(csv_path, header);

    vector<int> sizes = {10, 17, 15, 5, 51, 9};

    node<int> u = {.i=0, .f=0};
    node<int> v = {.i=0, .f=0};


    combined_list<long double, int> block_sims = block_similarity(fln, sim, sizes, N);
    out_path="../CSV/merged/similarity/tracto2016/zz_model/nodes_block_sim_4_r_6_3.csv";
    out_csv_2d(block_sims.nodes, out_path);
    out_path="../CSV/merged/similarity/tracto2016/zz_model/sim_block_sim_4_r_6_3.csv";
    out_csv_1d(block_sims.similarity, out_path);

    return 0; 
}