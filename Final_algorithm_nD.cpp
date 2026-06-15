#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <numeric>
#include <queue>
#include <string>
#include <unordered_set>
#include <cmath>
#include <chrono>
#include <utility>
#include "nanoflann.hpp"
#include <cstdlib>
#include <ctime>
#include "KDTreeVectorOfVectorsAdaptor.h"
#include <map>
#include <list>
#include <set>
#include <boost/container_hash/hash.hpp>
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <unordered_map>

using namespace nanoflann;
using namespace std;
using namespace std::chrono;

constexpr size_t dim = 10;
static double epsilon = pow(10.0,-14.0);
typedef KDTreeVectorOfVectorsAdaptor<vector<vector<double>>, double> my_kd_tree_t;

struct compare_tuple {
    bool operator()(const tuple <size_t, size_t, double, size_t> & tuple_1, const tuple <size_t, size_t, double, size_t> & tuple_2)
    {
        if (get<2>(tuple_1) != get<2>(tuple_2)){
            return get<2>(tuple_1) > get<2>(tuple_2);
        } else {
            if (get<0>(tuple_1) != get<0>(tuple_2)){
                return get<0>(tuple_1) > get<0>(tuple_2);
            } else {
                return get<1>(tuple_1) > get<1>(tuple_2);
            }
        }
    }
};

double l2_dist(vector<double> &v_1, vector<double>  &v_2){
    double squ_vec [dim]; // This is where we keep the squared values of the numbers
    double squ_to_add; //dynamically make sq_to_add
    for (size_t i = 0; i < dim; ++i){
        squ_to_add = pow(v_1[i] - v_2[i], 2.0);
        squ_vec[i] = squ_to_add;
    }

    //Let's sum the elements on our own
    double tot_sum = 0.0;
    for (auto squ: squ_vec){
        tot_sum = tot_sum + squ;
    }

    double result = pow(tot_sum,0.5);
    return result;
}

double l2_dist_2(vector<double> & v_1, vector<double> & v_2){
    vector<double> v_8 (dim);
    for (size_t i = 0; i < dim; ++i){
        v_8[i] = v_1[i]-v_2[i];
    }
    double sq_norm = inner_product(v_8.begin(), v_8.end(), v_8.begin(), 0.0);
    return sq_norm;
}

bool compare_the_pair(pair<size_t,double> pair_1, pair<size_t,double> pair_2){

    return pair_1.second < pair_2.second;
}

bool compare_result_vector (nanoflann::ResultItem<size_t, double> result_vec_1, nanoflann::ResultItem<size_t, double> result_vec_2){
    return result_vec_1.second < result_vec_2.second;
}

vector<size_t> find_all_neighbors(size_t ver_idx, vector<vector<double>> &point_matrix, size_t n) { // This is a function that will find the list of neighbors of ver_idx.
    size_t number_of_idx = n - ver_idx - 1;
    vector<pair<size_t,double>> idx_dist_return (number_of_idx);
    for (size_t i = ver_idx + 1; i < n; ++i){
        idx_dist_return[i-ver_idx - 1].first = i;
        idx_dist_return[i-ver_idx - 1].second = l2_dist(point_matrix[ver_idx], point_matrix[i]);
    }

    // Now we need to sort the vector
    sort(idx_dist_return.begin(), idx_dist_return.end(), compare_the_pair);
    // now we need to return the appropiate vector
    vector<size_t> vec_to_return (number_of_idx);
    for (size_t j = 0; j < number_of_idx; ++j){
        vec_to_return[j] = idx_dist_return[j].first;
    }
    return vec_to_return;
}

//======================================= CODE USED FOR READING CSV ============================================
vector<vector<double>> read_csv(string filename){
    // Reads a CSV file into a vector of <string, vector<int>> pairs where
    // each pair represents <column name, column values>

    // Create a vector of <string, int vector> pairs to store the result
    vector<vector<double>> result;

    // Create an input filestream
    ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw runtime_error("Could not open file");

    // Helper vars
    string line;
    double val;

    // Read data, line by line
    while(getline(myFile, line))
    {
        // Create a stringstream of the current line
        stringstream ss(line);

        // Keep track of the current column index
        vector<double> new_vector;
        // Extract each integer
        while(ss >> val){
            // cout << val << endl;
            // Add the current integer to the 'colIdx' column's values vector
            new_vector.push_back(val);

            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();


        }
        if (! new_vector.empty()) {
            result.push_back(new_vector);
            new_vector.clear();
        }
    }

    // Close file
    myFile.close();

    return result;
}

//======================================== CODE USED TO CALCULATE NUMBER OF CONNECTED COMPONENTS ====================================

size_t UF_find(size_t i, vector<size_t> & id){
    while (i != id[i])
    {
        id[i] = id[id[i]];
        i = id[i];
    }
    return i;
}

void UF_union(size_t p, size_t q, vector<size_t> & id, vector<size_t> & sz){
    size_t i = UF_find(p,id);
    size_t j = UF_find(q,id);
    if (i == j){
        return;
    }
    if (sz[i] < sz[j])
    {
        id[i] = j;
        sz[j] += sz[i];
    }
    else
    {
        id[j] = i;
        sz[i] += sz[j];
    }
}

class Edge_2 {

public:

    size_t i, j;  // Indices of endpoints in the point set

    double length;  // L_p distance



    Edge_2(size_t i, size_t j, double len){
        if (i < j){
            this->i = i;
            this->j = j;
            this->length = len;
        } else {
            this->j = i;
            this->i = j;
            this->length = len;
        }
    }

    // Comparison operator for sorting edges by length

    bool operator<(const Edge_2& other) const {

        if (length != other.length)

            return length < other.length;

        if (i != other.i)

            return i < other.i;

        return j < other.j;

    }



    // Equality operator for set operations

    bool operator==(const Edge_2& other) const {

        return (i == other.i && j == other.j) || (i == other.j && j == other.i);

    }

};

class Point_2 {


public:

    // Constructor

    explicit Point_2(const std::vector<double>& coords) : coordinates(coords) {}



    // Get dimension

    size_t dim() const {

        return coordinates.size();

    }



    // Get a coordinate

    double operator[](size_t i) const {

        return coordinates[i];

    }



    // Print the point

    friend std::ostream& operator<<(std::ostream& os, const Point_2& p) {

        os << "(";

        for (size_t i = 0; i < p.dim(); ++i) {

            os << p[i];

            if (i < p.dim() - 1) {

                os << ", ";

            }

        }

        os << ")";

        return os;

    }

    std::vector<double> coordinates;
};

// Calculate L_p distance between two points

double Lp_distance(const Point_2& p1, const Point_2& p2, double p) {

    if (p1.dim() != p2.dim()) {

        throw std::invalid_argument("Point_2s must have the same dimension");

    }



    double sum = 0.0;

    for (size_t i = 0; i < p1.dim(); ++i) {

        sum += std::pow(std::abs(p1[i] - p2[i]), p);

    }



    return std::pow(sum, 1.0 / p);

}

//===========Find RNG code==============================================================================
set<Edge_2> RNG_algorithm(const std::vector<Point_2>& V, double p, my_kd_tree_t & kd_tree) {

    size_t n = V.size();

    std::set<Edge_2> E;  // Final set of edges (RNG)



    // Phase 1: Construct a supergraph of RNG(V)

    for (size_t i = 0; i < n; ++i) {

        std::set<Edge_2> Ei;  // Set of edges adjacent to vi



        // Create potential edges from vi to all other points

        for (size_t j = 0; j < n; ++j) {

            if (i != j) {

                double dist = Lp_distance(V[i], V[j], p);

                Ei.insert(Edge_2(i, j, dist));

            }

        }



        // Process edges in Ei

        while (!Ei.empty()) {

            // Get minimal length edge (first element in the set)

            auto min_edge_it = Ei.begin();

            Edge_2 min_edge = *min_edge_it;



            // Remove minimal edge from Ei

            Ei.erase(min_edge_it);



            // Add minimal edge to E
            // Need to check if it is in E
            if (E.count(min_edge) == 0){
                E.insert(min_edge);
            }



            // Remove edges from Ei according to the algorithm

            std::set<Edge_2> Ei_new;

            for (const auto& edge : Ei) {
                if (min_edge.i == edge.i){
                    Edge_2 edge_1 = Edge_2(min_edge.j, edge.j, Lp_distance(V[min_edge.j], V[edge.j], p));
                    Edge_2 edge_2 = Edge_2(i, edge.j, Lp_distance(V[i], V[edge.j], p));

                    if (edge_2 < edge_1) {

                        Ei_new.insert(edge);

                    }
                } else if (min_edge.i == edge.j){
                    Edge_2 edge_1 = Edge_2(min_edge.j, edge.i, Lp_distance(V[min_edge.j], V[edge.i], p));
                    Edge_2 edge_2 = Edge_2(i, edge.i, Lp_distance(V[i], V[edge.i], p));
                    if (edge_2 < edge_1) {
                        Ei_new.insert(edge);

                    }
                } else if (min_edge.j == edge.j) {
                    Edge_2 edge_1 = Edge_2(min_edge.i, edge.i, Lp_distance(V[min_edge.i], V[edge.i], p));
                    Edge_2 edge_2 = Edge_2(i, edge.i, Lp_distance(V[i], V[edge.i], p));
                    if (edge_2 < edge_1) {

                        Ei_new.insert(edge);

                    }
                } else { // (min_edge.j == edge.i)
                    Edge_2 edge_1 = Edge_2(min_edge.i, edge.j, Lp_distance(V[min_edge.i], V[edge.j], p));
                    Edge_2 edge_2 = Edge_2(i, edge.j, Lp_distance(V[i], V[edge.j], p));

                    if (edge_2 < edge_1) {
                        Ei_new.insert(edge);

                    }
                }
            }

            Ei = Ei_new;

        }

    }

    cout << E.size() << endl;

    // Phase 2: Eliminate edges that don't belong to RNG(V)

    std::set<Edge_2> RNG;

    for (const auto& edge : E) {

        bool should_keep = true;

        // Let's use kd trees
        size_t a = edge.i;
        size_t b = edge.j;

        // Need to define the search radius
        double search_radius = pow(Lp_distance(V[a], V[b],p),p) + epsilon;
        vector<nanoflann::ResultItem<size_t, double>> radius_result_vector_a;
        vector<nanoflann::ResultItem<size_t, double>> radius_result_vector_b;

        nanoflann::RadiusResultSet<double> radius_resultSet_a(search_radius, radius_result_vector_a);
        nanoflann::RadiusResultSet<double> radius_resultSet_b(search_radius, radius_result_vector_b);
        vector<double> v_a = V[a].coordinates;
        vector<double> v_b = V[b].coordinates;
        kd_tree.index->findNeighbors(radius_resultSet_a, &v_a[0]);
        kd_tree.index->findNeighbors(radius_resultSet_b, &v_b[0]);

        // Need to intersect the two sets
        sort(radius_result_vector_a.begin(), radius_result_vector_a.end(), compare_result_vector);
        sort(radius_result_vector_b.begin(), radius_result_vector_b.end(), compare_result_vector);
        size_t n_ra = radius_result_vector_a.size() - 1;
        size_t n_rb = radius_result_vector_b.size() - 1;
        vector<size_t> r_a(n_ra);
        vector<size_t> r_b(n_rb);
        for (size_t j = 0; j < n_ra; ++j) {
            r_a[j] = radius_result_vector_a[j + 1].first;
        }

        for (size_t j = 0; j < n_rb; ++j) {
            r_b[j] = radius_result_vector_b[j + 1].first;
        }

        // We need to sort the vectors or else the intersection won't work
        sort(r_a.begin(), r_a.end());
        sort(r_b.begin(), r_b.end());

        vector<size_t> r_ab;
        set_intersection(r_a.begin(), r_a.end(), r_b.begin(), r_b.end(), back_inserter(r_ab));
        // At this stage r_ab should contain the indices of the intersection, though we need to
        // actually test they are in the intersection.

        // Check if any point is in the lune of the edge

        for (size_t k: r_ab) {

            if (k != edge.i && k != edge.j) {
                Edge_2 edge_ik = Edge_2(edge.i, k, Lp_distance(V[edge.i], V[k], p));
                Edge_2 edge_jk = Edge_2(edge.j, k, Lp_distance(V[edge.j], V[k], p));
                if (edge_ik < edge && edge_jk < edge) {

                    should_keep = false;

                    break;

                }

            }

        }



        if (should_keep) {

            RNG.insert(edge);

        }

    }

    return RNG;

}

vector<vector<size_t>> find_connected_components(vector<size_t>& lune_ab, vector<Point_2> V, Edge_2 edge_ab, vector<Edge_2>& our_edges, double p){
    size_t n_ps = lune_ab.size();
    // Need to actually search for the connected components
    vector<size_t> id(n_ps);
    vector<size_t> sz(n_ps);
    for (size_t qq = 0; qq < n_ps; ++qq) {
        id[qq] = qq;
        sz[qq] = 1;
    }
    size_t end_flag = 0;
    if (!our_edges.empty()){
        for (auto edge: our_edges){
            if (end_flag == 1){
                break;
            }
            Edge_2 edge_compare = Edge_2(lune_ab[edge.i], lune_ab[edge.j], Lp_distance(V[lune_ab[edge.i]], V[lune_ab[edge.j]],2));
            if (edge_compare < edge_ab) {
                UF_union(edge.i, edge.j, id, sz);
                if (sz[edge.i] == n_ps or sz[edge.j] == n_ps){
                    end_flag = 1;
                    break;
                }
            }
        }

        if (end_flag == 1){
            // in this case there is only one connected component.
            // just return point_set
            sort(lune_ab.begin(), lune_ab.end());
            vector<vector<size_t>> return_vec;
            return_vec.push_back(lune_ab);
            return return_vec;
        } else {
            for (size_t zz = 0; zz < n_ps; ++zz) {
                id[zz] = UF_find(zz, id);
            }

            //Need id_copy to find the number of connected components
            vector<size_t> id_copy = id;
            sort(id_copy.begin(), id_copy.end());
            id_copy.erase(unique(id_copy.begin(), id_copy.end()), id_copy.end());
            size_t n_components = id_copy.size();
            vector<vector<size_t>> return_vec (n_components);
            for (size_t i = 0; i < n_ps; ++i){
                for (size_t j = 0; j < n_components; ++j){
                    if (id[i] == id_copy[j]){
                        return_vec[j].push_back(lune_ab[i]);
                    }
                }
            }

            for (size_t i = 0; i < n_components; ++i){
                sort(return_vec[i].begin(), return_vec[i].end());
            }

            return return_vec;
        }
    } else {
        // Each point of lune_ab is its own connected component
        size_t n_components = lune_ab.size();
        vector<vector<size_t>> return_vec (n_components);
        for (size_t j = 0; j < n_components; ++j){
            return_vec[j].push_back(lune_ab[j]);
        }
    }
}
//===========================================================================================================================================================
// code for diameter

double diam (vector<size_t> vec, vector<vector<double>> & point_matrix){
    if (vec.size() == 2){
        return l2_dist(point_matrix[vec[0]], point_matrix[vec[1]]);
    } else {
        double side_1 = l2_dist(point_matrix[vec[0]], point_matrix[vec[1]]);
        double side_2 = l2_dist(point_matrix[vec[1]], point_matrix[vec[2]]);
        double side_3 = l2_dist(point_matrix[vec[0]], point_matrix[vec[2]]);
        vector<double> side_vec = {side_1, side_2, side_3};
        double max_val = *max_element(side_vec.begin(), side_vec.end());
        return max_val;
    }
}

int main (int argc, char** argv){
            auto start = high_resolution_clock::now();
            priority_queue<tuple<size_t, size_t, double, size_t>, vector<tuple<size_t, size_t, double, size_t>>, compare_tuple> min_heap;
            size_t num_neighbors = stoi(argv[2]);
            // size_t num_neighbors = 100;
            const double p = 2.0;
            cout << argv[1] << endl;
            vector<vector<double>> point_matrix = read_csv(argv[1]);
            string filename = argv[1];
            // string file_name = "r10_solid_donut_" + to_string(point_cloud_size) + "_" + to_string(round) + ".csv";
            // vector<vector<double>> point_matrix = read_csv(file_name);
            size_t n = point_matrix.size();

            if (num_neighbors == 0) {
                num_neighbors = pow(n, 0.5);
            }

            cout << "the value of num_neighbors is " << num_neighbors << endl;

            // typedef KDTreeVectorOfVectorsAdaptor<vector<vector<double>>, double> my_kd_tree_t; // make things more readable.
            my_kd_tree_t mat_index(dim, point_matrix, 10 /* max leaf */);

            // Need to find RNG(X)
            // Need to create a Point_2 vector.
            vector<Point_2> Points_2;
            for (auto point: point_matrix) {
                Points_2.push_back(Point_2(point));
            }
            // Now we can find RNG
            set<Edge_2> RNG_X = RNG_algorithm(Points_2, 2.0, mat_index);
            cout << RNG_X.size() << endl;
            size_t n_RNG_edges_3 = 0;
            //populate both vectors

            unordered_set<pair<size_t, size_t>, boost::hash<pair<size_t, size_t>>> possible_edges_3;
            // We need a map from the points of point_matrix to their respective vertex



            cout << n_RNG_edges_3 << endl;
            size_t total_death = RNG_X.size() - n + 1;

            // Node counters
            vector<size_t> zero_vector{0, 0, 0};
            vector<double> zero_vector_3d{0.0, 0.0, 0.0};
            size_t node_ctr_1 = 0;
            size_t node_ctr_g_1 = 0;
            // Need to define the number of neighbors we are going to look for.
            // int k_n = 100; // number of neighbors we are going to look for
            vector<vector<size_t>> N(n); // array storing all the nearest neighbors.
            // There will be n rows, and each row will have variable length.

            // We need something to store the bars. A vector of vectors will do the trick.
            vector<vector<double>> barcode_bars;
            // We need to make the AVL trees we will use in the reduction
            // Node* root_1 = nullptr;
            // Node* root_g_1 = nullptr;

            // What were nodes will be maps now
            unordered_map<size_t, vector<size_t>> root_1;
            unordered_map<size_t, vector<size_t>> root_g_1;

            size_t column_counter = 0;
            // First we need to make the point set. Let's use the function we defined above.

            cout << "point matrix is" << endl;
            // We need to create two sets of two maps for the one and two-simplices.
            map<vector<size_t>, size_t> one_simp_to_idx;
            map<size_t, vector<size_t>> idx_to_one_simp;

            // We need to build the kd-tree
            // This line of code is just to make things a bit easier. All we are doing is using
            // declaring that we will write KDTreeVectorOfVectorsAdaptor<vector<vector<double>>, double>
            // as my_kd_tree_t.

            // I think this is the part where we actually make the kd-tree.
            // my_kd_tree_t mat_index(dim, point_matrix, 10 /* max leaf */);

            // Number of neighbors we will initially search for

            // vectors to store the output.
            {
                vector<size_t> ret_indexes(num_neighbors);
                vector<double> out_dists_sqr(num_neighbors);
                nanoflann::KNNResultSet<double> resultSet(num_neighbors);
                resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

                // We need to create a minimum heap
                // elements of the minimum heap will be of the form tuple <int_t, int_t, double, int_t>

                cout << "populating N..." << endl;

                // We need to work on populating N
                for (size_t i = 0; i < n; ++i) {
                    // Need to populate the vector N[i]
                    vector<double> query_pt = point_matrix[i];
                    vector<size_t> ret_indexes(num_neighbors);
                    vector<double> out_dists_sqr(num_neighbors);
                    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]); // Need to initialise the result set
                    mat_index.index->findNeighbors(resultSet, &query_pt[0]); //This finds the neighbors
                    // Now we need to find the indices greater than i. To this end we have
                    // the following.
                    bool smallest_flag = false;
                    for (size_t j = 0; j < num_neighbors; ++j) {
                        if (ret_indexes[j] > i) {
                            N[i].push_back(ret_indexes[j]);
                            // We might as well add the smallest length to the minimum heap at this stage.
                            if (!smallest_flag) {
                                min_heap.push(tuple<size_t, size_t, double, size_t>{i, ret_indexes[j],
                                                                                    double(pow(out_dists_sqr[j], 0.5)),
                                                                                    0});
                                smallest_flag = true;
                            }
                        }
                    }

                    if (!smallest_flag and i != n - 1) {
                        //need to find all neighbors
                        N[i] = find_all_neighbors(i, point_matrix, n);
                        min_heap.push(tuple<size_t, size_t, double, size_t>{i, N[i][0],
                                                                            l2_dist(point_matrix[i],
                                                                                    point_matrix[N[i][0]]),
                                                                            0});
                    }

                    // Need to clear these so that they can be used again.
                    ret_indexes.clear();
                    out_dists_sqr.clear();

                }
            }


            cout << "beginning main program..." << endl;
            size_t one_simp_ctr = 0; // one simplex counter
            size_t death_counter = 0;
            size_t two_simp_ctr = 0;
            // size_t two_simp_ctr = 0; // two simplex counter
            while (death_counter < total_death) {
                // cout << one_simp_ctr << endl;
                // cout << death_counter << endl;
                size_t min_heap_flag = 0;
                tuple<size_t, size_t, double, size_t> tuple_now = min_heap.top();
                min_heap.pop();
                size_t a = get<0>(tuple_now);
                // cout << a << endl;
                size_t b = get<1>(tuple_now);
                // cout << b << endl;
                vector<size_t> lune_gen_vec = {a, b};
                Edge_2 edge_ab = Edge_2(a, b, Lp_distance(Points_2[a], Points_2[b], p));
                sort(lune_gen_vec.begin(), lune_gen_vec.end());
                one_simp_to_idx[vector<size_t>{a, b}] = one_simp_ctr;
                idx_to_one_simp[one_simp_ctr] = vector<size_t>{a, b};
                one_simp_ctr += 1;
                double r = get<2>(tuple_now);
                // cout << r << endl;
                size_t t = get<3>(tuple_now);
                // cout << t << endl;
                if (t + 2 > size_t(N[a].size())) { // If true we need to update N[a]
                    // Need to find the list of neighbors of a
                    //cout << "DID THIS ACTIVATE" << endl;
                    N[a] = find_all_neighbors(a, point_matrix, n);
                    if (size_t(N[a].size()) < t + 2) {
                        min_heap_flag = 1;
                    }
                }

                if (min_heap_flag == 0) {
                    //cout << "I SUSPECT THE PROBLEM IS HERE" << endl;
                    size_t b_dash = N[a][t + 1];
                    //cout << "MAYBE IT DID NOT LIKE THIS PART" << endl;
                    double r_dash = l2_dist(point_matrix[a], point_matrix[b_dash]);
                    //cout << "MAYBE IT DID NOT LIKE THIS PART" << endl;
                    size_t t_dash = t + 1;

                    // Now we need to add the new tuple to the heap
                    min_heap.push(tuple<size_t, size_t, double, size_t>{a, b_dash, r_dash, t_dash});
                }
                // cout << "doing a radius search..." << endl;
                // do a radius search

                // We will implement a quick check to see if we can avoid a lot of computation.
                vector<double> mid_point(dim);

                //Fill in the entries of mid_point.
                for (size_t i = 0; i < dim; ++i) {
                    mid_point[i] = (point_matrix[a][i] + point_matrix[b][i]) / 2.0;
                }

                double r_small = pow(0.5 * (2 - pow(3.0, 0.5)) * r, 2.0);
                vector<nanoflann::ResultItem<size_t, double>> radius_result_vector_small;


                double search_radius = pow(r, 2.0) + epsilon;
                // vector<pair<size_t, double>> match_dist;
                vector<double> query_pt_a = point_matrix[a];
                vector<double> query_pt_b = point_matrix[b];

                nanoflann::RadiusResultSet<double> radius_resultSet_small(r_small, radius_result_vector_small);
                mat_index.index->findNeighbors(radius_resultSet_small, &mid_point[0]);

                if (radius_result_vector_small.empty()) {
                    vector<nanoflann::ResultItem<size_t, double>> radius_result_vector_a;
                    vector<nanoflann::ResultItem<size_t, double>> radius_result_vector_b;
                    nanoflann::RadiusResultSet<double> radius_resultSet_a(search_radius, radius_result_vector_a);
                    nanoflann::RadiusResultSet<double> radius_resultSet_b(search_radius, radius_result_vector_b);

                    // cout << "using nanoflann to find radius..." << endl;

                    // radius_resultSet.init();
                    mat_index.index->findNeighbors(radius_resultSet_a, &query_pt_a[0]);

                    // cout << "calculating radius again..." << endl;
                    mat_index.index->findNeighbors(radius_resultSet_b, &query_pt_b[0]);

                    // cout << "Used nanoflann to find radius..." << endl;

                    // We need to sort the radius result vectors
                    sort(radius_result_vector_a.begin(), radius_result_vector_a.end(), compare_result_vector);
                    sort(radius_result_vector_b.begin(), radius_result_vector_b.end(), compare_result_vector);

                    size_t n_ra = radius_result_vector_a.size() - 1;
                    size_t n_rb = radius_result_vector_b.size() - 1;

                    // cout << n_ra << endl;
                    // cout << n_rb << endl;
                    if ((n_ra != 0) and (n_rb != 0)) {

                        vector<size_t> r_a(n_ra);
                        vector<size_t> r_b(n_rb);

                        for (size_t j = 0; j < n_ra; ++j) {
                            r_a[j] = radius_result_vector_a[j + 1].first;
                        }

                        for (size_t j = 0; j < n_rb; ++j) {
                            r_b[j] = radius_result_vector_b[j + 1].first;
                        }
                        // We need to sort the vectors or else the intersection won't work
                        sort(r_a.begin(), r_a.end());
                        sort(r_b.begin(), r_b.end());
                        // At this stage we have r_a and r_b. Now we need to find the intersection of r_a and r_b
                        vector<size_t> r_ab;
                        set_intersection(r_a.begin(), r_a.end(), r_b.begin(), r_b.end(), back_inserter(r_ab));
                        r_a.clear();
                        r_b.clear();

                        // Now we are at the stage where we find the number of connected components
                        size_t n_rab_2 = r_ab.size();
                        for (size_t j = n_rab_2 - 1; j == 0; --j) {
                            vector<size_t> vec_to_check_a = {a, r_ab[j]};
                            vector<size_t> vec_to_check_b = {b, r_ab[j]};
                            sort(vec_to_check_a.begin(), vec_to_check_a.end());
                            sort(vec_to_check_b.begin(), vec_to_check_b.end());
                            double r_11 = l2_dist(point_matrix[a], point_matrix[r_ab[j]]);
                            double r_22 = l2_dist(point_matrix[b], point_matrix[r_ab[j]]);
                            if (r_11 <= r && r_22 <= r) {
                                // Now we need to continue to check these
                                if (r_11 < r && r_22 < r) {
                                    // In this case it is definitely in the lune
                                    // Thus we continue
                                    continue;
                                } else {
                                    if (r_11 > r - epsilon && r_22 < r) {
                                        // In this case we need to check if sort([a r_ab[j]]) < [a b]
                                        // Now we need to compare vec_to_check and lune_gen_vec
                                        if (vec_to_check_a < lune_gen_vec) {
                                            continue;
                                        } else {
                                            r_ab.erase(r_ab.begin() + j);
                                        }
                                    } else if (r_11 < r && r_22 > r - epsilon) {
                                        // We just do the same thing as we did above
                                        // In this case we need to check if sort([a r_ab[j]]) < [a b]
                                        // Now we need to compare vec_to_check and lune_gen_vec
                                        if (vec_to_check_b < lune_gen_vec) {
                                            continue;
                                        } else {
                                            r_ab.erase(r_ab.begin() + j);
                                        }
                                    } else { // r_11 == r_2 && r_22 == r_2
                                        // In this case we need to check both the vectors
                                        if (vec_to_check_a < lune_gen_vec && vec_to_check_b < lune_gen_vec) {
                                            continue;
                                        } else {
                                            r_ab.erase(r_ab.begin() + j);
                                        }
                                    }
                                }
                            } else {
                                r_ab.erase(r_ab.begin() + j);
                            }
                        }

                        // cout << n_rab << endl;
                        size_t n_rab = r_ab.size();

                        if (n_rab == 1) { // only one point in the lune
                            column_counter += 1;
                            vector<size_t> two_simp_to_add{a, b, r_ab[0]};
                            sort(two_simp_to_add.begin(), two_simp_to_add.end());
                            size_t p_1 = one_simp_to_idx[vector<size_t>{two_simp_to_add[0], two_simp_to_add[1]}];
                            size_t p_2 = one_simp_to_idx[vector<size_t>{two_simp_to_add[1], two_simp_to_add[2]}];
                            size_t p_3 = one_simp_to_idx[vector<size_t>{two_simp_to_add[0], two_simp_to_add[2]}];
                            vector<size_t> vector_to_add = {p_1, p_2, p_3};
                            sort(vector_to_add.begin(), vector_to_add.end());
                            size_t l_value = *max_element(vector_to_add.begin(), vector_to_add.end());
                            // root_1 = insertNode(root_1, l_value, vector_to_add);
                            root_1[l_value] = vector_to_add;
                            node_ctr_1 += 1;
                            two_simp_ctr += 1;
                            vector_to_add.clear();
                        } else if (n_rab > 1) {
                            size_t n_conn_comp = 0;
                            vector<size_t> conn_comp_vector_dash;
                            if (n_rab == 2) {
                                double r_r = l2_dist(point_matrix[r_ab[0]], point_matrix[r_ab[1]]);
                                if (r_r < r) { // one connected component
                                    n_conn_comp = 1;
                                } else {
                                    n_conn_comp = 2;
                                    conn_comp_vector_dash = {r_ab[0], r_ab[1]};
                                }
                            } else { // We have at least three points meaning it makes sense to have a delaunay triangulation
                                //populate both vectors

                                // Before finding all connected components in the lune. Let us run a quick test
                                bool eye_flag = false; //Need this to test if there is a point in the eye shape
                                for (size_t j = 0; j < n_rab; ++j) {
                                    // size_t temp_idx = r_ab[rand()%n_rab];
                                    size_t temp_idx = r_ab[j];
                                    // Need to check the angle of the temporary idx
                                    vector<double> temp_u(dim);
                                    vector<double> temp_v(dim);
                                    for (size_t k = 0; k < dim; ++k) {
                                        temp_u[k] = point_matrix[a][k] - point_matrix[temp_idx][k];
                                        temp_v[k] = point_matrix[b][k] - point_matrix[temp_idx][k];
                                    }
                                    double norm_u = l2_dist(temp_u, zero_vector_3d);
                                    double norm_v = l2_dist(temp_v, zero_vector_3d);
                                    double cos_value =
                                            inner_product(temp_u.begin(), temp_u.end(), temp_v.begin(), 0.0) /
                                            (norm_u * norm_v);
                                    if (cos_value < -0.8660254038) {
                                        eye_flag = true;
                                        break;
                                    }
                                }
                                if (not eye_flag) {
                                    // Here we need to actually compute the connected components.


                                    std::vector<Point<dim>> points_2(n_rab);

                                    // Fill up points_2 with the entries from points
                                    for (size_t i = 0; i < n_rab; ++i) {
                                        for (size_t j = 0; j < dim; ++j) {
                                            points_2[i][j] = point_matrix[r_ab[i]][j];
                                        }
                                    }
                                    KdTreeSolver<dim> solver(points_2);
                                    std::vector<Edge> edges = solver.get_solution();

                                    vector<Edge_2> our_edges; // We will keep the output and convert it into Edge_2

                                    for (auto edge: edges) {
                                        Edge_2 our_edge = Edge_2(r_ab[edge.first], r_ab[edge.second],
                                                                 Lp_distance(Points_2[r_ab[edge.first]],
                                                                             Points_2[r_ab[edge.second]], p));
                                        Edge_2 our_edge_to_add = Edge_2(edge.first, edge.second,
                                                                        Lp_distance(Points_2[edge.first],
                                                                                    Points_2[edge.second], p));
                                        if (our_edge < edge_ab) {
                                            // We add our edge to our_edges
                                            our_edges.push_back(our_edge_to_add);
                                        }
                                    }
                                    vector<vector<size_t>> connected_components;
                                    if (our_edges.empty()) {
                                        // Each point of lune_ab is its own connected component
                                        size_t n_components = r_ab.size();
                                        vector<vector<size_t>> return_vec(n_components);
                                        for (size_t j = 0; j < n_components; ++j) {
                                            return_vec[j].push_back(r_ab[j]);
                                        }
                                        connected_components = return_vec;
                                    } else {
                                        connected_components = find_connected_components(r_ab, Points_2, edge_ab,
                                                                                         our_edges,
                                                                                         p);
                                    }
                                    // Now we have the connected components
                                    // We can now proceed to reduction if there is more than one connected component.
                                    n_conn_comp = connected_components.size();

                                    for (size_t u = 0; u < n_conn_comp; ++u) {
                                        conn_comp_vector_dash.push_back(connected_components[u][0]);
                                    }

                                } else {
                                    n_conn_comp = 1;
                                }
                            }


                            if (n_conn_comp == 1) { // In this case we essentially have a copy of the case above
                                // cout << " ACTIVATING THIS PART OF CODE" << endl;
                                column_counter += 1;
                                vector<size_t> two_simp_to_add = {a, b, r_ab[0]};
                                sort(two_simp_to_add.begin(), two_simp_to_add.end());
                                size_t p_1 = one_simp_to_idx[vector<size_t>{two_simp_to_add[0], two_simp_to_add[1]}];
                                size_t p_2 = one_simp_to_idx[vector<size_t>{two_simp_to_add[1], two_simp_to_add[2]}];
                                size_t p_3 = one_simp_to_idx[vector<size_t>{two_simp_to_add[0], two_simp_to_add[2]}];
                                vector<size_t> vector_to_add = {p_1, p_2, p_3};
                                sort(vector_to_add.begin(), vector_to_add.end());
                                size_t l_value = *max_element(vector_to_add.begin(), vector_to_add.end());
                                // cout << "L VALUE IS" << l_value << endl;
                                // Now we need to add this node to the AVL tree.
                                root_1[l_value] = vector_to_add;
                                node_ctr_1 += 1;
                                two_simp_ctr += 1;
                                // cout << "value of node_ctr_1 is " << node_ctr_1 << endl;
                                vector_to_add.clear();
                                two_simp_to_add.clear();
                            } else {
                                for (size_t c = 0; c < n_conn_comp; ++c) {
                                    // cout << "WE HAVE MORE THAN ONE CONNECTED COMPONENT" << endl;
                                    column_counter += 1;
                                    vector<size_t> two_simp_to_add = {a, b, conn_comp_vector_dash[c]};
                                    sort(two_simp_to_add.begin(), two_simp_to_add.end());
                                    two_simp_ctr += 1;
                                    size_t p_1 = one_simp_to_idx[vector<size_t>{two_simp_to_add[0],
                                                                                two_simp_to_add[1]}];
                                    size_t p_2 = one_simp_to_idx[vector<size_t>{two_simp_to_add[1],
                                                                                two_simp_to_add[2]}];
                                    size_t p_3 = one_simp_to_idx[vector<size_t>{two_simp_to_add[0],
                                                                                two_simp_to_add[2]}];
                                    vector<size_t> vector_to_add = {p_1, p_2, p_3};
                                    sort(vector_to_add.begin(), vector_to_add.end());
                                    size_t l_value = *max_element(vector_to_add.begin(), vector_to_add.end());
                                    size_t flag_4 = 0;
                                    size_t flag_5 = 0;
                                    //cout << l_value << endl;
                                    while (flag_4 == 0 and flag_5 == 0) {
                                        if (root_1.count(l_value) != 0) {
                                            // cout << "SEARCHING BINARY TREE 1" << endl;
                                            vector<size_t> new_vector = root_1[l_value];
                                            vector<size_t> symm_diff_vec;
                                            set_symmetric_difference(new_vector.begin(), new_vector.end(),
                                                                     vector_to_add.begin(), vector_to_add.end(),
                                                                     back_inserter(symm_diff_vec));
                                            vector_to_add = symm_diff_vec;
                                            symm_diff_vec.clear();
                                            if (vector_to_add.size() > 0) {
                                                l_value = *max_element(vector_to_add.begin(), vector_to_add.end());
                                            } else {
                                                flag_5 = 1;
                                            }
                                            // sort(vector_to_add.begin(), vector_to_add.end());
                                        } else if (root_g_1.count(l_value) != 0) {
                                            // cout << "SEARCHING BINARY TREE 2" << endl;
                                            vector<size_t> new_vector = root_g_1[l_value];
                                            vector<size_t> symm_diff_vec;
                                            set_symmetric_difference(new_vector.begin(), new_vector.end(),
                                                                     vector_to_add.begin(), vector_to_add.end(),
                                                                     back_inserter(symm_diff_vec));
                                            vector_to_add = symm_diff_vec;
                                            symm_diff_vec.clear();
                                            if (vector_to_add.size() > 0) {
                                                l_value = *max_element(vector_to_add.begin(), vector_to_add.end());
                                            } else {
                                                flag_5 = 1;
                                            }
                                            // sort(vector_to_add.begin(), vector_to_add.end());
                                        } else {
                                            // cout << "SEARCHING OVER" << endl;
                                            flag_4 = 1; // The column is reduced.
                                        }
                                    }
                                    if (flag_5 == 0) {
                                        root_g_1[l_value] = vector_to_add;
                                        node_ctr_g_1 += 1;
                                        // cout << "value of node_ctr_g_1 is " << node_ctr_g_1 << endl;
                                        // At this stage (l_value, column_counter) is a persistent pair
                                        double diam_one_simp = diam(idx_to_one_simp[l_value], point_matrix);
                                        double diam_two_simp = diam(two_simp_to_add, point_matrix);
                                        // cout << diam_one_simp << endl;
                                        // cout << diam_two_simp << endl;
                                        if (diam_one_simp - diam_two_simp != 0) {
                                            barcode_bars.push_back(vector<double>{diam_one_simp, diam_two_simp});
                                            death_counter += 1;
                                        }
                                    }
                                    vector_to_add.clear();
                                    two_simp_to_add.clear();
                                }
                                // printTree(root_g_1, "", true);
                            }


                        } else { // n_rab = 0
                            // cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
                        }

                        // At this stage it is safe to clear r_ab
                        // At this stage conn_comp_vector should contain one point from each connected component.

                        r_ab.clear();
                    }


                } else {
                    column_counter += 1;
                    vector<size_t> two_simp_to_add = {a, b, radius_result_vector_small[0].first};
                    sort(two_simp_to_add.begin(), two_simp_to_add.end());
                    size_t p_1 = one_simp_to_idx[vector<size_t>{two_simp_to_add[0], two_simp_to_add[1]}];
                    size_t p_2 = one_simp_to_idx[vector<size_t>{two_simp_to_add[1], two_simp_to_add[2]}];
                    size_t p_3 = one_simp_to_idx[vector<size_t>{two_simp_to_add[0], two_simp_to_add[2]}];
                    vector<size_t> vector_to_add = {p_1, p_2, p_3};
                    sort(vector_to_add.begin(), vector_to_add.end());
                    size_t l_value = *max_element(vector_to_add.begin(), vector_to_add.end());
                    // cout << "L VALUE IS" << l_value << endl;
                    // Now we need to add this node to the AVL tree.
                    // root_1 = insertNode(root_1, l_value, vector_to_add);
                    root_1[l_value] = vector_to_add;
                    node_ctr_1 += 1;
                    // cout << "value of node_ctr_1 is " << node_ctr_1 << endl;
                    vector_to_add.clear();
                    two_simp_to_add.clear();
                }
            }

            sort(barcode_bars.rbegin(), barcode_bars.rend());
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);

            auto time_taken = double(duration.count() / 1000000);

            // write to file
            string output_filename = "barcode_for_" + filename;
            ofstream out;
            out.open(output_filename);
            for (auto bar: barcode_bars) {
                out << bar[0] << ", " << bar[1] << endl;
            }
            out << "========================" << endl;
            out << n << ", " << time_taken << ", " << one_simp_ctr << ", " << column_counter << ", " << death_counter << ", " << endl;
            out.close();
    return 0;
}
