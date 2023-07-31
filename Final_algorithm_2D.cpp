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
#include <boost/container_hash/hash.hpp>
#include "delaunator.hpp"
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

using namespace nanoflann;
using namespace std;
using namespace std::chrono;

static size_t dim = 2;

// Function to generate random point cloud. (From nanoflann)
void generateRandomPointCloud(vector<vector<double>>& point_matrix, const size_t n_points, const size_t dimension, const double max_range = 10.0){
    // random_device rd;
    // uniform_real_distribution<double> dist(0.0,max_range);
    std::cout << "Generating " << n_points << " random points...";
    point_matrix.resize(n_points);
    for (size_t i = 0; i < n_points; i++){
        point_matrix[i].resize(dimension);
        for (size_t d = 0; d < dim; d++){
            point_matrix[i][d] = max_range *(rand() % 1000000) / 1000000;
        }
    }
    std::cout << "done\n";
}


struct compare_tuple {
    bool operator()(const tuple <size_t, size_t, double, size_t> & tuple_1, const tuple <size_t, size_t, double, size_t> & tuple_2)
    {
        return get<2>(tuple_1) > get<2>(tuple_2);
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
            cout << val << endl;
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

// Graph class represents a undirected graph
// using adjacency list representation
class Graph {
    // No. of vertices
    size_t V;

    // Pointer to an array containing adjacency lists
    list<size_t>* adj;

    // A function used by DFS
    void DFSUtil(size_t v, bool visited[]);

public:
    // Constructor
    explicit Graph(size_t V);

    void addEdge(size_t v, size_t w);
    vector<size_t> NumberOfconnectedComponents();
};

// Function to return the number of
// connected components in an undirected graph
vector<size_t> Graph::NumberOfconnectedComponents(){


    // Vector to store the single point of each connected component
    vector<size_t> conn_comp_points;
    // Mark all the vertices as not visited
    bool* visited = new bool[V];

    // To store the number of connected components
    size_t count = 0;
    for (size_t v = 0; v < V; v++)
        visited[v] = false;

    for (size_t v = 0; v < V; v++) {
        if (!visited[v]) { // This means you have a new connected component.
            conn_comp_points.push_back(v);
            DFSUtil(v, visited);
            count += 1;
        }
    }
    delete [] visited;
    return conn_comp_points;
}

void Graph::DFSUtil(size_t v, bool visited[])
{

    // Mark the current node as visited
    visited[v] = true;

    // Recur for all the vertices
    // adjacent to this vertex
    list<size_t>::iterator i;

    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
            DFSUtil(*i, visited);
}

Graph::Graph(size_t V){
    this->V = V;
    adj = new list<size_t>[V];
}

// Add an undirected edge
void Graph::addEdge(size_t v, size_t w){
    adj[v].push_back(w);
    adj[w].push_back(v);
}

//===========AVL TREE CODE==============================================================================
class Node {
public:
    size_t key; // This will be the l value
    Node *left;
    Node *right;
    int height;
    vector<size_t> column_vector; // This will be the column vector which stores what we want.
};

int max(int a, int b);

// Calculate height
int height(Node *N) {
    if (N == nullptr)
        return 0;
    return N->height;
}

int max(int a, int b) {
    return (a > b) ? a : b;
}

// New node creation
Node *newNode(size_t key, vector<size_t> & col_vec) {
    Node *node = new Node();
    node->key = key;
    node->left = nullptr;
    node->right = nullptr;
    node->height = 1;
    node->column_vector = col_vec;
    return (node);
}

// Rotate right
Node *rightRotate(Node *y) {
    Node *x = y->left;
    Node *T2 = x->right;
    // cout << "identified left anf right nodd" << endl;
    x->right = y;
    y->left = T2;
    y->height = max(height(y->left),
                    height(y->right)) +
                1;
    x->height = max(height(x->left),
                    height(x->right)) +
                1;
    return x;
}

// Rotate left
Node *leftRotate(Node *x) {
    Node *y = x->right;
    Node *T2 = y->left;
    y->left = x;
    x->right = T2;
    x->height = max(height(x->left),
                    height(x->right)) +
                1;
    y->height = max(height(y->left),
                    height(y->right)) +
                1;
    return y;
}

// Get the balance factor of each node
int getBalanceFactor(Node *N) {
    if (N == nullptr)
        return 0;
    return height(N->left) -
           height(N->right);
}

// Insert a node
Node *insertNode(Node *node, size_t key, vector<size_t> col_vec) {
    // Find the correct postion and insert the node
    if (node == nullptr){
        // cout << "creating root node..." << endl;
        return (newNode(key, col_vec));
    }
    if (key < node->key){
        // cout << "moving to left node..." << endl;
        node->left = insertNode(node->left, key, col_vec);
    }else if (key > node->key){
        // cout << "moving to right node..." << endl;
        node->right = insertNode(node->right, key, col_vec);
    } else {
        return node;
    }

    // Update the balance factor of each node and
    // balance the tree
    // cout << "updating balance factor" << endl;
    node->height = 1 + max(height(node->left),height(node->right));
    // cout << "getting balance factor" << endl;
    int balanceFactor = getBalanceFactor(node);
    // cout << "balance factor found" << endl;
    // cout << balanceFactor << endl;
    if (balanceFactor > 1) {
        if (key < node->left->key) {
            // cout << "performaing right rotate" << endl;
            return rightRotate(node);
        } else if (key > node->left->key) {
            // cout << "performing left right rotate" << endl;
            node->left = leftRotate(node->left);
            // cout << "performaing left rotate" << endl;
            return rightRotate(node);
        }
    }
    if (balanceFactor < -1) {
        if (key > node->right->key) {
            return leftRotate(node);
        } else if (key < node->right->key) {
            node->right = rightRotate(node->right);
            return leftRotate(node);
        }
    }
    return node;
}

// Node with minimum value
Node *nodeWithMimumValue(Node *node) {
    Node *current = node;
    while (current->left != nullptr)
        current = current->left;
    return current;
}

// Find a node
vector<size_t> FindNode(Node *root, size_t key) {
    // Find the node and delete it
    if (root == nullptr){
        return vector<size_t> {0,0,0};
    } else {
        Node* root_examined = root;
        int flag6 = 0; // When the flag is one you are done
        while (flag6 == 0){
            if (key < root_examined->key){
                root_examined = root_examined->left;
                if (root_examined == nullptr){
                    flag6 = 1;
                    return vector<size_t> {0,0,0};
                }
            } else if (key > root_examined->key) {
                root_examined = root_examined->right;
                if (root_examined == nullptr){
                    flag6 = 1;
                    return vector<size_t> {0,0,0};
                }
            } else {
                flag6 = 1;
                return root_examined-> column_vector;
            }
        }

    }
}

// Print the tree
void printTree(Node *root, string indent, bool last) {
    if (root != nullptr) {
        cout << indent;
        if (last) {
            cout << "R----";
            indent += "   ";
        } else {
            cout << "L----";
            indent += "|  ";
        }
        cout << root->key << endl;
        printTree(root->left, indent, false);
        printTree(root->right, indent, true);
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

int main (){
    auto start = high_resolution_clock::now();
    //First we start off with basic information needed for both computing RNG and persistence

    //size_t n = 100; // number of points

    double max_range = 10.0; //used in generating max point cloud
    const size_t num_neighbors = 10; ************* CHANGE NUMBER O CONNECTED COMPONENTS. 
    vector<vector<double>> point_matrix = read_csv("point_cloud_uniform.csv"); *****CHANGE FILE HERE 
    size_t n = point_matrix.size();
    cout << "point matrix successfully made" << endl;
    for (int i = 0; i < 100; ++i){
    	for (int j = 0; j < 2; ++j){
    		if (j == 0){
    			cout << point_matrix[i][j] << ",";
    		} else {
    			cout << point_matrix[i][j] << endl;
    		}
    	}
    }
    //generateRandomPointCloud(point_matrix, n, dim, max_range); // Generate the point cloud
    typedef KDTreeVectorOfVectorsAdaptor<vector<vector<double>>, double> my_kd_tree_t; // make things more readable.
    // Compute the number of bars in the barcode by computing the RNG
    //if (dim == 2) { // This is what we will do if the point cloud is dimension 2.

    //need a set for the unordered edges
    unordered_set<pair<int, int>, boost::hash<pair<int, int>>> possible_edges;
    /* x0, y0, x1, y1, ... */
    vector<double> coords(2 * n); // flattened version of point_matrix

    //populate both vectors
    for (size_t i = 0; i < 2 * n; ++i) {
        if (i % 2 == 0) {
            coords[i] = point_matrix[size_t(i / 2)][0];
        } else {
            coords[i] = point_matrix[size_t((i - 1) / 2)][1];
        }
    }

    //code for displaying the point cloud

    my_kd_tree_t mat_index(dim, point_matrix, 10 /* max leaf */);
    //triangulation happens here
    delaunator::Delaunator d(coords);
    for (size_t i = 0; i < d.triangles.size(); i += 3) {
        if (d.triangles[i] < d.triangles[i + 1]) {
            possible_edges.insert(pair<int, int>{d.triangles[i], d.triangles[i + 1]});
        } else {
            possible_edges.insert(pair<int, int>{d.triangles[i + 1], d.triangles[i]});
        }

        if (d.triangles[i + 1] < d.triangles[i + 2]) {
            possible_edges.insert(pair<int, int>{d.triangles[i + 1], d.triangles[i + 2]});
        } else {
            possible_edges.insert(pair<int, int>{d.triangles[i + 2], d.triangles[i + 1]});
        }

        if (d.triangles[i] < d.triangles[i + 2]) {
            possible_edges.insert(pair<int, int>{d.triangles[i], d.triangles[i + 2]});
        } else {
            possible_edges.insert(pair<int, int>{d.triangles[i + 2], d.triangles[i]});
        }


    }

    // From possible_edges obtain the edges of the RNG
    auto n_RNG_edges = size_t(possible_edges.size());
    for (auto edge: possible_edges) {
        int a = edge.first;
        int b = edge.second;
        vector<double> query_pt_a = point_matrix[a];
        vector<double> query_pt_b = point_matrix[b];
        double r_2 = pow(coords[2 * a] - coords[2 * b], 2.0) + pow(coords[2 * a + 1] - coords[2 * b + 1], 2.0);
        vector<nanoflann::ResultItem<size_t, double>> radius_result_vector_a;
        vector<nanoflann::ResultItem<size_t, double>> radius_result_vector_b;
        nanoflann::RadiusResultSet<double> radius_resultSet_a(r_2, radius_result_vector_a);
        nanoflann::RadiusResultSet<double> radius_resultSet_b(r_2, radius_result_vector_b);
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


            // cout << "created vectors for intersection..." << endl;

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

            r_a.clear();
            r_b.clear();

            if (!r_ab.empty()) {
                n_RNG_edges -= 1;
            }
        }
    }
    // }


    cout << n_RNG_edges << endl;
    size_t total_death = n_RNG_edges - n + 1;

    // Node counters
    vector<size_t> zero_vector {0,0,0};
    vector<size_t> order_added_1;
    size_t node_ctr_1 = 0;
    size_t node_ctr_g_1 = 0;
    // Need to define the number of neighbors we are going to look for.
    // int k_n = 100; // number of neighbors we are going to look for
    vector<vector<size_t>> N (n); // array storing all the nearest neighbors.
    // There will be n rows, and each row will have variable length.

    // We need something to store the bars. A vector of vectors will do the trick.
    vector<vector<double>> barcode_bars;
    // We need to make the AVL trees we will use in the reduction
    Node* root_1 = nullptr;
    Node* root_g_1 = nullptr;

    size_t column_counter = 0;
    // First we need to make the point set. Let's use the function we defined above.

    cout << "point matrix is" << endl;
    /*
    for (size_t i = 0; i < n; ++i){
        for (size_t j = 0; j < dim; ++j){
            if (j == dim - 1){
                cout << point_matrix[i][j];
            } else {
                cout << point_matrix[i][j] << ", ";
            }
        }
        cout << "" << endl;
    }
    */


    // We need to create two sets of two maps for the one and two-simplices.
    map<vector<size_t>, size_t> one_simp_to_idx;
    map<size_t, vector<size_t>> idx_to_one_simp;
    map<vector<size_t>, size_t> two_simp_to_idx;
    map<size_t, vector<size_t>> idx_to_two_simp;

    // We need to build the kd-tree
    // This line of code is just to make things a bit easier. All we are doing is using
    // declaring that we will write KDTreeVectorOfVectorsAdaptor<vector<vector<double>>, double>
    // as my_kd_tree_t.

    // I think this is the part where we actually make the kd-tree.
    // my_kd_tree_t mat_index(dim, point_matrix, 10 /* max leaf */);

    // Number of neighbors we will initially search for

    // vectors to store the output.
    vector<size_t> ret_indexes(num_neighbors);
    vector<double> out_dists_sqr(num_neighbors);
    nanoflann::KNNResultSet<double> resultSet(num_neighbors);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

    // We need to create a minimum heap
    // elements of the minimum heap will be of the form tuple <int_t, int_t, double, int_t>
    priority_queue<tuple<size_t, size_t, double, size_t>, vector<tuple<size_t, size_t, double, size_t>>, compare_tuple> min_heap;

    cout << "populating N..." << endl;

    // We need to work on populating N
    for (size_t i = 0; i < n; ++i){
        // Need to populate the vector N[i]
        vector<double> query_pt = point_matrix[i];
        vector<size_t> ret_indexes(num_neighbors);
        vector<double> out_dists_sqr(num_neighbors);
        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]); // Need to initialise the result set
        mat_index.index->findNeighbors(resultSet, &query_pt[0]); //This finds the neighbors
        // Now we need to find the indices greater than i. To this end we have
        // the following.
        bool smallest_flag = false;
        for (size_t j = 0; j < num_neighbors; ++j){
            if (ret_indexes[j] > i){
                N[i].push_back(ret_indexes[j]);
                // We might as well add the smallest length to the minimum heap at this stage.
                if (!smallest_flag){
                    min_heap.push(tuple <size_t, size_t, double, size_t> {i, ret_indexes[j], double(pow(out_dists_sqr[j],0.5)), 0});
                    smallest_flag = true;
                }
            }
        }
        // Need to clear these so that they can be used again.
        ret_indexes.clear();
        out_dists_sqr.clear();
    }
    // Display N
    /*
    for (size_t w = 0; w < n; ++w){
        for (size_t x : N[w]){
            cout << x << ", " ;
        }
        cout << "" << endl;
    }
    */

    // Display minimum heap

    // Here we find k_d. The number of barcodes we expect to find.
    // For now just set k_d to some number
    // size_t k_d = ???
    // We need to use the other program in order to get this number.
    // For now we will just check a certain number of edges like in the matlab program.

    // size_t n_edges_checked = 2000;

    // Now we can begin with the main program

    cout << "beginning main program..." << endl;
    // ****************MAIN PROGRAM STARTS HERE***************
    size_t one_simp_ctr = 0; // one simplex counter
    size_t death_counter = 0;
    // size_t two_simp_ctr = 0; // two simplex counter
    while (death_counter < total_death){
        size_t min_heap_flag = 0;
        // cout << i << endl;
        tuple<size_t, size_t, double, size_t> tuple_now = min_heap.top();
        min_heap.pop();
        size_t a = get<0>(tuple_now);
        // cout << a << endl;
        size_t b = get<1>(tuple_now);
        // cout << b << endl;
        one_simp_to_idx[vector<size_t> {a,b}] = one_simp_ctr;
        idx_to_one_simp[one_simp_ctr] = vector<size_t> {a,b};
        one_simp_ctr += 1;
        double r = get<2>(tuple_now);
        // cout << r << endl;
        size_t t = get<3>(tuple_now);
        // cout << t << endl;
        if (t + 2 > size_t(N[a].size())){ // If true we need to update N[a]
            // Need to find the list of neighbors of a
            //cout << "DID THIS ACTIVATE" << endl;
            N[a] = find_all_neighbors(a, point_matrix,n);
            if (size_t(N[a].size()) < t + 2){
                min_heap_flag = 1;
            }
        }

        if (min_heap_flag == 0){
            //cout << "I SUSPECT THE PROBLEM IS HERE" << endl;
            size_t b_dash = N[a][t + 1];
            //cout << "MAYBE IT DID NOT LIKE THIS PART" << endl;
            double r_dash = l2_dist(point_matrix[a],point_matrix[b_dash]);
            //cout << "MAYBE IT DID NOT LIKE THIS PART" << endl;
            size_t t_dash = t + 1;

            // Now we need to add the new tuple to the heap
            min_heap.push(tuple<size_t, size_t, double, size_t> {a, b_dash, r_dash, t_dash});
        }
        // cout << "doing a radius search..." << endl;
        // do a radius search
        double search_radius = pow(r,2.0);
        // vector<pair<size_t, double>> match_dist;
        vector<double> query_pt_a = point_matrix[a];
        vector<double> query_pt_b = point_matrix[b];


        // cout << "creating nanoflann structures..." << endl;

        vector<nanoflann::ResultItem<size_t, double>> radius_result_vector_a;
        vector<nanoflann::ResultItem<size_t, double>> radius_result_vector_b;
        nanoflann::RadiusResultSet<double> radius_resultSet_a(search_radius,radius_result_vector_a);
        nanoflann::RadiusResultSet<double> radius_resultSet_b(search_radius,radius_result_vector_b);

        // cout << "using nanoflann to find radius..." << endl;

        // radius_resultSet.init();
        mat_index.index-> findNeighbors(radius_resultSet_a, &query_pt_a[0]);

        // cout << "calculating radius again..." << endl;
        mat_index.index-> findNeighbors(radius_resultSet_b, &query_pt_b[0]);

        // cout << "Used nanoflann to find radius..." << endl;

        // We need to sort the radius result vectors
        sort(radius_result_vector_a.begin(), radius_result_vector_a.end(), compare_result_vector);
        sort(radius_result_vector_b.begin(), radius_result_vector_b.end(), compare_result_vector);


        size_t n_ra = radius_result_vector_a.size() - 1;
        size_t n_rb = radius_result_vector_b.size() - 1;

        // cout << n_ra << endl;
        // cout << n_rb << endl;
        if ((n_ra != 0) and (n_rb != 0)){

            vector<size_t> r_a (n_ra);
            vector<size_t> r_b (n_rb);


            // cout << "created vectors for intersection..." << endl;

            for (size_t j = 0; j < n_ra; ++j){
                r_a[j] = radius_result_vector_a[j+1].first;
            }

            for (size_t j = 0; j < n_rb; ++j){
                r_b[j] = radius_result_vector_b[j+1].first;
            }
            // We need to sort the vectors or else the intersection won't work
            sort(r_a.begin(), r_a.end());
            sort(r_b.begin(), r_b.end());
            // At this stage we have r_a and r_b. Now we need to find the intersection of r_a and r_b

            // cout << "intersecting the vectors" << endl;
            /*size_t n_rab;
            if (n_ra >= n_rb){
                n_rab = n_ra;
            } else {
                n_rab = n_rb;
            }
            */
            vector<size_t> r_ab;
            /*
            cout << "-----------------------" << endl;
            for (size_t u = 0; u < n_ra; ++u){
                cout << r_a[u] << ", " ;
            }
            cout << "" << endl;

            for (size_t u = 0; u < n_rb; ++u){
                cout << r_b[u] << ", " ;
            }
            cout << "" << endl;
            cout << "------------------------" << endl;
            */
            set_intersection(r_a.begin(), r_a.end(), r_b.begin(), r_b.end(), back_inserter(r_ab));

            //At this stage r_ab should contain the intersection between a and b.

            // Print the contents of r_ab for testing purposes

            /*
            cout << "***********************************" << endl;
            if (r_ab.size() > 0){
                for (size_t z = 0; z < size_t(r_ab.size()); ++z){
                    cout << r_ab[z] << endl;
                }
                cout << "***********************************" << endl;
            }
            */
            r_a.clear();
            r_b.clear();

            // Now we are at the stage where we find the number of connected components
            size_t n_rab = r_ab.size();

            // cout << n_rab << endl;

            if (n_rab == 1){ // only one point in the lune
                column_counter += 1;
                vector<size_t> two_simp_to_add {a,b,r_ab[0]};
                sort(two_simp_to_add.begin(), two_simp_to_add.end());
                two_simp_to_idx[two_simp_to_add] = column_counter;
                idx_to_two_simp[column_counter] = two_simp_to_add;
                size_t p_1 = one_simp_to_idx[ vector<size_t> {two_simp_to_add[0],two_simp_to_add[1]}];
                size_t p_2 = one_simp_to_idx[ vector<size_t> {two_simp_to_add[1],two_simp_to_add[2]}];
                size_t p_3 = one_simp_to_idx[ vector<size_t> {two_simp_to_add[0],two_simp_to_add[2]}];
                //cout << p_1 << endl;
                //cout << p_2 << endl;
                //cout << p_3 << endl;
                vector<size_t> vector_to_add = {p_1, p_2, p_3};
                sort(vector_to_add.begin(), vector_to_add.end());
                //cout << vector_to_add.size()<< endl;
                size_t l_value = *max_element(vector_to_add.begin(), vector_to_add.end());
                //cout << l_value << endl;
                // printTree(root_1, "", true);
                /*
                for (size_t i_2 = 0; i_2 < order_added_1.size(); ++i_2){
                    cout << order_added_1[i_2] << ", ";
                }
                cout << "" << endl;
                */
                // cout << "----------------------------" << endl;
                // Now we need to add this node to the AVL tree.
                // cout << "What about here?" << endl;
                root_1 = insertNode(root_1, l_value, vector_to_add);
                // printTree(root_1, "", true);
                order_added_1.push_back(l_value);
                node_ctr_1 += 1;
                // cout << "value of node_ctr_1 is " << node_ctr_1 << endl;
                // cout << "is it here?" << endl;
                // cout << r_ab[0] << endl;
                // cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^";
                vector_to_add.clear();
            } else if (n_rab > 1){
                // cout << "making graph..." << endl;
                Graph prox_graph(n_rab);
                // We need to add the edges to the graph
                for (size_t k = 1; k < n_rab; ++k){
                    for (size_t y = 0; y < k; ++y){
                        double r_r = l2_dist(point_matrix[r_ab[y]], point_matrix[r_ab[k]]);
                        if (r_r < r){
                            prox_graph.addEdge(y,k);
                            // cout << "edge added..." << endl;
                        }
                    }
                }
                // cout << "finding connected components..." << endl;
                vector<size_t> conn_comp_vector_dash = prox_graph.NumberOfconnectedComponents(); // actually gives you a vector
                // with a point from each connected component.
                // need to change the indices to their corresponding values in r_ab.
                // cout << "we are ok" << endl;
                size_t n_conn_comp = conn_comp_vector_dash.size();

                // cout << "finding conn_comp_vector..." << endl;
                for (size_t k = 0; k < n_conn_comp; ++k){
                    size_t replace_value = conn_comp_vector_dash[k];
                    conn_comp_vector_dash[k] = r_ab[replace_value];
                    //cout << conn_comp_vector_dash[k] << ", ";
                }
                //cout << "" << endl;
                // cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^";
                if (n_conn_comp == 1){ // In this case we essentially have a copy of the case above
                    // cout << " ACTIVATING THIS PART OF CODE" << endl;
                    column_counter += 1;
                    vector<size_t> two_simp_to_add = {a,b,r_ab[0]};
                    sort(two_simp_to_add.begin(), two_simp_to_add.end());
                    two_simp_to_idx[two_simp_to_add] = column_counter;
                    idx_to_two_simp[column_counter] = two_simp_to_add;
                    size_t p_1 = one_simp_to_idx[ vector<size_t> {two_simp_to_add[0],two_simp_to_add[1]}];
                    size_t p_2 = one_simp_to_idx[ vector<size_t> {two_simp_to_add[1],two_simp_to_add[2]}];
                    size_t p_3 = one_simp_to_idx[ vector<size_t> {two_simp_to_add[0],two_simp_to_add[2]}];
                    vector<size_t> vector_to_add = {p_1, p_2, p_3};
                    sort(vector_to_add.begin(), vector_to_add.end());
                    size_t l_value = *max_element(vector_to_add.begin(), vector_to_add.end());
                    // cout << "L VALUE IS" << l_value << endl;
                    // Now we need to add this node to the AVL tree.
                    root_1 = insertNode(root_1, l_value, vector_to_add);
                    order_added_1.push_back(l_value);
                    node_ctr_1 += 1;
                    // cout << "value of node_ctr_1 is " << node_ctr_1 << endl;
                    vector_to_add.clear();
                    two_simp_to_add.clear();
                } else {
                    for (size_t c = 0; c < n_conn_comp; ++c){
                        // cout << "WE HAVE MORE THAN ONE CONNECTED COMPONENT" << endl;
                        column_counter += 1;
                        vector<size_t> two_simp_to_add = {a,b,conn_comp_vector_dash[c]};
                        sort(two_simp_to_add.begin(), two_simp_to_add.end());
                        two_simp_to_idx[two_simp_to_add] = column_counter;
                        idx_to_two_simp[column_counter] = two_simp_to_add;
                        size_t p_1 = one_simp_to_idx[ vector<size_t> {two_simp_to_add[0],two_simp_to_add[1]}];
                        size_t p_2 = one_simp_to_idx[ vector<size_t> {two_simp_to_add[1],two_simp_to_add[2]}];
                        size_t p_3 = one_simp_to_idx[ vector<size_t> {two_simp_to_add[0],two_simp_to_add[2]}];
                        vector<size_t> vector_to_add = {p_1, p_2, p_3};
                        sort(vector_to_add.begin(), vector_to_add.end());
                        size_t l_value = *max_element(vector_to_add.begin(), vector_to_add.end());
                        size_t flag_4 = 0;
                        size_t flag_5 = 0;
                        //cout << l_value << endl;
                        while (flag_4 == 0 and flag_5 == 0) {
                            if (FindNode(root_1, l_value) != zero_vector ){
                                // cout << "SEARCHING BINARY TREE 1" << endl;
                                vector<size_t> new_vector = FindNode(root_1, l_value);
                                vector<size_t> symm_diff_vec;
                                set_symmetric_difference(new_vector.begin(), new_vector.end(), vector_to_add.begin(), vector_to_add.end(), back_inserter(symm_diff_vec));
                                vector_to_add = symm_diff_vec;
                                symm_diff_vec.clear();
                                if (vector_to_add.size() > 0){
                                    l_value = *max_element(vector_to_add.begin(), vector_to_add.end());
                                } else {
                                    flag_5 = 1;
                                }
                                // sort(vector_to_add.begin(), vector_to_add.end());
                            } else if (FindNode(root_g_1, l_value) != zero_vector) {
                                // cout << "SEARCHING BINARY TREE 2" << endl;
                                vector<size_t> new_vector = FindNode(root_g_1, l_value);
                                vector<size_t> symm_diff_vec;
                                set_symmetric_difference(new_vector.begin(), new_vector.end(), vector_to_add.begin(), vector_to_add.end(), back_inserter(symm_diff_vec));
                                vector_to_add = symm_diff_vec;
                                symm_diff_vec.clear();
                                if (vector_to_add.size() > 0){
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
                            root_g_1 = insertNode(root_g_1, l_value, vector_to_add);
                            node_ctr_g_1 += 1;
                            // cout << "value of node_ctr_g_1 is " << node_ctr_g_1 << endl;
                            // At this stage (l_value, column_counter) is a persistent pair
                            double diam_one_simp = diam(idx_to_one_simp[l_value], point_matrix);
                            double diam_two_simp = diam(idx_to_two_simp[column_counter], point_matrix);
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



    }

    sort(barcode_bars.rbegin(), barcode_bars.rend());
    for (auto bar: barcode_bars){
        cout << bar[0] << ", " << bar[1] << endl;
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: " << double(duration.count() / 1000000) << " seconds" << endl;
    vector<double> test_vec;
    cout << test_vec.max_size() << endl;
    return 0;
}

