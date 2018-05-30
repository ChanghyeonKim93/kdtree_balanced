#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>

inline double distance_euclidean(const double* _point_q, const double* _point_r, const int _ndim) {
	double temp = 0;
	for (int i = 0; i < _ndim; i++) {
		temp += (_point_q[i] - _point_r[i])*(_point_q[i] - _point_r[i]);
	}
	return sqrt(temp);
}

inline double distance_manhattan(const double* _point_q, const double* _point_r, const int _ndim) {
	double temp = 0;
	for (int i = 0; i < _ndim; i++) {
		temp += abs(_point_q[i] - _point_r[i]);
	}
	return temp;
}

typedef struct _PointNode{
  double* point;
  int index;
} PointNode;

typedef struct _Node{
  double* point;
  _Node* left;
  _Node* right;

  std::vector<double*>    leafData;
  std::vector<PointNode*> pointNodes;
  int numOfPoints;
  bool isLeaf;
} Node;

class BKDTree{

/// Constructor & destructor
public:
  BKDTree(){}
  BKDTree(const std::vector<std::vector<double>>& _points, const int& _binSize){
    this->npoints   = _points.size();
    this->ndim      = _points[0].size();
    this->binSize   = _binSize;
    this->max_depth = (int)round(log2(this->npoints / this->binSize));
  }
  ~BKDTree(){
    printf("B k-d tree is deleted !\n");
  };

/// Member variables for class
private:
  int npoints;
  int ndim;
  int binSize;
  int max_depth;


/// querying nearest neighbor
public:
  void search_NN_dfs(Node* node, const double* _point_q, const int _depth, const int _ndim, int* _index);


/// For debug
public:
  void print_points(const double* _point, const int _ndim, const bool _isLeaf);
  void print_tree(const Node* node, const int _ndim, const int _depth);
  void print_leaf(const Node* node, const int _ndim, const int _depth);
  int  verify_tree(const Node* node, const int _ndim, const int _depth);


/// For kdtree build ( public functions )
public:
  Node* create_tree(double** _points, const int _npoints, const int _ndim, const int _max_depth);


/// For kdtree build (intrinsic-only functions)
private:
  void initialize_reference(double** _points, double** reference, const int _npoints);
  PointNode* new_point_node(double* _point, int _index);
  Node* new_node(double* _point);
  void merge_sort(double** reference, double** temporary, const int low, const int high, const int _cur_dim, const int _ndim);
  void insert_leaf_data(Node* node, double* _point, const int _depth, const int _ndim, const int _index);
  int remove_duplicates(double** reference, const int _npoints, const int _cur_dim, const int _ndim);
  double super_key_compare(const double* _point_a, const double* _point_b, const int _cur_dim, const int _ndim);
  Node* build_tree(double*** references, double** temp, const int start, const int end, const int _ndim, const int _depth, const int _max_depth);
};



/* ====================================================================================
 * ====================================================================================
 * =================                  IMPLEMENTATION                  =================
 * ====================================================================================
 * ==================================================================================== */

 BKDTree::PointNode* new_point_node(double* _point, int _index) {
 	PointNode* point_node;

 	if( (point_node = (PointNode*)malloc(sizeof(PointNode)))  == NULL) {
 		printf("error allocating new PointNode ! \n");
 		exit(1);
 	}

 	point_node->point = _point;
 	point_node->index = _index;

 	return point_node;
 }


 BKDTree::Node* new_node(double* _point) {
 	Node* node;

 	if ((node = (Node*)malloc(sizeof(Node))) == NULL) {
 		printf("error allocating new Node ! \n");
 		exit(1);
 	}
 	// initialize the data w/o leafData
 	node->point = _point;
 	node->left  = node->right = NULL;

 	node->leafData.resize(0,0);
 	node->pointNodes.resize(0,0);
 	node->numOfPoints = 0;
 	node->isLeaf      = false;

 	return node;
 }


 BKDTree::void initialize_reference(double** _points, double** reference) {
 	for (int i = 0; i<this->npoints; i++) {
 		reference[i] = _points[i];
 	}
 }

 BKDTree::double super_key_compare(const double* _point_a, const double* _point_b, const int _cur_dim) {
 	double diff = 0;
 	for (int i = 0; i< this->ndim; i++) {
 		int r = i + _cur_dim;
 		r = (r < this->ndim) ? r : r - this->ndim;
 		diff = _point_a[r] - _point_b[r];
 		if (diff != 0) break; // diff�� 0�� �ƴϸ� ����
 	}
 	return diff;
 }
