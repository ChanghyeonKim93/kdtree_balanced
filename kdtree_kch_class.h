#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>

inline double distance_euclidean(const double* _point_q, const double* _point_r, const int _ndim) {
	double temp = 0;
	for (int i = 0; i < _ndim; i++) {
		temp += (_point_q[i] - _point_r[i])*(_point_q[i] - _point_r[i]);
	}
	//std::cout<<"dist:"<<sqrt(temp)<<std::endl;
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
	BKDTree(const std::vector<std::vector<double>>& _points, const int& _binSize, const double& _dist_thres){
		this->npoints   = _points.size();
		this->ndim      = _points[0].size();
		this->binSize   = _binSize;
		this->dist_thres= _dist_thres;
		this->maxDepth = (int)round(log2(this->npoints / this->binSize));
		std::cout<<"\nB k-d tree handler initialization - npoints : "<<this->npoints<<", ndim : "<<this->ndim<<", binSize : "<<this->binSize<<", max depth : "<<this->maxDepth<<std::endl<<std::endl;
	}
	~BKDTree(){
		printf("\nB k-d tree handler is deleted !\n\n");
	}

	/// Member variables for class
private:
	int npoints;
	int ndim;
	int binSize;
	int maxDepth;
	double dist_thres;

	/// querying nearest neighbor
public:
	void search_NN_dfs(Node* node, const double* _point_q, const int _depth, int* _index);
	void search_NN_dfs_index(Node* _node, const std::vector<std::vector<double>>& _points_q_vec, std::vector<int> _index_vec);

	/// For debug
public:
	void print_points(const double* _point, const bool _isLeaf);
	void print_tree(const Node* node,  const int _depth);
	void print_leaf(const Node* node,  const int _depth);
	int  verify_tree(const Node* node, const int _depth);

	/// For kdtree build ( public functions )
public:

	Node* create_tree(double** _points);


	/// For kdtree build (intrinsic-only functions)
private:
	// void point_allocate(const std::vector<std::vector>>& _points_vec);
	void initialize_reference(double** _points, double** reference);
	PointNode* new_point_node(double* _point, int _index);
	Node* new_node(double* _point);
	void merge_sort(double** reference, double** temporary, const int low, const int high, const int _cur_dim);
	void insert_leaf_data(Node* node, double* _point, const int _depth, const int _index);
	int remove_duplicates(double** reference, const int _cur_dim);
	double super_key_compare(const double* _point_a, const double* _point_b, const int _cur_dim);
	Node* build_tree(double*** references, double** temp, const int start, const int end, const int _depth);
};


/* ====================================================================================
* ====================================================================================
* =================                  IMPLEMENTATION                  =================
* ====================================================================================
* ==================================================================================== */

PointNode* BKDTree::new_point_node(double* _point, int _index) {
	PointNode* point_node;

	if( (point_node = (PointNode*)malloc(sizeof(PointNode)))  == NULL) {
		printf("error allocating new PointNode ! \n");
		exit(1);
	}

	point_node->point = _point;
	point_node->index = _index;

	return point_node;
}


Node* BKDTree::new_node(double* _point) {
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


void BKDTree::initialize_reference(double** _points, double** reference) {
	for (int i = 0; i<this->npoints; i++) {
		reference[i] = _points[i];
	}
}


double BKDTree::super_key_compare(const double* _point_a, const double* _point_b, const int _cur_dim) {
	double diff = 0;
	for (int i = 0; i< this->ndim; i++) {
		int r = i + _cur_dim;
		r = (r < this->ndim) ? r : r - this->ndim;
		diff = _point_a[r] - _point_b[r];
		if (diff != 0) break; // diff�� 0�� �ƴϸ� ����
	}
	return diff;
}


void BKDTree::merge_sort(double** reference, double** temporary, const int low, const int high, const int _cur_dim) {
	int i, j, k;

	if (high > low) {
		const int mid = (high + low) / 2;
		merge_sort(reference, temporary, low,     mid,  _cur_dim);
		merge_sort(reference, temporary, mid + 1, high, _cur_dim);

		for (i = mid + 1; i > low; i--) {
			temporary[i - 1] = reference[i - 1];
		}

		for (j = mid; j < high; j++) {
			temporary[mid + high - j] = reference[j + 1];
		}
		for (k = low; k <= high; k++) {
			reference[k] = (BKDTree::super_key_compare(temporary[i], temporary[j], _cur_dim) < 0) ? temporary[i++] : temporary[j--];
		}
	}
}


int BKDTree::remove_duplicates(double** reference, const int _cur_dim) {
	int end = 0;
	for (int j = 1; j < this->npoints; j++) {
		double compare = BKDTree::super_key_compare(reference[j], reference[j - 1], _cur_dim);
		if (compare < 0) {
			printf("merge sort failure: super_key_compare(ref[%d], ref[%d]), current lvl=(%d), compare value = %lf\n", j, j - 1, _cur_dim, compare);
			exit(1);
		}
		else if (compare > 0) reference[++end] = reference[j];
	}
	return end;
}


Node* BKDTree::build_tree(double*** references, double** temp, const int start, const int end,  const int _depth) {

	Node* node = (Node*) malloc(sizeof(Node));
	int axis = _depth % this->ndim;
	//printf("current depth : %d, max depth : %d\n", _depth, _max_depth);

	if (_depth < this->maxDepth) { // max_depth ������ node�� �߰��Ѵ�.
		if (end == start) {
			node = BKDTree::new_node(references[0][end]);
		}
		else if (end == start + 1) {
			node        = BKDTree::new_node(references[0][start]);
			node->right = BKDTree::new_node(references[0][end]);
		}
		else if (end == start + 2) {
			node        = BKDTree::new_node(references[0][start + 1]);
			node->left  = BKDTree::new_node(references[0][start]);
			node->right = BKDTree::new_node(references[0][end]);
		}
		else if (end > start + 2) {
			const int median = (start + end) / 2;
			node        = BKDTree::new_node(references[0][median]);

			for (int i = start; i <= end; i++) {
				temp[i] = references[0][i];
			}

			int lower, upper, lowerSave, upperSave;

			for (int i = 1; i < this->ndim; i++) {
				lower = start - 1;
				upper = median;
				for (int j = start; j <= end; j++) {
					double compare = BKDTree::super_key_compare(references[i][j], node->point, axis);
					if (compare < 0) {
						references[i - 1][++lower] = references[i][j];
					}
					else if (compare > 0) {
						references[i - 1][++upper] = references[i][j];
					}
				}
				lowerSave = lower;
				upperSave = upper;
			}
			for (int i = start; i <= end; i++) {
				references[this->ndim - 1][i] = temp[i];
			}
			node->left  = BKDTree::build_tree(references, temp, start, lower, _depth + 1);
			node->right = BKDTree::build_tree(references, temp, median + 1, upper, _depth + 1);
		}
	}
	else { // In case of the leaf node.
		double* temp = (double*)malloc(sizeof(double)*this->ndim);
		for (int i = 0; i<this->ndim; i++) temp[i] = 0;
		node = BKDTree::new_node(temp);
		node->left = node->right = NULL;
		node->isLeaf = true;
		//printf("Leaf exit\n");
	}

	return node; // tree root �� �� ���̴�.
}


void BKDTree::insert_leaf_data(Node* node, double* _point, const int _depth, const int _index) {// in case of reaching to a leaf, make the leaf node as a leaf node containing several points.

	PointNode* point_node = (PointNode*) malloc(sizeof(PointNode));

	if (node->isLeaf) { // In case of leaf node, push_back the new point into the node.
		node->leafData.push_back(_point);
		point_node = BKDTree::new_point_node(_point, _index);

		node->pointNodes.push_back(point_node);
		//std::cout << "Point node - index : " << point_node->index << std::endl;
		node->numOfPoints++ ;
		//std::cout << "reach the leaf, " << node << " , " << node->numOfPoints << "," << std::endl;
		// memory problem? Linux shows no problems.
	}

	else { // traveling nodes until reaching the leaf node.
		free(point_node); // point_node is not used. Therefore, point_node is deleted.
		if (_point[_depth] <= node->point[_depth]) { // Go to left child
			BKDTree::insert_leaf_data(node->left,  _point, (_depth + 1) % this->ndim, _index);
		}
		else { // Go to right child
			BKDTree::insert_leaf_data(node->right, _point, (_depth + 1) % this->ndim, _index);
		}
	}
}


Node* BKDTree::create_tree(double** _points) {

	double*** references = (double***)malloc(this->ndim*sizeof(double**));
	double** temp = (double**)malloc(this->npoints*sizeof(double));
	for (int i = 0; i < this->ndim; i++) {   // dimension ��ŭ merge �ؼ� references�� �����Ѵ�.
		references[i] = (double**)malloc(this->npoints*sizeof(double*));
		BKDTree::initialize_reference(_points, references[i]);
		BKDTree::merge_sort(references[i], temp, 0, this->npoints - 1, i);
	}

	int *end = (int*)malloc(this->ndim*sizeof(int));
	for (int i = 0; i< this->ndim; i++) { // �ߺ� �� �� �� �������ش�.
		end[i] = BKDTree::remove_duplicates(references[i], i);
	}

	for (int i = 0; i < this->ndim - 1; i++) {
		for (int j = i + 1; j < this->ndim; j++) {
			if (end[i] != end[j]) {
				printf("reference removal error in create_tree\n");
				exit(1);
			}
		}
	} // ���������� ������ ����.

	Node* root = BKDTree::build_tree(references, temp, 0, end[0], 0); // ���Ⱑ �����ΰ� ?

	for (int j = 0; j < this->npoints; j++) {
		BKDTree::insert_leaf_data(root, _points[j], 0, j);
		//std::cout<<j<<std::endl;
	}

	for (int i = 0; i < this->ndim; i++) {
		free(references[i]);
	}
	free(references);
	free(temp);

	return root;
}


void BKDTree::search_NN_dfs(Node* node, const double* _point_q, const int _depth, int* _index) {
	PointNode* currPointNode;

	if (node->isLeaf) { // In case of leaf node, push_back the new point into the node.
		double minDist = 1000000000;
		*_index = -1;

		for (int k = 0; k<node->numOfPoints; k++) {
			//std::cout<<node->pointNodes[k]->index<<std::endl;
			double currDist = distance_euclidean(_point_q, node->pointNodes[k]->point, this->ndim);
			//double currDist = distance_manhattan(_point_q, node->pointNodes[k]->point, this->ndim);
			if (currDist < minDist) {
				minDist = currDist;
				currPointNode = node->pointNodes[k];
				*_index = node->pointNodes[k]->index;

				if(minDist<=this->dist_thres) return;
			}
			//std::cout<<*_index<<std::endl;
		}
	}
	else { // traveling nodes until reaching the leaf node.
		if (_point_q[_depth] <= node->point[_depth]) { // Go to left child
			BKDTree::search_NN_dfs(node->left,  _point_q, (_depth + 1) % this->ndim,_index);
		}
		else { // Go to right child
			BKDTree::search_NN_dfs(node->right, _point_q, (_depth + 1) % this->ndim,_index);
		}
	}
	return ;
}

void BKDTree::search_NN_dfs_index(Node* _node, const std::vector<std::vector<double>>& _points_q_vec, std::vector<int> _index_vec){

	int* indexTemp = new int;
	const double** _points_q;
	int npoints_q = _points_q_vec.size();
	int ndim_q    = _points_q_vec[0].size();
	std::vector<int> index_vec;
	index_vec.reserve(npoints_q);

	if(ndim_q != this->ndim){
		printf("Error - query point dimension does not match reference points ! \n");
		exit(1);
	}

	double* _point_q = (double*)malloc(sizeof(double)*ndim_q);
	for(int i = 0; i < npoints_q; i++){
		for(int j = 0; j < ndim_q; j++) _point_q[j] = _points_q_vec[i][j];
		BKDTree::search_NN_dfs(_node, _point_q, 0, indexTemp);
		index_vec.push_back(*indexTemp);
	}


	delete indexTemp;
	_index_vec.swap(index_vec);

	return; // void return
}




int BKDTree::verify_tree(const Node* node, const int _depth)
{
	int count = 1;
	if (node->point == NULL) {
		printf("point is null\n");
		exit(1);
	}
	// The partition cycles as x, y, z, w...
	int axis = _depth % this->ndim;

	if (node->left != NULL) {
		if (node->left->point[axis] > node->point[axis]) {
			//printf("child is > node!\n");
			//exit(1);
		}
		if (BKDTree::super_key_compare(node->left->point, node->point, axis) >= 0) {
			//printf("child is >= node!\n");
			//exit(1);
		}
		count += BKDTree::verify_tree(node->left, _depth + 1);
	}
	if (node->right != NULL) {
		if (node->right->point[axis] < node->point[axis]) {
			//printf("child is < node!\n");
			//exit(1);
		}
		if (BKDTree::super_key_compare(node->right->point, node->point, axis) <= 0) {
			//printf("child is <= node!\n");
			//exit(1);
		}
		count += BKDTree::verify_tree(node->right, _depth + 1);
	}

	return count; // the number of nodes.
}



void BKDTree::print_points(const double* _point, const bool _isLeaf) {
	printf("[ %0.0f,", _point[0]);
	for (int i = 1; i<this->ndim - 1; i++) printf("%0.0f,", _point[i]);
	printf("%0.0f ]", _point[this->ndim - 1]);
	if (_isLeaf == true) {
		printf("< LEAF >");
	}
	else {
		printf("<nd>");
	}
}

void BKDTree::print_tree(const Node* node, const int _depth) {
	if (node != NULL) {
		BKDTree::print_tree(node->right, _depth + 1);

		for (int i = 0; i<_depth; i++) {
			printf("             ");
		}
		BKDTree::print_points(node->point, node->isLeaf);
		printf("\n");
		BKDTree::print_tree(node->left, _depth + 1);
	}
}

void BKDTree::print_leaf(const Node* node, const int _depth) {
	if (node != NULL) {
		BKDTree::print_leaf(node->right, _depth + 1);
		if (node->isLeaf) {
			for (int i = 0; i < node->numOfPoints; i++) {
				std::cout << "coordinate : [";
				for (int j = 0; j< this->ndim; j++) {
					std::cout << node->pointNodes[i]->point[j] << ",";
				}
				std::cout << " ], index : " << node->pointNodes[i]->index << std::endl;
			}
		}
		BKDTree::print_leaf(node->left, _depth + 1);
	}
}
