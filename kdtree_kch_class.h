#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>
//typedef std::vector<double> Point;

inline double distance_euclidean(const double* _point_q, const double* _point_r, const int _nDims) {
	double temp = 0;
	for (int i = 0; i < _nDims; i++) temp += (_point_q[i] - _point_r[i])*(_point_q[i] - _point_r[i]);
	//std::cout<<"dist:"<<sqrt(temp)<<std::endl;
	return sqrt(temp);
}

inline double distance_manhattan(const double* _point_q, const double* _point_r, const int _nDims) {
	double temp = 0;
	for (int i = 0; i < _nDims; i++) temp += fabs(_point_q[i] - _point_r[i]);
	return temp;
}

typedef struct _PointNode {
	double* point;
	int index;
} PointNode;

typedef struct _Node {
	double* point;
	_Node* left;
	_Node* right;

	std::vector<double*>    leafData;
	std::vector<PointNode*> pointNodes;
	int numOfPoints;
	bool isLeaf;
} Node;

class BKDTree {

    /// Constructor & destructor
public:
	BKDTree(){}
	BKDTree(const std::vector<std::vector<double>>& _points_vec, const int& _binSize, const double& _distThres);
	~BKDTree();

	/// Member variables for class
public:
	int nPoints;
	int nDims;
	int binSize;
	int maxDepth;
	double distThres;

	Node* treeRootNode;
	std::vector<Node*> nodesPtrs;

	/// querying nearest neighbor
public:
	void search_NN_dfs(Node* node, const double* _point_q, const int _depth, int* _index);
	void search_NN_dfs_index(Node* _node, const std::vector<std::vector<double>>& _points_q_vec, std::vector<int>& _index_vec);

	void kdtree_nearest_neighbor_approximate(const std::vector<std::vector<double>>& _points_q_vec,  std::vector<int>& _index_vec);

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
	Node* build_tree_recursively(double*** references, double** temp, const int start, const int end, const int _depth);
	void delete_malloc();
};


/**
 * ----------------------------------------------------------------------------
 * 
 * 
 *                                IMPLEMENTATION
 * 
 * 
 * ----------------------------------------------------------------------------
 */

/**
 * Method name            : BKDTree()
 * Function of the method : Creates a BKDTree ( depth-first-search version ) with the provided vector data.
 * 
 * @param _points_vec     : a std::vector<std::vector<double>> containing the point data.
 *                          the number of points and the dimensionality are inferred by the data.
 * @param _binSize        : the size of the leaf bin. recommended range is from 5 to 20.
 * @param _distThres      : distance threshold value for rejecting matches with high potential of long matching.
 */
BKDTree::BKDTree(const std::vector<std::vector<double>>& _points_vec, const int& _binSize, const double& _distThres){
		this->nPoints   = _points_vec.size();
		this->nDims     = _points_vec[0].size();
		this->binSize   = _binSize;
		this->distThres = _distThres;
		this->maxDepth  = (int)ceil(log2(  ceil((double)this->nPoints / (double)this->binSize) ) );
		std::cout<<"\nB k-d tree handler initialization - nPoints : "<<this->nPoints<<", nDims : "<<this->nDims<<", binSize : "<<this->binSize<<", max depth : "<<this->maxDepth<<std::endl<<std::endl;
		
		double** points = (double**)malloc(this->nPoints*sizeof(double*)); // change the data type.
		for (int i = 0; i < this->nPoints; i++)
		{
			points[i] = (double*)malloc(this->nDims*sizeof(double));
			for (int j = 0; j < this->nDims; j++) points[i][j] = _points_vec[i][j];
		}

		this->treeRootNode = this->create_tree(points); // allocate the treeRootNode address

		this->nodesPtrs.reserve(0);
		
		// free the memory.
		for(int i= 0; i<this->nPoints; i++) free(points[i]);
		free(points);
}

BKDTree::~BKDTree(){
	BKDTree::delete_malloc();
	printf("\nB k-d tree handler is deleted !\n\n");
}

PointNode* BKDTree::new_point_node(double* _point, int _index) {
	PointNode* point_node;

	if( (point_node = (PointNode*)malloc(sizeof(PointNode)))  == NULL) {
		printf("error allocating new PointNode ! \n");
		exit(1);
	}

	point_node->point = (double*)malloc(sizeof(double)*this->nDims);
	for(int i = 0; i < this->nDims; i++) point_node->point[i] = _point[i];
	point_node->index = _index;

	return point_node;
}



void BKDTree::initialize_reference(double** _points, double** reference) {
	for (int i = 0; i<this->nPoints; i++) reference[i] = _points[i];
}


double BKDTree::super_key_compare(const double* _point_a, const double* _point_b, const int _cur_dim) {
	double diff = 0;
	for (int i = 0; i < this->nDims; i++) 
	{
		int r = i + _cur_dim;
		r = (r < this->nDims) ? r : r - this->nDims; // circular
		diff = _point_a[r] - _point_b[r];
		if (diff != 0) break; // 만약, 두개가 같지 않으면 ( 즉, 다르면 ) 그냥 패스 ~
	}
	return diff; // 만약, 모든 차원의 숫자가 같으면 0이 나오고 / a가 크면 (+) / b가 크면 (-).
}


void BKDTree::merge_sort(double** reference, double** temporary, const int low, const int high, const int _cur_dim) {
	int i, j, k;

	if (high > low) 
	{
		const int mid = (high + low) / 2;
		merge_sort(reference, temporary, low,     mid,  _cur_dim);
		merge_sort(reference, temporary, mid + 1, high, _cur_dim);

		for (i = mid + 1; i > low; i--)
		{
			temporary[i - 1] = reference[i - 1];
		}
		for (j = mid; j < high; j++)
		{
			temporary[mid + high - j] = reference[j + 1];
		}
		for (k = low; k <= high; k++)
		{
			reference[k] = (BKDTree::super_key_compare(temporary[i], temporary[j], _cur_dim) < 0) ? temporary[i++] : temporary[j--];
		}
	}
}


int BKDTree::remove_duplicates(double** reference, const int _cur_dim) {
	int end = 0;
	for (int j = 1; j < this->nPoints; j++) 
	{
		double compare = BKDTree::super_key_compare(reference[j], reference[j - 1], _cur_dim);
		if (compare < 0) 
		{
			printf("merge sort failure: super_key_compare(ref[%d], ref[%d]), current lvl=(%d), compare value = %lf\n", j, j - 1, _cur_dim, compare);
			exit(1);
		}
		else if (compare > 0) reference[++end] = reference[j];
	}
	return end;
}


Node* BKDTree::build_tree_recursively(double*** references, double** temp, const int start, const int end,  const int _depth) {
	// references : already sorted ! 
	// initial values : start ( 0 ), end ( npoints - 1 )
	Node* node = (Node*) malloc(sizeof(Node)); // root node
	int axis = (_depth < this->nDims) ? _depth : _depth - this->nDims; 
	
	//printf("current depth : %d, max depth : %d, axis : %d\n", _depth, this->maxDepth,axis);

	if(_depth < this->maxDepth) // max_depth ������ node�� �߰��Ѵ�.
	{
		if (end == start) // 1개 남았을 때
		{
			node        = BKDTree::new_node(references[0][end]);
			this->nodesPtrs.push_back(node);
		}
		else if (end == start + 1) // 2개 남았을 때
		{
			node        = BKDTree::new_node(references[0][start]);
			node->right = BKDTree::new_node(references[0][end]);

			this->nodesPtrs.push_back(node);
			this->nodesPtrs.push_back(node->right);
		}
		else if (end == start + 2) // 3개 남았을 때
		{
			node        = BKDTree::new_node(references[0][start + 1]);
			node->left  = BKDTree::new_node(references[0][start]);
			node->right = BKDTree::new_node(references[0][end]);

			this->nodesPtrs.push_back(node);
			this->nodesPtrs.push_back(node->right);
			this->nodesPtrs.push_back(node->left);
		}
		else if (end >  start + 2) // 4개 이상 남았을 때
		{
			const int median = (start + end) / 2; // 내림된다. ex) (0+5)/2 = 2.5 (x) / 2 (o) , (0+4)/2 = 2
			node        = BKDTree::new_node(references[0][median]);

			for (int i = start; i <= end; i++)
			{
				temp[i] = references[0][i];
			}

			int lower, upper, lowerSave, upperSave;

			for (int i = 1; i < this->nDims; i++)
			{
				lower = start - 1;
				upper = median;

				for (int j = start; j <= end; j++)
				{
					double compare = BKDTree::super_key_compare(references[i][j], node->point, axis);
					if (compare < 0) // 만약, references[i][j]가 현재 중심점인 node->point 보다 모든 차원에서 작다면, (한 차원에서만 작아도)
					{
						references[i - 1][++lower] = references[i][j]; // 
					}
					else if (compare > 0) // 만약, references[i][j]가 현재 중심점인 node->point 보다 모든 차원에서 크다면, (한 차원에서만 커도)
					{
						references[i - 1][++upper] = references[i][j];
					}
				}
				lowerSave = lower;
				upperSave = upper;
			}

			for (int i = start; i <= end; i++)
			{
				references[this->nDims - 1][i] = temp[i];
			}
			node->left  = BKDTree::build_tree_recursively(references, temp, start, lower, _depth + 1);
			node->right = BKDTree::build_tree_recursively(references, temp, median + 1, upper, _depth + 1);

			this->nodesPtrs.push_back(node);
			//this->nodesPtrs.push_back(node->left);
			//this->nodesPtrs.push_back(node->right);
		}
	}
	else // In case of the leaf node.
	{ 
		double* temp = (double*)malloc(sizeof(double)*this->nDims);
		for (int i = 0; i < this->nDims; i++) temp[i] = -3;
		node = BKDTree::new_node(temp);
		this->nodesPtrs.push_back(node);
		
		node->left   = NULL;
		node->right  = NULL;
		node->isLeaf = true;
		//printf("Leaf exit\n");
	}

	return node; // tree root �� �� ���̴�.
}




Node* BKDTree::create_tree(double** _points) {
	double*** references = (double***)malloc(this->nDims*sizeof(double**));
	double**  temp       = (double**)malloc(this->nPoints*sizeof(double*));
	for (int i = 0; i < this->nDims; i++) // dimension
	{   
		references[i] = (double**)malloc(this->nPoints*sizeof(double*));
		BKDTree::initialize_reference(_points, references[i]);
		BKDTree::merge_sort(references[i], temp, 0, this->nPoints - 1, i);
	}

	int *end = (int*)malloc(this->nDims*sizeof(int));
	for (int i = 0; i< this->nDims; i++) // 트리 구성 중, 완전 동일한 점 삭제.
	{ 
		end[i] = BKDTree::remove_duplicates(references[i], i);
	}

	for (int i = 0; i < this->nDims - 1; i++) // 만약, 병합정렬 된 모든 차원에서의 end가 같지 않으면, 제대로 removal이 이루어지지 않은 것이므로 에러발생.
	{
		for (int j = i + 1; j < this->nDims; j++) 
		{
			if (end[i] != end[j]) 
			{
				printf("reference removal error in create_tree\n");
				exit(1);
			}
		}
	} 

	Node* root = BKDTree::build_tree_recursively(references, temp, 0, end[0], 0); //


	for (int i = 0; i < this->nPoints; i++) 
	{
		BKDTree::insert_leaf_data(root, _points[i], 0, i);
	}
	
	for (int i = 0; i < this->nDims; i++) 
	{
		free(references[i]);
	}
	free(references);
	free(temp);

	return root;
}


void BKDTree::search_NN_dfs(Node* node, const double* _point_q, const int _depth, int* _index) {
	PointNode* currPointNode;

	if (node->isLeaf) // In case of leaf node, push_back the new point into the node.
	{ 
		double minDist = 1000000000;
		*_index = -1;

		for (int k = 0; k<node->numOfPoints; k++) 
		{
			//std::cout<<node->pointNodes[k]->index<<std::endl;
			double currDist = distance_euclidean(_point_q, node->pointNodes[k]->point, this->nDims);
			//double currDist = distance_manhattan(_point_q, node->pointNodes[k]->point, this->nDims);
			if (currDist < minDist) 
			{
				minDist = currDist;
				currPointNode = node->pointNodes[k];
				*_index = node->pointNodes[k]->index;

				if(minDist<=this->distThres) return;
			}
			//std::cout<<*_index<<std::endl;
		}
	}
	else // traveling nodes until reaching the leaf node.
	{ 
		if (_point_q[_depth] <= node->point[_depth]) // Go to left child
		{ 
			BKDTree::search_NN_dfs(node->left,  _point_q, (_depth + 1) % this->nDims,_index);
		}
		else  // Go to right child
		{	
			BKDTree::search_NN_dfs(node->right, _point_q, (_depth + 1) % this->nDims,_index);
		}
	}
	return ;
}

void BKDTree::search_NN_dfs_index(Node* _node, const std::vector<std::vector<double>>& _points_q_vec, std::vector<int>& _index_vec){

	int* indexTemp = new int;
	const double** _points_q;
	int nPoints_q = _points_q_vec.size();
	int nDims_q    = _points_q_vec[0].size();
	std::vector<int> index_vec;
	index_vec.reserve(nPoints_q);

	if(nDims_q != this->nDims)
	{
		printf("Error - query point dimension does not match reference points ! \n");
		exit(1);
	}

	double* _point_q = (double*)malloc(sizeof(double)*this->nDims);
	for(int i = 0; i < nPoints_q; i++)
	{
		for(int j = 0; j < this->nDims; j++) _point_q[j] = _points_q_vec[i][j];
		BKDTree::search_NN_dfs(_node, _point_q, 0, indexTemp);
		index_vec.push_back(*indexTemp);
		_index_vec.push_back(*indexTemp);
		// std::cout<<index_vec[i]<<std::endl;
	}


	delete indexTemp;
	//_index_vec.swap(index_vec);

	return; // void return
}


void BKDTree::kdtree_nearest_neighbor_approximate(const std::vector<std::vector<double>>& _points_q_vec,  std::vector<int>& _index_vec){
	this->search_NN_dfs_index(this->treeRootNode, _points_q_vec, _index_vec);
}

void BKDTree::delete_malloc()
{
	// delete all malloc s in the nodes. Following, delete nodes.
	int numNodes = this->nodesPtrs.size();
	printf("num nodes : %d\n",numNodes);
	
	for(int i=0; i < numNodes; i++)
	{
		for(int j=0;j<this->nodesPtrs[i]->pointNodes.size(); j++) // delete pointNodes.
		{
			free(this->nodesPtrs[i]->pointNodes[j]->point);
			free(this->nodesPtrs[i]->pointNodes[j]);
		}
		free(this->nodesPtrs[i]->point);
		free(this->nodesPtrs[i]); // finally, delete nodes.
	}
	
	printf("All memories of nodes and pointNodes are returned !\n");
}


int BKDTree::verify_tree(const Node* node, const int _depth)
{
	int count = 1;
	if (node->point == NULL) {
		printf("point is null\n");
		exit(1);
	}
	// The partition cycles as x, y, z, w...
	int axis = _depth % this->nDims;

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
	printf("[ %0.0lf,", _point[0]);
	for (int i = 1; i<this->nDims - 1; i++) printf("%0.0lf,", _point[i]);
	printf("%0.0lf ]", _point[this->nDims - 1]);
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
				for (int j = 0; j< this->nDims; j++) {
					std::cout << node->pointNodes[i]->point[j] << ",";
				}
				std::cout << " ], index : " << node->pointNodes[i]->index << std::endl;
			}
		}
		BKDTree::print_leaf(node->left, _depth + 1);
	}
}
