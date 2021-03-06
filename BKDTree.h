#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

inline double distance_euclidean(const std::vector<double>& _point_q, const std::vector<double>& _point_r, const int nDims, const bool& useGradient) {
	double temp = 0.0;
	double dSpatial = 0.0;

	if( nDims == 4 && useGradient == true)
	{
		temp = (_point_q[2]*_point_r[2] + _point_q[3]*_point_r[3] - 1.0);
		dSpatial += (_point_q[0] - _point_r[0])*(_point_q[0] - _point_r[0]) + (_point_q[1] - _point_r[1])*(_point_q[1] - _point_r[1]);
		return dSpatial + 0.175*temp*temp;
	} else {
		for (int i = 0; i < nDims; i++) temp += (_point_q[i] - _point_r[i])*(_point_q[i] - _point_r[i]);
		return sqrt(temp);
	}
	//printf("dist:%0.3lf\n",sqrt(temp));
}


inline double distance_manhattan(const std::vector<double>& _point_q, const std::vector<double>& _point_r, const int _nDims) {
	double temp = 0;
	for (int i = 0; i < _nDims; i++) temp += fabs(_point_q[i] - _point_r[i]);
	return temp;
}

double** double_array_allocator_2d(int npoints, int ndim)
{
	double** ptr = new double*[npoints];
	for(int i = 0; i < npoints; i++) ptr[i] = new double[ndim];
	return ptr;
	printf("this->points : allocated. Memory size : %d\n", sizeof(double)*npoints*ndim);
}

void double_array_delete_2d(double** ptr, int npoints)
{
	for(int i = 0; i < npoints; i++) delete[] ptr[i];
	delete[] ptr;
	//printf("this->points : deleted.\n");
}


void double_array_delete_3d(double*** ptr, int ndims, int npoints)
{
	printf("ndims : %d, npoints : %d\n", ndims, npoints);

	for(int i = 0; i < ndims; i++)
	{
		for(int j = 0; j < npoints; j++)
		{
			delete[] ptr[i][j];
		}
		delete[] ptr[i];
	}
	delete[] ptr;
	printf("[INFO] DELETION - this->references : deleted.\n");
}


typedef struct _PointNode {
	std::vector<double> point;
	int index;
} PointNode;

typedef struct _Node {
	_Node* left;
	_Node* right;
	std::vector<double> refPoint;

	// leaf node only.
	std::vector<PointNode*> pointNodes;
	int numOfPoints;
	bool isLeaf;
} Node;

class BKDTree {
/// Constructor & destructor
public:
	BKDTree(){};
	BKDTree(const std::vector<std::vector<double>>& _points_vec, const int& _binSize, const double& _distThres);
	~BKDTree();

/// Member variables for class
public:
	int nPoints;
	int nDims;
	int binSize;
	int maxDepth;
	int nNodes;
	double distThres;
	double** points;

	Node* treeRootNode;
	std::vector<Node*> nodesPtrs; // carry all addresses of nodes.

/// querying nearest neighbor
public:
	void search_NN_depth_first_search(Node* node, const std::vector<double>& _point_q, const int _depth, int* _index);
	void search_NN_depth_first_search_multiple(Node* _node, const std::vector<std::vector<double>>& _points_q_vec, std::vector<int>& _index_vec);
	void kdtree_nearest_neighbor_approximate(const std::vector<std::vector<double>>& _points_q_vec,  std::vector<int>& _index_vec);
	
/// For kdtree build ( public functions )
public:
	Node* create_tree(double** _points);

public:
	void print_points(const std::vector<double>& _point, const bool _isLeaf);
	void print_tree(const Node* node, const int _depth);
	void print_leaf(const Node* node, const int _depth);


/// For kdtree build (intrinsic-only functions)
private:
	void initialize_reference(double** _points, double** reference);
	PointNode* new_point_node(double* _point, int _index);
	Node* new_node(double* _point);
	Node* new_node_leaf(const double& dummyNum);
	void merge_sort(double** reference, double** temporary, const long low, const long high, const int _cur_dim);
	void insert_leaf_data(Node* node, double* _point, const int _depth, const int _index);
	int remove_duplicates(double** reference, const int _cur_dim);
	double super_key_compare(const double* _point_a, const double* _point_b, const int _cur_dim);
	double super_key_compare(const double* _point_a, double* _point_b, const int _cur_dim);

	Node* build_tree_recursively(double*** references, double** temp, const int start, const int end, const int _depth);
	void delete_node_recursively(Node* node);
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

		this->nNodes    = 0;

		printf("\n[INFO] (IN CLASS) Balanced k-d tree handler initialization - nPoints:%d, ",this->nPoints);
		printf("nDims:%d, ",       this->nDims);
		printf("binSize:%d, ",     this->binSize);
		printf("max Depth:%d\n\n", this->maxDepth);
		printf("[INFO] Expected memory consumption (dynamic allocations) : %0.0lf [KBytes]\n", 0.001*( (double)sizeof(double)*(double)this->nDims*(double)this->nPoints + (double)sizeof(double)*(double)this->nDims*(double)this->nPoints*(double)this->nDims ));

		
		this->points = NULL;
		this->points = double_array_allocator_2d( this->nPoints, this->nDims ); // change the data type.

		for (int i = 0; i < this->nPoints; i++)
			for (int j = 0; j < this->nDims; j++) this->points[i][j] = _points_vec[i][j];
		
		this->treeRootNode = this->create_tree(this->points); // Create a tree and allocate the treeRootNode address
}

BKDTree::~BKDTree(){
	double_array_delete_2d(this->points, this->nPoints);
	printf("[INFO] DELETION - Total # of existing nodes: %d\n", this->nNodes);
	delete_node_recursively(this->treeRootNode);
	printf("[INFO] DELETION - Balanced k-d tree handler is deleted !\n\n");
}

Node* BKDTree::new_node(double* _point) { // send the node address 
	Node* node = NULL;
	if ( (node = new Node() ) == NULL) 
	{
		printf("error allocating new Node ! \n");
		exit(1);
	}
	// initialize the data
	node->refPoint.reserve(this->nDims);
	node->left     = NULL;
	node->right    = NULL;

	// allocating the new point by deep copying.
	for(int i = 0; i<this->nDims; i++) node->refPoint.push_back(_point[i]);

	//node->pointNodes.resize(0,0);
	node->numOfPoints = 0;
	node->isLeaf      = false;
 
	return node;
}

PointNode* BKDTree::new_point_node(double* _point, int _index) {
	PointNode* pointNode = NULL;
	if( (pointNode = new PointNode() ) == NULL)
	{
		printf("error allocating new PointNode ! \n");
		exit(1);
	}

	// initialize the data
	for(int i = 0; i < this->nDims; i++) pointNode->point.push_back(_point[i]);
	
	pointNode->index = _index;

	return pointNode;
}

Node* BKDTree::new_node_leaf(const double& inputNum) { // send the node address 
	Node* node = NULL;
	if ( (node = new Node() ) == NULL)
	{
		printf("error allocating new Node ! \n");
		exit(1);
	}
	// initialize the data
	node->refPoint.reserve(this->nDims);
	node->left     = NULL;
	node->right    = NULL;

	// allocating the new point by deep copying.
	for(int i = 0; i < this->nDims; i++) node->refPoint.push_back(inputNum);

	// node->pointNodes.resize(0,0);
	node->numOfPoints = 0;
	node->isLeaf      = true;

	return node;
}


void BKDTree::initialize_reference(double** _points, double** reference) {
	for (int j = 0; j < this->nPoints; j++)
		for(int k = 0; k < this->nDims; k++) reference[j][k] = _points[j][k];			
}


double BKDTree::super_key_compare(const double* _point_a, const double* _point_b, const int _cur_dim) {
	double diff = 0;
	for (int i = 0; i < this->nDims; i++) 
	{
		int r = i + _cur_dim;
		r = (r < this->nDims) ? r : r - this->nDims; // circular
		diff = _point_a[r] - _point_b[r];
		if (diff != 0) break; // 만약, ?�개�?같�? ?�으�?( �? ?�르�?) 그냥 ?�스 ~
	}
	return diff; // 만약, 모든 차원???�자�?같으�?0???�오�?/ a�??�면 (+) / b�??�면 (-).
}


double BKDTree::super_key_compare(const double* _point_a, double* _point_b, const int _cur_dim) {
	double diff = 0;
	for (int i = 0; i < this->nDims; i++) 
	{
		int r = i + _cur_dim;
		r = (r < this->nDims) ? r : r - this->nDims; // circular
		diff = _point_a[r] - _point_b[r];
		if (diff != 0) break; // 만약, ?�개�?같�? ?�으�?( �? ?�르�?) 그냥 ?�스 ~
	}
	return diff; // 만약, 모든 차원???�자�?같으�?0???�오�?/ a�??�면 (+) / b�??�면 (-).
}

void BKDTree::merge_sort(double** reference, double** temporary, const long low, const long high, const int _cur_dim) {
	long i, j, k;
	double* pointTemp = NULL;
	if (high > low)
	{
		const long mid = low + ( (high - low) >> 1 );
		merge_sort(reference, temporary, low,     mid,  _cur_dim);
		merge_sort(reference, temporary, mid + 1, high, _cur_dim);

		for (i = mid + 1; i > low; i--)
		{
			for(int m=0;m < this->nDims; m++) temporary[i - 1][m] = reference[i - 1][m];
			//temporary[i - 1] = reference[i - 1];
		}
		for (j = mid; j < high; j++)
		{
			for(int m=0;m < this->nDims; m++) temporary[mid + high - j][m] = reference[j + 1][m];
			//temporary[mid + high - j] = reference[j + 1];
		}
		for (k = low; k <= high; k++)
		{
			pointTemp = (BKDTree::super_key_compare(temporary[i], temporary[j], _cur_dim) < 0) ? temporary[i++] : temporary[j--];
			for(int m = 0; m < this->nDims; m++) reference[k][m] = pointTemp[m];
			//reference[k] = (BKDTree::super_key_compare(temporary[i], temporary[j], _cur_dim) < 0) ? temporary[i++] : temporary[j--];
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
		else
		{
			std::cout<< "Duplicated !!"<<std::endl;
		}
	}
	return end;
}


void BKDTree::insert_leaf_data(Node* node, double* _point, const int _depth, const int _index) {// in case of reaching to a leaf, make the leaf node as a leaf node containing several points.
	
	if (node->isLeaf == true) // In case of leaf node, push_back the new point into the node.
	{ 
		PointNode* point_node = NULL ;//= (PointNode*) malloc(sizeof(PointNode));
		point_node = BKDTree::new_point_node(_point, _index);
		node->pointNodes.push_back(point_node);
		node->numOfPoints++;

		//std::cout << "reach the leaf, " << node << " , " << node->numOfPoints << ", size : " << node->pointNodes.size() << std::endl;
		// memory problem? Linux shows no problems.
	}
	else // traveling nodes until reaching the leaf node.
	{
		if (_point[_depth] <= node->refPoint[_depth]) // Go to left child
		{
			BKDTree::insert_leaf_data(node->left,  _point, (_depth + 1) % this->nDims, _index);
		} else {// Go to right child
			BKDTree::insert_leaf_data(node->right, _point, (_depth + 1) % this->nDims, _index);
		}
	}
}


Node* BKDTree::build_tree_recursively(double*** references, double** temp, const int start, const int end,  const int _depth) {
	// references : already sorted ! 
	// initial values : start ( 0 ), end ( npoints - 1 )
	Node* node = NULL; // root node
	//int axis = (_depth < this->nDims) ? _depth : _depth - this->nDims; 
	int axis = _depth % this->nDims;

	//std::cout<<"axis: "<<axis<<std::endl;
	//printf("current depth : %d, max depth : %d, axis : %d\n", _depth, this->maxDepth,axis);

	if(_depth < this->maxDepth) // If the depth does not reach to max_depth.
	{
		if (end == start) // 1�??�았????		
        {
			std::cout << "end == start" << std::endl;
			node        = BKDTree::new_node_leaf(-1);
			this->nNodes += 1;
		}
		else if (end == start + 1) // 2�??�았????		
        {
			//std::cout<<"end == start+1"<<std::endl;
			//for(int i =0; i < this->nDims; i++) std::cout<<references[0][start][i]<<", ";
			//std::cout<<std::endl;
			
			node        = BKDTree::new_node(references[0][start]);
			node->right = BKDTree::new_node_leaf(-2);
			node->left  = BKDTree::new_node_leaf(-2);
			
			//std::cout<<node<<", "<<node->right<<","<<node->left<<std::endl;

			this->nNodes += 3;
		}
		else if (end == start + 2) // 3�??�았????		
        {
			//std::cout<<"end == start+2"<<std::endl;
			//for(int i = 0; i < this->nDims; i++) std::cout << node->right->refPoint[i] << ", ";
			//std::cout<<std::endl;

			node        = BKDTree::new_node(references[0][start + 1]);
			node->right = BKDTree::new_node_leaf(-3.0);
			node->left  = BKDTree::new_node_leaf(-3.0);

			this->nNodes += 3;
		}
		else if (end >  start + 2) // 4�??�상 ?�았????		
        {
			//std::cout<<"end > start+2"<<std::endl;
			const int median = start + ( (end - start) / 2); // ?�림?�다. ex) (0+5)/2 = 2.5 (x) / 2 (o) , (0+4)/2 = 2
			node = BKDTree::new_node(references[0][median]);

			this->nNodes += 1;

			//std::cout<<"reference median: "<<std::endl;
			//for(int k = 0; k < this->nDims; k++) std::cout<<references[0][median][k]<<std::endl;

			for (int i = start; i <= end; i++)
			{
				//for(int m = 0; m < this->nDims; m++)
				//{
				//	temp[i][m] = references[0][i][m];
				//}
				temp[i] = references[0][i]; //  ?�기�?문제 ?�어 보임.
				//std::cout<<"temp add.: "<<temp[i]<<", ref. add. : "<<references[0][i]<<std::endl;
			}

			int lower, upper, lowerSave, upperSave;

			for(int i = 1; i < this->nDims; i++)
			{
				lower = start - 1;
				upper = median;
				
				for (int j = start; j <= end; j++)
				{
					double compare = BKDTree::super_key_compare(references[i][j], references[0][median], axis);
					//std::cout<<"cur dim: "<<i <<" ref : "<<references[i][j][axis]<<", current : "<<pointTemp[axis] <<" COMPARE : "<<compare<<" === ";
					if (compare < 0) // 만약, references[i][j]�??�재 중심?�인 node->point 보다 모든 차원?�서 ?�다�? (??차원?�서�??�아??
					{
						//std::cout<<"   down  "<<std::endl;
						references[i - 1][++lower] = references[i][j]; // 
					}
					else if (compare > 0) // 만약, references[i][j]�??�재 중심?�인 node->point 보다 모든 차원?�서 ?�다�? (??차원?�서�?커도)
					{
						//std::cout<<"   up    "<<std::endl;
						references[i - 1][++upper] = references[i][j];
					}
					else
					{	
						//std::cout<<"              SAME  "<<std::endl;
						//std::cout<<references[i][j][axis]<<", "<<pointTemp[axis]<<std::endl;
					}
				}

				// check the new indices for the reference array.
				if(lower < start || lower >= median)
				{
					std::cout<<"1. incorrect range for lower at depth = "<<_depth<<": start = "<<start<<"  lower = "<<lower<<"  median = "<<median<<std::endl;
					exit(1);
				}
				if( upper <= median || upper > end)
				{
					std::cout<<"2. incorrect range for upper at depth = "<<_depth<<": median = "<<median<<"  upper = "<<upper<<"  end = "<<end<<std::endl;
					exit(1);
				}
				if( i > 1 && lower != lowerSave)
				{
					std::cout<<"3. lower = "<<lower<<" != lowerSave = " <<lowerSave<<std::endl;
				}
				if(i > 1 && upper != upperSave)
				{
					std::cout<<"4. upper = "<<upper<<" != upperSave = "<<upperSave<<std::endl;
				}
				lowerSave = lower;
				upperSave = upper;
			}

			for (int i = start; i <= end; i++)
			{
				//for(int m = 0; m < this->nDims; m++)
				//{
				//	references[this->nDims - 1][i][m] = temp[i][m];
				//}
				references[this->nDims - 1][i] = temp[i];
			}
			//std::cout<<" LEFT start: "<<start<<", lower: "<<lower<<", depth: "<<_depth<<std::endl;
			node->left  = BKDTree::build_tree_recursively(references, temp, start, lower, _depth + 1);
			//std::cout<<" RIGHT median+1: "<<median+1<<", upper: "<<upper<<", depth: "<<_depth<<std::endl;
			node->right = BKDTree::build_tree_recursively(references, temp, median + 1, upper, _depth + 1);
		}
	}
	else // In case of the leaf node.
	{
		//printf(" [LEAF NODE], depth: %d\n",_depth);
		//printf("Leaf exits\n");

		node = BKDTree::new_node_leaf(-9);
		
		this->nNodes += 1;
	}

	return node; // tree root �� �� ���̴�.
}

Node* BKDTree::create_tree(double** _points) {
	double*** references = new double**[this->nDims];
	double*** referencesFixedOrder = new double**[this->nDims];
	double**  temp       = new double*[this->nPoints]; 
	double**  tempFixedOrder = new double*[this->nPoints];
	//double*** references = (double***)malloc(this->nDims*sizeof(double**));
	//double**  temp       = (double**)malloc(this->nPoints*sizeof(double*));

	for (int i = 0; i < this->nPoints; i++)
	{
		temp[i] = new double[this->nDims];
		tempFixedOrder[i] = temp[i];
	}

	for (int i = 0; i < this->nDims; i++) // dimension
	{   
		references[i] = new double*[this->nPoints];
		referencesFixedOrder[i] = new double*[this->nPoints];
		for (int j = 0; j < this->nPoints; j++)
		{
			references[i][j] = new double[this->nDims];
			referencesFixedOrder[i][j] = references[i][j];
		}

		BKDTree::initialize_reference(_points, references[i]);
		BKDTree::merge_sort(references[i], temp, (long)0, (long)this->nPoints - 1, i);
	}
	
	// display sorted array.
	for(int i = 0; i < this->nDims; i++)  
	{
		//std::cout<<"  current dim : "<<i<<std::endl;
		//for(int j = 0; j < this->nPoints; j++) printf("%lf\n", references[i][j][i]);
		//printf("\n");
	}
	

	int* end = new int[this->nDims];
	for (int i = 0; i < this->nDims; i++) // ?�리 구성 �? ?�전 ?�일??????��.
	{ 
		end[i] = BKDTree::remove_duplicates(references[i], i);
	}

	for (int i = 0; i < this->nDims - 1; i++) // 만약, 병합?�렬 ??모든 차원?�서??end�?같�? ?�으�? ?��?�?removal???�루?��?�??��? 것이�?�� ?�러발생.
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

	Node* root = NULL;
	root = BKDTree::build_tree_recursively(references, temp, 0, end[0], 0); //

	// insert point data into leaf nodes.
	for (int i = 0; i < this->nPoints; i++) 
	{
		/*std::cout<<"point : ";
		for (int j = 0; j < this->nDims; j++)
		{
			std::cout<<_points[i][j]<<", ";
		}
		std::cout<<std::endl;*/
		
		BKDTree::insert_leaf_data(root, _points[i], 0, i);
	}
	
	// delete memories.
	double_array_delete_3d(referencesFixedOrder, this->nDims, this->nPoints);
	for(int i = 0; i < this->nDims; i++) delete[] references[i];
	delete[] references;

	printf("[INFO] DELETION - this->references : deleted.\n");
	//double_array_delete_2d(temp, this->nPoints);
	for(int i = 0; i < this->nPoints; i++) delete[] tempFixedOrder[i];
	delete[] tempFixedOrder;
	delete[] temp;

	delete[] end;
	printf("[INFO] in create_tree : memories are deleted.\n");
	
	return root;
}


void BKDTree::delete_node_recursively(Node* node) {
	// delete all malloc s in the nodes. Following, delete nodes.
	//std::vector<double> nullTemp1, nullTemp2, nullTemp3;

	if(node->isLeaf != true) {
		delete_node_recursively(node->left);
		delete_node_recursively(node->right);
		//node->left->refPoint.swap(nullTemp1);
		//node->right->refPoint.swap(nullTemp2);
		//node->refPoint.swap(nullTemp3);
		
		//delete node->left;
		//delete node->right;
		//delete node;
	} else { // in the leaf node,
		for(int i = 0; i < node->numOfPoints; i++) delete node->pointNodes[i]; // delete all 'pointNodes' in the node. 
	}
	delete node;
	// printf("All memories of nodes and pointNodes are returned !\n");
}


void BKDTree::search_NN_depth_first_search(Node* node, const std::vector<double>& _point_q, const int _depth, int* _index) {
	if (node->isLeaf == true) // In case of leaf node, push_back the new point into the node.
	{ 
		double minDist = 1e30;
		*_index = -1;

		for (int k = 0; k<node->numOfPoints; k++)
		{
			//std::cout<<node->pointNodes[k]->index<<std::endl;
			double currDist = distance_euclidean(_point_q, node->pointNodes[k]->point, this->nDims, false);
			//double currDist = distance_manhattan(_point_q, node->pointNodes[k]->point, this->nDims);
			if (currDist < minDist)
			{
				minDist = currDist;
				*_index = node->pointNodes[k]->index;

				//if (minDist <= this->distThres) return;
			}
			//std::cout<<*_index<<std::endl;
		}
		if(minDist >= this->distThres) *_index = -1;
	}
	else // traveling nodes until reaching the leaf node.
	{ 
		if (_point_q[_depth] <= node->refPoint[_depth]) // Go to left child
		{ 
			BKDTree::search_NN_depth_first_search(node->left,  _point_q, (_depth + 1) % this->nDims,_index);
		}
		else // Go to right child
		{	
			BKDTree::search_NN_depth_first_search(node->right, _point_q, (_depth + 1) % this->nDims,_index);
		}
	}
	return ;
}



void BKDTree::search_NN_depth_first_search_multiple(Node* _node, const std::vector<std::vector<double>>& _points_q_vec, std::vector<int>& _index_vec){

	int* indexTemp = new int();
	int nPointsQuery = _points_q_vec.size();
	int nDimsQuery    = _points_q_vec[0].size();
	
	//std::vector<int> index_vec;
	//index_vec.reserve(nPointsQuery);

	if(nDimsQuery != this->nDims)
	{
		printf("Error - query point dimension does not match reference points ! \n");
		exit(1);
	}

	for(int i = 0; i < nPointsQuery; i++)
	{
		BKDTree::search_NN_depth_first_search(_node, _points_q_vec[i], 0, indexTemp);
		//index_vec.push_back(*indexTemp);
		_index_vec.push_back(*indexTemp);
	}

	delete indexTemp;
}


void BKDTree::kdtree_nearest_neighbor_approximate(const std::vector<std::vector<double>>& _points_q_vec,  std::vector<int>& _index_vec){
	this->search_NN_depth_first_search_multiple(this->treeRootNode, _points_q_vec, _index_vec);
}


void BKDTree::print_points(const std::vector<double>& _point, const bool _isLeaf) {

	printf("[ %0.0lf,", _point[0]);
	for (int i = 1; i<this->nDims - 1; i++) printf("%0.0lf,", _point[i]);
	printf("%0.0lf ]", _point[this->nDims - 1]);

	if (_isLeaf == true) {
		printf(" < LF >");
	}
	else {
		printf(" < NE >");
	}
}

void BKDTree::print_tree(const Node* node, const int _depth) {
	if (node != NULL) {
		if(node->right != NULL) BKDTree::print_tree(node->right, _depth + 1);

		for (int i = 0; i < _depth; i++) {
			printf("                       ");
		}
		//BKDTree::print_points(aa, node->isLeaf);
		BKDTree::print_points(node->refPoint, node->isLeaf);
		printf("\n");
		if(node->left != NULL)		BKDTree::print_tree(node->left, _depth + 1);
	}
}

void BKDTree::print_leaf(const Node* node, const int _depth) {
	if (node != NULL) {
		BKDTree::print_leaf(node->right, _depth + 1);
		if (node->isLeaf) {
			for (int i = 0; i < node->numOfPoints; i++) {
				std::cout << "coordinate : [";
				for (int j = 0; j < this->nDims; j++) {
					std::cout << node->pointNodes[i]->point[j] << ",";
				}
				std::cout << " ], index : " << node->pointNodes[i]->index << std::endl;
			}
		}
		BKDTree::print_leaf(node->left, _depth + 1);
	}
}
