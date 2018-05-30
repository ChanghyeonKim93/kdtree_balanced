#include <iostream>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <stdlib.h>

typedef double* Point;
typedef std::vector<Point> Points;
typedef std::vector<int> Indices;

struct PointNode{
    double* point;
    int index;
};

struct Node{
    double* point;
    Node* left;
    Node* right;

    Points leafData;
    std::vector<PointNode*> pointNodes;
    int numOfPoints;
    bool isLeaf;
};

PointNode* new_point_node(double* _point, int _index){
    PointNode* point_node;
    if( (point_node = (PointNode*)malloc(sizeof(PointNode)))==NULL){
        printf("error allocating new PointNode ! \n");
        exit(1);
    }

    point_node->point = _point;
    point_node->index = _index;

    return point_node;
}

Node* new_node(double* _point){
    Node* node;

    if( ( node = (Node*)malloc(sizeof(Node))) == NULL){
        printf("error allocating new Node ! \n");
        exit(1);
    }
    // initialize the data w/o leafData
    node->point = _point;

    node->left = node->right = NULL;

    node->leafData.resize(0);
    node->pointNodes.resize(0);
    node->numOfPoints = 0;
    node->isLeaf = false;

    return node;
}

void initialize_reference(double** _points, double** reference, const int _npoints){
    for(int i=0; i<_npoints; i++){
        reference[i] = _points[i];
    }
}

double super_key_compare(const double* _point_a, const double* _point_b, const int _cur_dim, const int _ndim){
    double diff = 0;
    for(int i = 0; i< _ndim; i++){
        int r = i + _cur_dim;
        r = ( r < _ndim) ? r : r - _ndim;
        diff = _point_a[r] - _point_b[r];
        if( diff != 0 ) break; // diff가 0이 아니면 멈춤
    }
    return diff;
}

void merge_sort(double** reference, double** temporary, const int low, const int high, const int _cur_dim, const int _ndim){
    int i,j,k;

    if(high>low){
        const int mid = (high + low)/2;
        merge_sort(reference, temporary, low,   mid,  _cur_dim, _ndim);
        merge_sort(reference, temporary, mid+1, high, _cur_dim, _ndim);

        for(i = mid+1; i > low; i--){
            temporary[i-1] = reference[i-1];
        }

        for(j = mid; j < high; j++){
            temporary[mid+high-j] = reference[j+1];
        }
        for(k = low; k <= high; k++){
            reference[k] = (super_key_compare(temporary[i], temporary[j], _cur_dim, _ndim) < 0 ) ? temporary[i++] : temporary[j--];
        }
    }
}

int remove_duplicates(double** reference, const int _npoints, const int _cur_dim, const int _ndim){
    int end = 0;
    for(int j = 1; j < _npoints; j++){
        double compare = super_key_compare(reference[j], reference[j-1], _cur_dim, _ndim);
        if( compare < 0){
            printf("merge sort failure: super_key_compare(ref[%d], ref[%d]), current lvl=(%d), compare value = %d\n", j, j-1, _cur_dim, compare);
            exit(1);
        }
        else if(compare > 0) reference[++end] = reference[j];
    }

    return end;
}

Node* build_tree(double*** references, double** temp, const int start, const int end, const int _ndim, const int _depth, const int _max_depth){

    Node* node;
    int axis = _depth % _ndim;
    //printf("current depth : %d, max depth : %d\n", _depth, _max_depth);

    if(_depth < _max_depth){ // max_depth 까지만 node를 추가한다.
        if(end == start){
            node = new_node( references[0][end]);
        }
        else if(end == start + 1){
            node = new_node(references[0][start]);
            node->right = new_node(references[0][end]);
        }
        else if(end == start + 2){
            node = new_node(references[0][start+1]);
            node->left = new_node(references[0][start]);
            node->right = new_node(references[0][end]);
        }
        else if(end > start + 2){
            const int median = (start + end ) / 2;
            node = new_node(references[0][median]);

            for (int i = start; i <=end; i++){
                temp[i] = references[0][i];
            }

            int lower, upper, lowerSave, upperSave;

            for(int i = 1; i < _ndim; i++){
                lower = start - 1;
                upper = median;
                for(int j = start; j <=end; j++) {
                    double compare = super_key_compare(references[i][j], node->point, axis, _ndim);
                    if( compare < 0){
                        references[i-1][++lower] = references[i][j];
                    }
                    else if( compare > 0){
                        references[i-1][++upper] = references[i][j];
                    }
                }
                lowerSave = lower;
                upperSave = upper;
            }
            for(int i = start; i<= end; i++){
                references[_ndim-1][i] = temp[i];
            }
            node->left  = build_tree(references, temp, start,    lower, _ndim, _depth + 1, _max_depth);
            node->right = build_tree(references, temp, median+1, upper, _ndim, _depth + 1, _max_depth);
        }
    }
    else{
        double* temp = (double*)malloc(sizeof(double)*_ndim);
        for(int i = 0; i<_ndim; i++) temp[i]=0;
        node = new_node(temp);
        node->left = node->right = NULL;
        node->isLeaf = true;
        //printf("Leaf exit\n");
    }

    return node; // tree root 가 될 것이다.
}

void insert_leaf_data(Node* node, double* _point, const int _depth, const int _ndim, const int _index){// in case of reaching to a leaf, make the leaf node as a leaf node containing several points.

    if(node->isLeaf){ // In case of leaf node, push_back the new point into the node.
        node->leafData.push_back(_point);

        PointNode* point_node;
        point_node = (PointNode*)malloc(sizeof(PointNode));
        point_node->point = _point;
        point_node->index = _index;
        node->pointNodes.push_back(point_node);

        node->numOfPoints++;
        //std::cout<<"reach the leaf, "<<node<<" , "<<node->numOfPoints<<"," <<node->pointNodes.size()<<" : ";
    }

    else{ // traveling nodes until reaching the leaf node.
        if( _point[_depth] <= node->point[_depth] ){ // Go to left child
            insert_leaf_data(node->left,  _point, ( _depth + 1 ) % _ndim, _ndim, _index);
        }
        else{ // Go to right child
            insert_leaf_data(node->right, _point, ( _depth + 1 ) % _ndim, _ndim, _index);
        }
    }
}

Node* create_tree(double** _points, const int _npoints, const int _ndim, const int _max_depth){

    double*** references = (double***)malloc(_ndim*sizeof(double**));
    double** temp = (double**) malloc(_npoints*sizeof(double));
    for(int i = 0; i < _ndim; i++){   // dimension 만큼 merge 해서 references에 저장한다.
        references[i] = (double**)malloc(_npoints*sizeof(double*));
        initialize_reference(_points, references[i], _npoints);
        merge_sort(references[i], temp, 0, _npoints-1, i, _ndim);
    }

    int *end = (int*)malloc(_ndim*sizeof(int));
    for(int i = 0; i< _ndim; i++){ // 중복 된 것 을 삭제해준다.
        end[i] = remove_duplicates(references[i], _npoints, i, _ndim);
    }

    for(int i = 0; i < _ndim - 1; i++){
        for(int j = i + 1; j < _ndim; j++){
            if(end[i] != end[j]){
                printf("reference removal error in create_tree\n");
                exit(1);
            }
        }
    } // 여기까지는 문제가 없다.

    Node* root = build_tree(references, temp, 0, end[0], _ndim, 0, _max_depth); // 여기가 문제인가 ?

    for(int j = 0; j < _npoints; j++){
        insert_leaf_data(root, _points[j], 0, _ndim, j);
        //std::cout<<j<<std::endl;
    }

    for(int i = 0; i < _ndim; i++){
        free(references[i]);
    }
    free(references);
    free(temp);

    return root;
}

double distance_euclidean(const double* _point_q, const double* _point_r, const int _ndim){
    double temp = 0 ;
    for(int i = 0; i < _ndim; i++){
	temp += (_point_q[i]-_point_r[i])*(_point_q[i]-_point_r[i]);
    }
    return sqrt(temp);
}

int search_NN_dfs(Node* node, const double* _point_q, const int _depth, const int _ndim){
    PointNode* curr_point_node;
    int min_index=0;
    if(node->isLeaf){ // In case of leaf node, push_back the new point into the node.
	double min_dist = 100000000000;
	int min_index = -1;

	for(int k=0;k<node->numOfPoints; k++){
	    //std::cout<<node->pointNodes[k]->index<<std::endl;
	    double curr_dist = distance_euclidean(_point_q, node->pointNodes[k]->point, _ndim);
	    if(curr_dist < min_dist){
		min_dist  = curr_dist;
		min_index = k;
		curr_point_node = node->pointNodes[k];
	    }
	}
       	//std::cout<<"Searching is reaching the leaf. : "<<node<<", min dist idx : "<<min_index<<std::endl;
    }

    else{ // traveling nodes until reaching the leaf node.
        if( _point_q[_depth] <= node->point[_depth] ){ // Go to left child
            search_NN_dfs(node->left,  _point_q, ( _depth + 1 ) % _ndim, _ndim);
        }
        else{ // Go to right child
            search_NN_dfs(node->right, _point_q, ( _depth + 1 ) % _ndim, _ndim);
        }
    }

    return min_index;
}

int verify_tree(const Node* node, const int _ndim, const int _depth)
{
    int count = 1 ;
    if (node->point == NULL) {
    	printf("point is null\n");
    	exit(1);
    }
    // The partition cycles as x, y, z, w...
    int axis = _depth % _ndim;

    if (node->left != NULL) {
    	if (node->left->point[axis] > node->point[axis]) {
            //printf("child is > node!\n");
            //exit(1);
        }
        if (super_key_compare(node->left->point, node->point, axis, _ndim) >= 0) {
            //printf("child is >= node!\n");
            //exit(1);
        }
    	count += verify_tree(node->left, _ndim, _depth + 1);
    }
    if (node->right != NULL) {
    	if (node->right->point[axis] < node->point[axis]) {
            //printf("child is < node!\n");
            //exit(1);
        }
        if (super_key_compare(node->right->point, node->point, axis, _ndim) <= 0) {
            //printf("child is <= node!\n");
            //exit(1);
        }
    	count += verify_tree(node->right, _ndim, _depth + 1);
    }

    return count; // the number of nodes.
}

void print_points(const double* _point, const int _ndim, const bool _isLeaf) {
    printf("[ %0.0f,", _point[0]);
    for (int i=1; i<_ndim-1; i++) printf("%0.0f,", _point[i]);
    printf("%0.0f ]", _point[_ndim-1]);
    if(_isLeaf == true){
        printf("< LEAF >");
    }
    else{
        printf("<nd>");
    }
}

void print_tree(const Node* node, const int _ndim, const int _depth) {
    if (node!=NULL) {
    	print_tree(node->right, _ndim, _depth+1);

    	for (int i=0; i<_depth; i++){
            printf("             ");
    	}
    	print_points(node->point, _ndim, node->isLeaf);
    	printf("\n");
    	print_tree(node->left, _ndim, _depth+1);
    }
}

void print_leaf(const Node* node, const int _ndim, const int _depth) {
    if(node!=NULL){
        print_leaf(node->right, _ndim,_depth+1);
        if(node->isLeaf){
            for(int i = 0; i < node->numOfPoints; i++){
                std::cout<<"coordinate : [";
                for(int j = 0; j< _ndim; j++){
                    std::cout<<node->pointNodes[i]->point[j]<<",";
                }
                std::cout<<" ], index : "<<node->pointNodes[i]->index<<std::endl;
            }
        }
        print_leaf(node->left, _ndim, _depth+1);
    }
}
