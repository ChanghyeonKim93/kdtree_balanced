#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "kdtree_kch.h"

#include <random>
#include <ctime>
#include <functional>

typedef double* Point;
typedef std::vector<Point> Points;

int main(){
    clock_t start, finish;
    std::mt19937 engine((unsigned int)time(NULL));
    std::uniform_real_distribution<> distribution(10,515);
    auto generator = std::bind(distribution,engine);


    int numOfPoints = 5000;
    int ndim = 3;

    std::vector<std::vector<double>> points_vec;
    for(int i=0; i<numOfPoints; i++){
        std::vector<double> temp;
        for(int j= 0; j<ndim; j++){
            temp.push_back(generator());
	    //temp.push_back(j);
        }
        points_vec.push_back(temp);
    }

    double** points = (double**)malloc(numOfPoints*sizeof(double*));

    for(int i = 0; i < numOfPoints; i++){
        points[i] = (double*)malloc(ndim*sizeof(double));
        for(int j=0; j< ndim; j++){
            points[i][j] = points_vec[i][j];
        }
    }

    // tree generation
    int binSize = 200;
    int max_depth = round(log2(numOfPoints/binSize));
    std::cout<<"max depth : "<<max_depth<<std::endl;
    Node* tree_root;
    //int max_depth = 4;


    tree_root=create_tree(points,numOfPoints,ndim, max_depth); // 5000 points : 6 ms

    start  = clock();
    for(int i = 0; i < numOfPoints; i++){ // 5000 queries : 0.1 ms . Therefore, 30 iterations corresponds to 3 ms
	int a= search_NN_dfs(tree_root, points[i], 0, ndim);
	for(int j = 0; j< ndim; j++){
	    //std::cout<<"min dist point : "<<point_node->point[j]<<",";
	}
	for(int j = 0; j< ndim; j++){
	    //std::cout<<points[i][j]<<",";
	}
	//std::cout<<std::endl;
	
    }
    finish = clock();

    std::cout<<"Elapsed time : "<< (finish - start)/1000000.0 <<std::endl;
    std::cout<<"The number of nodes : "	<<verify_tree(tree_root, ndim, 0)<<std::endl;
    //
    //print_tree(tree_root,ndim,0);
    //print_leaf(tree_root,ndim,0);

    for(int i =0; i < numOfPoints; i++){
        free(points[i]);
    }
    free(points);
    std::cout<<"Program is running !"<<std::endl;
    return 0;
}