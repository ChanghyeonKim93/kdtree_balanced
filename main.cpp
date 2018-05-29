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
    std::mt19937 engine((unsigned int)time(NULL));
    std::uniform_real_distribution<> distribution(10,515);
    auto generator = std::bind(distribution,engine);


    int numOfPoints = 153;
    int ndim = 3;

    std::vector<std::vector<double>> points_vec;
    for(int i=0; i<numOfPoints; i++){
        std::vector<double> temp;
        for(int j= 0; j<ndim; j++){
            temp.push_back(generator());
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
    int binSize = 5;
    int max_depth = round(log2(numOfPoints/binSize));
    std::cout<<"max depth : "<<max_depth<<std::endl;
    Node* tree_root;
    //int max_depth = 4;


    tree_root=create_tree(points,numOfPoints,ndim, max_depth);
    std::cout<<"The number of nodes : "<<verify_tree(tree_root, ndim, 0)<<std::endl;
    //
    print_tree(tree_root,ndim,0);
    print_leaf(tree_root,ndim,0);

    for(int i =0; i < numOfPoints; i++){
        free(points[i]);
    }
    free(points);
    std::cout<<"Program is running !"<<std::endl;
    return 0;
}
