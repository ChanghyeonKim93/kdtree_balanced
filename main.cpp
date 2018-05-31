#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
//#include "kdtree_kch.h"
#include "kdtree_kch_class.h"

#include <random>
#include <ctime>
#include <functional>

typedef double* Point;
typedef std::vector<Point> Points;

int main() {

  clock_t start, finish;
	std::mt19937 engine((unsigned int)time(NULL));
	std::uniform_real_distribution<> distribution(1, 1000);
	auto generator = std::bind(distribution, engine);


	int numOfPoints  = 5000;
  int numOfQueries = 1000;
	int ndim =4;

	std::vector<std::vector<double>> points_vec;
	for (int i = 0; i < numOfPoints; i++) {
		std::vector<double> temp;
		for (int j = 0; j < ndim; j++) {
			temp.push_back(generator());
		}
		points_vec.push_back(temp);
	}

  std::vector<std::vector<double>> points_q_vec;
  for (int i = 0; i < numOfQueries; i++) {
    std::vector<double> temp;
    for (int j = 0; j < ndim; j++) {
      temp.push_back(generator());
      //temp.push_back(j);
    }
    points_q_vec.push_back(temp);
  }

	double** points = (double**)malloc(numOfPoints*sizeof(double*));

	for (int i = 0; i < numOfPoints; i++) {
		points[i] = (double*)malloc(ndim*sizeof(double));
		for (int j = 0; j< ndim; j++) {
			points[i][j] = points_vec[i][j];
		}
	}

	// tree generation
	int binSize = 20; // bin size �� �� 10~50 ������ �����ѵ�.
	int max_depth = (int)round(log2(numOfPoints / binSize));
  double dist_thres = 5;
  std::cout << "max depth : " << max_depth << std::endl;

  start = clock();
  BKDTree* kdtree_handler = new BKDTree(points_vec, binSize,dist_thres);
  Node* tree_root;
  tree_root = kdtree_handler->create_tree(points);
  finish = clock();
  std::cout << "Tree gen elapsed time : " << (finish - start) / 1000.0 << "[ms]"<< std::endl;


  start = clock();
  std::vector<int> indexVec;
  indexVec.reserve(numOfPoints);
  for(int j = 0; j < 50 ; j++){
    int* index = new int;
    kdtree_handler->search_NN_dfs_index(tree_root, points_q_vec,indexVec);
    /*for(int i = 0; i<numOfQueries;i++){
      kdtree_handler->search_NN_dfs(tree_root, points[j], 0, index);
    }*/
  }
  finish = clock();

	std::cout << "Search elapsed time : " << (finish - start) / 1000.0 << "[ms]"<< std::endl;
	//
	//print_tree(tree_root,ndim,0);
	//print_leaf(tree_root,ndim,0);


	for (int i = 0; i < numOfPoints; i++) {
		free(points[i]);
	}
	free(points);
  delete kdtree_handler;
	std::cout << "Program is running !" << std::endl;
	return 0;
}
