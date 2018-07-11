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
  printf("\n----------------------------\n");
  printf("-                          -\n");
  printf("-      program starts.     -\n");
  printf("-                          -\n");
  printf("----------------------------\n\n\n");

  clock_t start, finish;
	std::mt19937 engine((unsigned int)time(NULL));
	std::uniform_real_distribution<> distribution(1, 1000);
	auto generator = std::bind(distribution, engine);


	int numOfPoints  = 6000;
  int numOfQueries = 1000;
	int ndim = 4;

	std::vector<std::vector<double>> points_vec;
	for (int i = 0; i < numOfPoints; i++) 
  {
		std::vector<double> temp;
		for (int j = 0; j < ndim; j++) temp.push_back(generator());
		points_vec.push_back(temp);
	}

  std::vector<std::vector<double>> points_q_vec;
  for (int i = 0; i < numOfQueries; i++) 
  {
    std::vector<double> temp;
    for (int j = 0; j < ndim; j++) temp.push_back(points_vec[i][j]+0.05*(generator()-500));
    points_q_vec.push_back(temp);
  }

	double** points = (double**)malloc(numOfPoints*sizeof(double*));
	for (int i = 0; i < numOfPoints; i++) 
  {
		points[i] = (double*)malloc(ndim*sizeof(double));
		for (int j = 0; j< ndim; j++) points[i][j] = points_vec[i][j];
	}

	// tree generation
	int binSize = 20; // bin size �� �� 10~50 ������ �����ѵ�.
	int maxDepth = (int)ceil( log2(  ceil( (double)numOfPoints / (double)binSize ) ) );
  double distThres = 5;
  std::cout<< "max depth : " << maxDepth << std::endl;
  std::cout<< "Required # of bins: "<<ceil( (double)numOfPoints / (double)binSize )<<", Total capacity: "<<pow(2.0,maxDepth)*binSize<<std::endl;

  start = clock();
  BKDTree* bkdTree2 = new BKDTree(points_vec, binSize, distThres);
  //Node* tree_root;
  //tree_root = kdtree_handler->create_tree(points);
  finish = clock();
  std::cout << "Tree gen elapsed time : " << (finish - start) / 1000.0 << "[ms]"<< std::endl;


  start = clock();
  std::vector<int> indexVec;
  indexVec.reserve(numOfPoints);
  bkdTree2->kdtree_nearest_neighbor_approximate(points_q_vec, indexVec);
  finish = clock();


  /*for(int j = 0; j < 500 ; j++){
    int* index = new int;
    kdtree_handler->search_NN_dfs_index(tree_root, points_q_vec,indexVec);
    //for(int i = 0; i<numOfQueries;i++){
    //  kdtree_handler->search_NN_dfs(tree_root, points[j], 0, index);
    //}
  }*/

  std::cout<<"index size:"<<indexVec.size()<<std::endl;
  for(int i = 0; i < indexVec.size(); i++) std::cout<<indexVec[i]<<std::endl;

	std::cout << "Search elapsed time : " << (finish - start) / 1000.0 << "[ms]"<< std::endl;
	//
	//bkdTree2->print_tree(bkdTree2->treeRootNode,0);
	//bkdTree2->print_leaf(bkdTree2->treeRootNode,0);


	for (int i = 0; i < numOfPoints; i++) {
		free(points[i]);
	}
	free(points);
  delete bkdTree2;
	std::cout << "Program is running !" << std::endl;
	return 0;
}
