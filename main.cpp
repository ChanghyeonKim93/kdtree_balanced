#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "BKDTree.h"
//#include "kdtree_kch.h"
//#include "kdtree_kch_class.h"

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
	std::uniform_real_distribution<> distribution(1.0, 1000.0);
	auto generator = std::bind(distribution, engine);


	int numOfPoints  = 200;
  int numOfQueries = 1;
	int ndim         = 4;

	int binSize      = 1; // bin size �� �� 10~50 ������ �����ѵ�.
	int maxDepth     = (int)ceil( log2(  ceil( (double)numOfPoints / (double)binSize ) ) );
  double distThres = 20;

	std::vector<std::vector<double>> points_vec;
	for (int i = 0; i < numOfPoints; i++)
  {
		std::vector<double> temp;
		for (int j = 0; j < ndim; j++)
    {
      temp.push_back(generator());
      std::cout<<temp[j]<<", ";
    }
    std::cout<<std::endl;
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
  std::cout<< "max depth : " << maxDepth << std::endl;
  std::cout<< "Required # of bins: "<<ceil( (double)numOfPoints / (double)binSize )<<", Total capacity: "<<pow(2.0,maxDepth)*binSize<<std::endl;

  start  = clock();
  BKDTree* tree = new BKDTree(points_vec, binSize, distThres);
  finish = clock();

  std::cout << "Tree gen elapsed time : " << (finish - start) / 1000.0 << "[ms]"<< std::endl;


  start  = clock();
  std::vector<int> indexVec;
  indexVec.reserve(numOfPoints);
  //tree->kdtree_nearest_neighbor_approximate(points_q_vec, indexVec);
  //std::cout<<"index size:"<<indexVec.size()<<std::endl;
  //for(int i = 0; i < indexVec.size(); i++) std::cout<< indexVec[i] <<std::endl;
  finish = clock();

	std::cout << "Search elapsed time : " << (finish - start) / 1000.0 << "[ms]"<< std::endl;
	//
  tree->print_tree(tree->treeRootNode,0);
	//tree->print_leaf(tree->treeRootNode,0);


	for (int i = 0; i < numOfPoints; i++) {
		free(points[i]);
	}
	free(points);

  delete tree;
	std::cout << "Program is running !" << std::endl;
	return 0;
}
