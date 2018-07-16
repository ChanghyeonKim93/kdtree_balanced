#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <random>
#include <ctime>
#include <functional>
#include "BKDTree.h"

void print_start();

int main() {
  print_start();

  clock_t start, finish;
	std::mt19937 engine((unsigned int)time(NULL));
	std::uniform_real_distribution<> distribution(1.0, 1000.0);
	auto generator = std::bind(distribution, engine);

	int numOfPoints  = 5000;
  int numOfQueries = 1;
	int ndim         = 4;

	int binSize      = 10; // bin size �� �� 10~50 ������ �����ѵ�.
  double distThres = 20;
	int maxDepth     = (int)ceil( log2(  ceil( (double)numOfPoints / (double)binSize ) ) );

	std::vector<std::vector<double>> points_vec;
	for (int i = 0; i < numOfPoints; i++)
  {
		std::vector<double> temp;
		for (int j = 0; j < ndim; j++)
    {
      temp.push_back(generator());
    }
		points_vec.push_back(temp);
	}

  std::vector<std::vector<double>> points_q_vec;
  for (int i = 0; i < numOfQueries; i++) 
  {
    std::vector<double> temp;
    for (int j = 0; j < ndim; j++) temp.push_back(points_vec[i][j]+0.05*(generator()-500));
    points_q_vec.push_back(temp);
  }


	// tree generation
  std::cout << "max depth : " << maxDepth << std::endl;
  std::cout << "Required # of bins: " << ceil( (double)numOfPoints / (double)binSize ) <<", Total capacity: " << pow(2.0,maxDepth)*binSize << std::endl;

  start  = clock();
  BKDTree* tree = NULL;
  tree = new BKDTree(points_vec, binSize, distThres);
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

  // tree->print_tree(tree->treeRootNode,0);
	//tree->print_leaf(tree->treeRootNode,0);

  delete tree;
	std::cout << "Program is running !" << std::endl;
	return 0;
}





/** 
 * 
 * functions
 * 
 */

void print_start()
{
  printf("\n----------------------------\n");
  printf("-                          -\n");
  printf("-      program starts.     -\n");
  printf("-                          -\n");
  printf("----------------------------\n\n\n");
}
