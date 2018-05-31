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
	std::uniform_real_distribution<> distribution(10, 515);
	auto generator = std::bind(distribution, engine);


	int numOfPoints = 5000;
	int ndim = 4;

	std::vector<std::vector<double>> points_vec;
	for (int i = 0; i<numOfPoints; i++) {
		std::vector<double> temp;
		for (int j = 0; j<ndim; j++) {
			temp.push_back(generator());
		}
		points_vec.push_back(temp);
	}

  std::vector<std::vector<double>> points_q_vec;
  for (int i = 0; i<numOfPoints/10; i++) {
    std::vector<double> temp;
    for (int j = 0; j<ndim; j++) {
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
	int binSize = 50; // bin size �� �� 10~50 ������ �����ѵ�.
	int max_depth = (int)round(log2(numOfPoints / binSize));
	std::cout << "max depth : " << max_depth << std::endl;

  BKDTree* kdtree_handler = new BKDTree(points_vec, binSize);
  Node* tree_root;
  tree_root = kdtree_handler->create_tree(points);


  start = clock();

  std::vector<int> indexVec;
  indexVec.reserve(numOfPoints);
  for(int j = 0; j < 30 ; j++){
    kdtree_handler->search_NN_dfs_index(tree_root, points_q_vec,indexVec);
  } 

	/*Node* tree_root;

	tree_root = create_tree(points, numOfPoints, ndim, max_depth); // 5000 points : 6 ms

for(int j = 0; j<10; j++){
for (int i = 0; i < numOfPoints/15 ; i++) { // 5000 queries : 0.1 ms . Therefore, 30 iterations corresponds to 3 ms
		int* a = new int;
		search_NN_dfs(tree_root, points[i], 0, ndim,a);
		//if(i != *a) std::cout<<" real index : "<<i<<",  estimated index : "<<*a<<std::endl;
		for (int j = 0; j< ndim; j++) {
			//std::cout<<"min dist point : "<<point_node->point[j]<<",";
		}
		for (int j = 0; j< ndim; j++) {
			//std::cout<<points[i][j]<<",";
		}
		//std::cout<<std::endl;

	}

}*/

	finish = clock();

	std::cout << "Elapsed time : " << (finish - start) / 1000000.0 << std::endl;
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
