#include "BKDTree.h"
#include "mex.h"
#include "matrix.h" //isNaN/isinf
#include <vector>
void retrieve_tree( const mxArray* matptr, BKDTree* & tree){
    // retrieve pointer from the MX form
    double* pointer0 = mxGetPr(matptr);
    // check that I actually received something
    if( pointer0 == NULL )
        mexErrMsgTxt("vararg{1} must be a valid k-D tree pointer\n");
    // convert it to "long" datatype (good for addresses)
    intptr_t pointer1 = (intptr_t) pointer0[0];
    // convert it to "KDTree"
    tree = (BKDTree*) pointer1;
    // check that I actually received something
    if( tree == NULL )
        mexErrMsgTxt("vararg{1} must be a valid k-D tree pointer\n");
    if( tree -> nDims <= 0 )
        mexErrMsgTxt("the k-D tree must have k>0"); 
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){   
	// check input
	if( nrhs != 1 )
		mexErrMsgTxt("A address to tree should be passed.\n");
    
    // retrieve the tree pointer
    BKDTree* tree;

    retrieve_tree( prhs[0], tree ); 
    delete tree;
    printf("The tree is deleted.\n");
}