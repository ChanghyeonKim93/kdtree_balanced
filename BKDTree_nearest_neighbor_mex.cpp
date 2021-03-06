#include "BKDTree.h"
#include "mex.h"

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
void retrieve_data( const mxArray* matptr, double*& data, int& npoints, int& ndims){	
	// retrieve pointer from the MX form
    data = mxGetPr(matptr);
    // check that I actually received something
    if( data == NULL )
        mexErrMsgTxt("vararg{2} must be a [kxN] matrix of data\n");
    // retrieve amount of points
    npoints = mxGetM(matptr);
    ndims   = mxGetN(matptr);
}
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){
	// check number of arguments
	if( nrhs!=2 )
		mexErrMsgTxt("This function requires 2 arguments\n");
	if( !mxIsNumeric(prhs[0]) )
		mexErrMsgTxt("varargin{0} must be a valid kdtree pointer\n");
	if( !mxIsNumeric(prhs[1]) )
		mexErrMsgTxt("varargin{1} must be a query set of points\n");
		
    // retrieve the tree pointer
    BKDTree* tree;

    retrieve_tree( prhs[0], tree ); 
    // retrieve the query data
    double* query_data;
    int npoints, ndims;
    retrieve_data( prhs[1], query_data, npoints, ndims );
    // printf("query size: %dx%d\n", npoints, ndims);
    
    // check dimensions
    if( ndims != tree->nDims ) 
    	mexErrMsgTxt("vararg{1} must be a [Nxk] matrix of N points in k dimensions\n");   
    
    // npoints x 1 indexes in output
    plhs[0] = mxCreateDoubleMatrix(npoints, 1, mxREAL);
    double* indexes = mxGetPr(plhs[0]);
    //plhs[1] = mxCreateDoubleMatrix(npoints, 1, mxREAL);
    //double* dists   = mxGetPr(plhs[1]);

    // cout << "nindexes: " << mxGetM(plhs[0]) << "x" << mxGetN(plhs[0]) << endl;
    
    // execute the query FOR EVERY point storing the index
    std::vector<std::vector<double>> points_q_vec;
    for(int i =0; i<npoints; i++){
        std::vector< double > query;
        for(int j=0;j<ndims;j++) query.push_back(query_data[i+j*npoints]); 
        points_q_vec.push_back(query);
    }

    std::vector<int> indexVec;
    indexVec.reserve(0);
    tree->kdtree_nearest_neighbor_approximate(points_q_vec, indexVec);
    for(int i=0; i<npoints; i++)
    {
        //printf("%d\n",indexVec[i]);
        indexes[i] = indexVec[i]+1; //M-idx
    }
}