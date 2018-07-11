% BKDTREE_NEAREST_NEIGHBOR_MEX query a bkd-tree for nearest neighbor
%
% SYNTAX
% [idxs,dst] = bkdtree_nearest_neighbor_mex( tree, P )
%
% INPUT PARAMETERS
%   tree: a pointer to the previously constructed k-d tree
%   P: a set of N k-dimensional points stored in a 
%      NxK matrix (i.e. each row is a point). For each of these
%      points a kd-tree query is executed and the index
%      of the closest point is stored in the n-th position in the
%      output
%
% OUTPUT PARAMETERS
%   idxs: a column vector of scalars that index the point database.
%         In k-th position the index of the point in the database
%         closest to P(k,:) can be found.
%

% Copyright (c) 2018 Changhyeon Kim
% All Rights Reserved
% email: rlackd93@snu.ac.kr
% $Revision: 1.0$  Created on: 2018/07/12
