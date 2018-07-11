% BKDTREE_BUILD construct a bkd-tree from a point cloud
%
% SYNTAX
% tree = bkdtree_build(p)
%
% INPUT PARAMETERS
%   P: a set of N k-dimensional points stored in a 
%      NxK matrix. (i.e. each row is a point)
%
% OUTPUT PARAMETERS
%   tree: a pointer to the created data structure
%
% DESCRIPTION
% Given a point set p, builds a k-d tree as specified in [1] 
% with a preprocessing time of O(d N logN), N number of points, 
% d the dimensionality of a point
% 
%
%
% Copyright (c) 2018 Changhyeon Kim
% All Rights Reserved
% email: rlackd93@snu.ac.kr
% $Revision: 1.0$  Created on: 2018/07/12