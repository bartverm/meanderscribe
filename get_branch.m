function branch_idx=get_branch(cnod,conn) 
% returns indices of current node and all its parents, incl root node
%
%   BRANCH_IDX = get_branch(CNOD, CONN) computes all nodes parent of CNOD
%    based on the tree given in CONN (from scale_space_tree).
%
%   See also: scale_space_tree

% The MIT License (MIT)
% 
% Copyright (c) 2016 Bart Vermeulen
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

check_conn(conn);                                                          % Check validity of CONN
assert(isscalar(cnod) && cnod > 0 && mod(cnod,1)==0 &&...       
    cnod<=numel(conn),'CNOD should be a scalar pointing to a node in CONN')% Check validity of CNOD


branch_idx=[];                                                             % Initialize vector with indices of all nodes in the branch
while (conn(cnod)~=0)                                                      % while current node is not the root node
    branch_idx=[branch_idx; cnod];                             %#ok<AGROW>   store it (Vector growth is not really significant here)
    cnod=conn(cnod);                                                       % continue with parent node
end
branch_idx=[branch_idx; cnod];                                             % store root node