function nc=nchild(conn)
% Counts direct children of current node
%   
%   NC = nchild( CONN ) return the number of direct children for each node
%   in CONN. CONN is a pointer to parent vector, typically obtained from
%   scale_space_tree.
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

%% Check input
check_conn(conn);

%% Count children for each node in tree
nc=accumarray(conn+1,conn,[numel(conn)+1 1],@numel,0);                     % count children for each node
nc(1)=[];                                                                  % remove root position
