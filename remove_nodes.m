function conn=remove_nodes(conn,frm)
% Remove nodes from a directed tree (they become isolated root nodes)
%
%   CONN = remove_nodes(CONN, FRM) removes the nodes given in FRM from the
%   tree specified in CONN. This tree is typically generated with
%   scale_space_tree. The child nodes of the removed nodes are connected to
%   the parent of the removed nodes. The removed nodes are made root nodes.
%   They end up being root nodes with no children. 
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

% Check input
check_conn(conn);
assert(isnumeric(frm) &&...
    all(isfinite(frm)) &&...
    all(frm>0) &&...
    all(mod(frm,1)==0) &&...
    all(frm<=numel(conn)),...
    'FRM should specify the elements in CONN to be removed');

% remove nodes
for cn=1:numel(frm)
    conn(conn==frm(cn),1)=conn(frm(cn),1);                                 % set new connection for child of current node
    conn(frm(cn),1)=0;                                                     % make current node a root node
end