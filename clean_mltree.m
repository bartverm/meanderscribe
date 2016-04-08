function conn=clean_mltree(conn, meander_id)
% Removes all nodes from the tree with scale smaller than the meander scale
%
%   CONN = clean_mltree( CONN, MEANDER_ID ) removes all children of the
%   meander nodes from the tree specified in CONN (from scale_space_tree). 
%   Meander nodes are identified with MEANDER_ID (from detect_meanders)
%   containing indices of the meander nodes in CONN.
%
%   This function implements part of Section 2.4 in Vermeulen et al. (2016)
%
%   See also: scale_space_tree, detect_meanders
%
%   References:
%   Vermeulen, B., A. J. F. Hoitink, G. Zolezzi, J. D. Abad, and R. Aalto 
%       (2016), Multi-scale structure of meanders, Geophys. Res. Lett., 43,
%       doi:10.1002/2016GL068238.

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

%% Input checking
check_conn(conn)
assert(isnumeric(meander_id) &&...
    all(isfinite(meander_id)) &&...
    all(meander_id>0) &&...
    all(mod(meander_id,1)==0) &&...
    all(meander_id<=numel(conn)),'Invalid meander id''s')

%% Identify nodes that are children of meander nodes
frm=[];
for cp=1:size(conn,1)
    if conn(cp)==0, continue, end  % Skip root nodes
    if any(conn==cp), continue, end % Skip if not lowest level node
    branch=flipud(get_branch(cp,conn)); % Get branch (order from root to lowest level)
    mnode=intersect(branch,meander_id); % Look for meander nodes in branch  
    if ~isempty(mnode) % if there is a meander node in branch (can only be one!)
        mnode_idx=find(branch==mnode); % id in branch of meander node
        frm=[frm; setdiff(branch(mnode_idx+1:end),frm)]; %#ok<AGROW> % mark all nodes below meander node to be removed
    end
end

%% Remove them from the tree
conn=remove_nodes(conn,frm);

