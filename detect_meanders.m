function idx_out=detect_meanders(wave,conn,peak_row,peak_col)
% Detects meanders from power image and scale space tree
%
%   MEANDER_ID = det_meand_min(WAVE, CONN, PEAK_ROW, PEAK_COL) Given the
%   spectrum power WAVE (from Torrence and Compo), the scale space tree 
%   specified in CONN and the rows and columns of power peaks in PEAK_ROW 
%   and PEAK_COL, respectively, the function returns the indices of the 
%   nodes in CONN identified as meanders.
%
%   This function implements parts of Section 2.4 of Vermeulen, et al. 2016
%   
%   See also: scale_space_tree, find_peak_in_pol
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

pwr=wave.^2;
% Check input
check_conn(conn)
assert(ismatrix(pwr) && isnumeric(pwr) && all(isfinite(pwr(:))),...
    'POWER should be a numeric, finite, 2D matrix.')

assert(isnumeric(peak_row) && isnumeric(peak_col) &&...
    isequal(size(peak_row),size(peak_col)) &&...
    all(isnan(peak_row) | mod(peak_row,1)==0) &&...
    all(isnan(peak_col) | mod(peak_col,1)==0) &&...
    all(isnan(peak_row) | peak_row <= size(pwr,1)) &&...
    all(isnan(peak_col) | peak_col <= size(pwr,2)),...
    'Invalid PEAK_ROW or PEAK_COL')
    

ml_out=[];                                                                 % vector to contain previously found multiple loop nodes
ml_iter=true;                                                              % True unless no meander nodes are found in branches
while ml_iter                                                              % While a meander node is ancestor of anothoer meander node
    idx_out=[];                                                            % vector to contain ids of meander nodes
    ml_iter=false;                                                         % set iteration control to false (i.e. do not iterate unless a meander node is found in branch, which is checked later)
    for cp=1:size(conn,1)                                                  % loop over all nodes in the tree
        if conn(cp)==0, continue, end                                      % Skip root nodes
        if any(conn==cp), continue, end                                    % Skip non-leaf nodes
        branch=get_branch(cp,conn);                                        % Get current branch
        branch=setdiff(branch,ml_out);                                     % Remove multiple loop nodes from branch (Note that this function sorts the branch putting the root on top!)
        if isempty(branch), continue,end                                   % if nothing remains, continue to next branch
        br_pwr=pwr(sub2ind(size(pwr),peak_row(branch),peak_col(branch)));  % Get power on branch nodes    
        [~,br_m_id]=max(br_pwr);                                           % Get index of meander node (i.e. node with highest power) in branch 
        max_idx=branch(br_m_id);                                           % Get index of meander node (i.e. node with highest power) in tree  
        if ~isempty(max_idx) && ~any(idx_out==max_idx)                     % If we didn't know this meander node (since one meander node can belong to several branches)
            idx_out=[idx_out; max_idx];                        %#ok<AGROW>   store it
        end
        if br_m_id>1                                                       % if meander node is not the root of the current branch (i.e. it has parents)
            new_ml=branch(1:br_m_id-1);                                    % get all parents of the meander node in the current branch (they all are multiple loop nodes)
            ml_out=[ml_out; new_ml];                           %#ok<AGROW> % store all detected multiple loop nodes
            ml_iter=true;                                                  % Iterate!
        end
    end
end
return


    


