function varargout=scale_space_tree(wave)
% Construct a scale space tree from the zero crossings of the wavelet
%
%   CONN = scale_space_tree( WAVE ) returns the ternary scale space tree
%   generated based on the zero crossings of the wavelet spectrum. 
%   WAVE is a scale space image with each row representing a scale and the
%   column a position in space or time. The topmost row correpsonds to the
%   smallest scale. This matrix is typically made with the Torrence and
%   Compo package.
%   CONN is a column vector in which each row corresponds to a tree node. 
%   Root nodes have a value of zero. All other nodes have a value pointing 
%   to the row in CONN of their parent node
%
%   [CONN, REGIONS] = scale_space_tree( WAVE ) returns the matrix REGIONS.
%   This matrix contains the same number of rows as CONN. For each row in
%   CONN the corresponding row in REGIONS gives the scale boundaries of the
%   node and the spatial boundaries derived from where the zero crossing
%   lines leave the scale space plane. The first column is the smallest
%   bounding period, the second column is the larges bounding period, the
%   third column is the spatial coordinate where the region starts 
%   and the last column is the spatial coordinate where the region ends.
%   Use the function plot_witkin to plot this matrix.
%
%  [CONN, REGIONS, POLY] = scale_space_tree( WAVE ) returns a vector of
%  cells with the same size as CONN containing the coordinates of a polygon 
%  bounding the region of the corresponding node in CONN. Unlike REGIONS
%  these coordinates follow the zero crossing lines bounding the region.
%  The coordinates are given as (row, col)
%
%  [CONN, REGIONS, POLY, ZC_LINES] = scale_space_tree( WAVE ) returns all
%  zero crossing lines detected. Each cell element in ZC_LINES contains an
%  Nx2 matrix containing the coordinates of the zero crossing lines given
%  as (row, col).
%
%  [..., ZC_SIGN] = scale_space_tree( WAVE ) returns the sign of each zero
%  crossing line in ZC_LINES
%
%  This functions implements Section 2.3 of Vermeulen, et al. (2016)
%   
%  See also: detect_meanders, find_peak_in_pol, plot_witkin
%
%  References:
%  Vermeulen, B., A. J. F. Hoitink, G. Zolezzi, J. D. Abad, and R. Aalto 
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

%% Prepare input
wave = flipud(wave);                                                       % put largest scale on top of the matrix (easier when processing from large to small scales)

%% Find zero crossing locations in spectrum (function changes sign), which mark the edges in the s-direction of the 2d-segments
zcr=[diff(wave>0,1,2), wave(:,1)*0];                                       % Find zero crossing (where wavelet spectrum changes sign). This gives -1 when going from positive to negative, 1 for negative to positive, 0 otherwise

% To properly bound entire spectrum we add zero crossings at the right and
% left edge of the spectrum (i.e maximum and minimum s). Since positive to
% negative (-1) and negative to positive (1) transitions always alternate,
% we add -1 if nearest zero crossing was 1 and vice-versa

for cr=1:size(zcr,1)                                                       % for all scales (rows)
    % left side
    if zcr(cr,1)==0                                                        % if there is no zero-crossing at the beginning of the spectrum (s==0)
        zcr(cr,1)=-zcr(cr,find(zcr(cr,:)~=0,1,'first'));                   % Add crossing with opposite sign of nearest one
    end
    % right side
    if zcr(cr,end)==0                                                      % if there is no zero-crossing at the end of the spectrum (s==max(s))
        zcr(cr,end)=-zcr(cr,find(zcr(cr,:)~=0,1,'last'));                  % Add zero crossing with opposite sign of nearest one
    end 
end

% At this stage we detected all zero crossing points, but we still did not
% connect them! We still need to connect the points to get the
% zero-crossing lines.

%% Find singular points (where new zero-crossing lines emerge), which mark the edges in scale-direction of the 2d segments
nzeros=sum(abs(zcr),2);                                                    % Total number of zero_crossing lines emerging at each row
sing_pnts_rw=find([1; diff(nzeros)>0]);                                    % Find row (period) where singular points occur (where number of zero crossing changes). Note that first row/period is considered as singular one, such that it marks the start of a segment

% compute total number of nodes for variale initialization
nnew_zc_lines_lateral=sum(abs([0 0;diff(zcr(:,[1 end]),1,1)])/2,2);        % Number of new zero_crossing lines emerging from lateral singular points at each row
nnew_zc_lines=[0; diff(nzeros)];                                           % Total number of zero crossing lines emerging at each row
nnew_zc_lines_internal=nnew_zc_lines-nnew_zc_lines_lateral;                % Number of new zero_crossing lines emerging from internal singular points
nnew_internal_nodes=nnew_zc_lines_internal+nnew_zc_lines_internal/2;       % Number of nodes emerging at internal singular points (3 for each pair of zc lines
nnew_lateral_nodes=nnew_zc_lines_lateral;                                  % Number of nodes emerging at lateral singular points (1 for each zc line)
n_nodes_initial=nzeros(1)-1;                                               % Number of nodes at largest scale
n_nodes=n_nodes_initial+sum(nnew_internal_nodes+nnew_lateral_nodes);       % Total number of nodes in the tree
n_zclines=nzeros(1)+sum(nnew_zc_lines);                                    % Total number of zero crossings in the tree

%% Make zero-crossing lines from zero-crossing locations
% Now we build the tree based on zero-crossing lines and singular points.
ml_idx=zeros(1,size(wave,2));                                              % vector to store 'claimed' sections (same size of s). Starts with zeros, since all first nodes are root nodes.
loc_x=nan(size(wave));                                                     % Matrix to store s location of zero crossings (i.e. the s-location of the zero crossing line in smallest period)

% Initialize variable
[xtr_colmin, xtr_colmax, xtr_rwmin, xtr_rwmax, xtr_conn]=...
    deal(nan(n_nodes,1));                                                  
xtr_poly=cell(n_nodes,1);
zc_lines=cell(n_zclines,1);
zc_sign=nan(n_zclines,1);

% Initialize counters
cc=0; % counts the segments in the tree
last_zc=0; % counts the zero crossing lines

for cs=1:numel(sing_pnts_rw)                                               % loop over all singular periods (i.e. periods where number of zero crossing changes)
    %
    % Determine s-location of zero crossing. In a matrix with the same size
    % as the spectrum this location is stored in each value belonging to a
    % zero-crossing line
    %
    zercr_col=find(zcr(sing_pnts_rw(cs),:)~=0);                            % find the zero crossings at current period
    znew_col=find(isnan(loc_x(sing_pnts_rw(cs),zercr_col)));               % Find the new zero crossings (not processed at already at larger periods, which therefore still don't have an s-location assigned to them)
    col_fin=znew_col*nan;                                                  % make a nan vector to store ending column, i.e. s location where the zero crossing line ends (in smallest periods)
    for cl=1:numel(znew_col)                                               % for each new zero crossing (we will now find the s-location)
        [lrw,lcol]=get_zcline(zcr,sing_pnts_rw(cs),zercr_col(znew_col(cl))); % track the zero crossing line (see function fnd_line for details)
        last_zc=last_zc+1;                                                 % increase zc-line counter
        zc_lines{last_zc}=[size(zcr,1)+1-lrw,lcol];                        % store current zc-line
        zc_sign(last_zc)=zcr(sing_pnts_rw(cs),zercr_col(znew_col(cl)));    % store sign of zc-line
        col_fin(cl)=lcol(end);                                             % store ending column for current zc-line
        loc_x(sub2ind(size(loc_x),lrw,lcol))=lcol(end);                    % New zero crossing lines get ending column assinged to them
    end
    %% Build the tree (Section 2.3 paragraph 3)
    % Build the tree. The important variable to construct is xtr_conn. 
    % It contains the connectivity of the tree. It is a vector in
    % which each element corresponds to a node (which corresponds to a 2d
    % segment) in the tree. Each element contains a number which points to
    % the parent node. If an element contains a zero it is a root-node.
    % Example: [0 1 2 0 2 4] means: first node is a root node, second is
    % connected to first, third to second (this node is a leaf, i.e. no 
    % children), fourth is a root node, fifth is connected to the second 
    % (and is a leaf node), sixth is connected to fourth (and is a leaf 
    % node). All other xtr-variables contain the bounding edges of the 
    % 2d-portion of the plane
    for cn=1:numel(zercr_col)-1                                            % For all zero crossing at current period except the last one (i.e. all left boundaries of segments)
        if any(znew_col==cn | znew_col==cn+1)                              % if current or next node is new
            cc=cc+1;                                                       % we found a new node, thus increment node counter
            xtr_rwmin(cc)=sing_pnts_rw(cs);                                % period where node starts
            xtr_colmin(cc)=loc_x(sub2ind(size(loc_x),sing_pnts_rw(cs),...
                zercr_col(cn)));                                           % left side column (minimum s location) for current segment at smallest period
            xtr_colmax(cc)=loc_x(sub2ind(size(loc_x),sing_pnts_rw(cs),...
                zercr_col(cn+1)));                                         % right side column (maximum s location) for node
            xtr_idx=round((xtr_colmin(cc)+xtr_colmax(cc))/2);              % The column for the node is taken equal to the center between minimum and maximum column 
            xtr_conn(cc)=ml_idx(xtr_idx);                                  % Find parent node
            ml_idx(xtr_colmin(cc):xtr_colmax(cc))=cc;                      % Claim the s-portion for this node
            if xtr_conn(cc)~=0                                             % If parent is not root, make the node boundaries for the parent (now we now where parent ends!)
                xtr_rwmax(xtr_conn(cc))=sing_pnts_rw(cs);                  %period (row) where node/segment ends
                %% Make a polygon bounding the node/segment (we do this for the parent segment since we do not know where current node will end!)
                [lline_rw, lline_col]=find((...                            
                    loc_x==xtr_colmin(xtr_conn(cc))));                     % Get left zero crossing line
                [lline_rw,srti]=sort(lline_rw);                            % Sort point (rows) in ascending period order
                lline_col=lline_col(srti);                                 % Sort point (columns) in ascending period order
                [rline_rw, rline_col]=find((...
                    loc_x==xtr_colmax(xtr_conn(cc))));                     % Get right zero crossing line
                [rline_rw,srti]=sort(rline_rw);                            % Sort point (rows) in ascending period order
                rline_col=rline_col(srti);                                 % Sort point (columns) in ascending period order
                lline_col(lline_rw<xtr_rwmin(xtr_conn(cc))|...
                    lline_rw>xtr_rwmax(xtr_conn(cc)))=[];                  % keep only part of left line column-coordinates (s) between upper and lower bound of segment (minimum and maximum row/period)
                lline_rw(lline_rw<xtr_rwmin(xtr_conn(cc))|...
                    lline_rw>xtr_rwmax(xtr_conn(cc)))=[];                  % keep only part of left line row-coordinates (period) between upper and lower bound of segment (minimum and maximum row/period)
                rline_col(rline_rw<xtr_rwmin(xtr_conn(cc))|...
                    rline_rw>xtr_rwmax(xtr_conn(cc)))=[];                  % keep only part of right line column-coordinates (s) between upper and lower bound of segment (minimum and maximum row/period)
                rline_rw(rline_rw<xtr_rwmin(xtr_conn(cc))|...
                    rline_rw>xtr_rwmax(xtr_conn(cc)))=[];                  % keep only part of right line row-coordinates (period) between upper and lower bound of segment (minimum and maximum row/period)
                xtr_poly{xtr_conn(cc)}=...
                    [size(zcr,1)+1-[lline_rw;flipud(rline_rw);lline_rw(1)],...
                    [lline_col;flipud(rline_col);lline_col(1)]];           % construct the bounding polygon of segment parent to the current one
            end % if xtr_conn(cc)~=0 
        end % if any(znew_col==cn | znew_col==cn+1)
    end % for cn=1:numel(zercr_col)-1  
end % loop over all singular periods (i.e. periods where number of zero crossing changes)

%% Make lower bounds and polygon for leaves (smallest scale segments)
xtr_rwmax(isnan(xtr_rwmax))=size(zcr,1);                                   % add lower boundary for leave nodes found at earlier singular point
fleaf=find(xtr_rwmax==size(zcr,1));                                        % get all leave nodes

% TODO: code below is the same as code above in the for loop, make it an
% internal function
for cc=1:numel(fleaf)
    % make the leaf polygons
    [lline_rw, lline_col]=find((loc_x==xtr_colmin(fleaf(cc))));            % Get left zero crossing line
    [lline_rw,srti]=sort(lline_rw);                                        % Sort point (rows) in ascending period order
    lline_col=lline_col(srti);                                             % Sort point (columns) in ascending period order
    [rline_rw, rline_col]=find((loc_x==xtr_colmax(fleaf(cc))));            % Get right zero crossing line
    [rline_rw,srti]=sort(rline_rw);                                        % Sort point (rows) in ascending period order
    rline_col=rline_col(srti);                                             % Sort point (columns) in ascending period order
    lline_col(lline_rw<xtr_rwmin(fleaf(cc))|...
        lline_rw>xtr_rwmax(fleaf(cc)))=[];                                 % keep only part of left line column-coordinates (s) between upper and lower bound of segment (minimum and maximum row/period)
    lline_rw(lline_rw<xtr_rwmin(fleaf(cc))|...
        lline_rw>xtr_rwmax(fleaf(cc)))=[];                                 % keep only part of left line row-coordinates (period) between upper and lower bound of segment (minimum and maximum row/period)
    rline_col(rline_rw<xtr_rwmin(fleaf(cc))|...
        rline_rw>xtr_rwmax(fleaf(cc)))=[];                                 % keep only part of right line column-coordinates (s) between upper and lower bound of segment (minimum and maximum row/period)
    rline_rw(rline_rw<xtr_rwmin(fleaf(cc))|...
        rline_rw>xtr_rwmax(fleaf(cc)))=[];                                 % keep only part of right line row-coordinates (period) between upper and lower bound of segment (minimum and maximum row/period)
    xtr_poly{fleaf(cc)}=...
        [size(zcr,1)+1-[lline_rw;flipud(rline_rw);lline_rw(1)],...
        [lline_col;flipud(rline_col);lline_col(1)]];                       % construct the bounding polygon of segment parent to the current one
end

%% Generate output as requested
if nargout>0
    varargout{1}=xtr_conn;
end
if nargout > 1
   varargout{2} = [size(zcr,1)+1-xtr_rwmin,...
    size(zcr,1)+1-xtr_rwmax, xtr_colmin, xtr_colmax];
end
if nargout>2
    varargout{3}=xtr_poly;
end
if nargout>3
    varargout{4}=zc_lines;
end
if nargout>4
    varargout{5}=zc_sign;
end

