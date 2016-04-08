function [xi, yi]=get_centers(conn,peak_row,peak_col,period,ds,x,y)
% Computes the center of curving sections in the multiple loop tree
%
% [XI,YI] = get_centers(CONN, PEAK_ROW, PEAK_COL, PERIOD, DS, X, Y) 
%   computes a coordinate near the center of the curving section in the 
%   tree given in CONN (from scale_space_tree), having period at ROW-th 
%   PERIOD, at the location given in COL. DS is the sampling distance of 
%   the curve given in X,Y.
%
%   See also: scale_space_tree, plottree
%
%   This function can aid to produce a figure similar to the lowest panel
%   of Figure 5 in Vermeulen, et al. (2016)
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

check_conn(conn);
assert(isnumeric(peak_row) && isnumeric(peak_col) &&...
    isequal(size(conn),size(peak_row),size(peak_col)) && ...
    all(isnan(peak_row) | mod(peak_row,1)==0) &&...
    all(isnan(peak_col) | mod(peak_col,1)==0) &&...
    all(isnan(peak_row) | peak_row <= numel(period)) &&...
    all(isnan(peak_col) | peak_col <= numel(x)),'Wrong ROW and/or COL')

assert(isnumeric(period) && isvector(period) && all(isfinite(period)),'Wrong PERIOD')

assert(isnumeric(ds) && isscalar(ds) && isfinite(ds) && ds>0,'Wrong DS')

assert(isnumeric(x) && isnumeric(y) && isequal(size(x),size(y)),'Wrong X and/or Y')

       

[xi,yi]=deal(nan(size(conn,1),1));                                         % Initialize output
fgood=find(all(isfinite([peak_row,peak_col]),2));                                    % Get all tree nodes with power peak
for cn=fgood'                                                              % Loop over all valid nodes
    per=period(peak_row(cn))/2;                                                 % Get period of half-meander
    span=round(per/4/ds);                                                  % Compute span, over which circumcenter is computed (only mid-part, since edges can deviate from overall curving)
    beg_idx=max(1,peak_col(cn)-span);                                           % Index of starting point
    end_idx=min(numel(x),peak_col(cn)+span);                                    % Index of mid point
    mid_idx=round((end_idx+beg_idx)/2);                                    % Index of end point
    x1=x(beg_idx);                                                         % Get x coordinate of starting point
    y1=y(beg_idx);                                                         % Get y coordinate of starting point
    x2=x(mid_idx);                                                         % Get x coordinate of mid point
    y2=y(mid_idx);                                                         % Get y coordinate of mid point
    x3=x(end_idx);                                                         % Get x coordinate of end point
    y3=y(end_idx);                                                         % Get y coordinate of end point
    t=triangulation([1 2 3],[x1;x2;x3],[y1;y2;y3]);                        % Make triangulation of three points
    cc=t.circumcenter(1);                                                  % Compute circumcenter of three points
    xi(cn)=cc(1);                                                          % Assign x coordinate of circumcenter to xi
    yi(cn)=cc(2);                                                          % Assign y coordinate of circumcenter to yi
    l=sqrt((x2-xi(cn))^2+(y2-yi(cn))^2);                                   % Compute distance from mid-point to circumcenter
    rvec=[xi(cn)-x2; yi(cn)-y2]./l;                                        % Make unit vector from mid-point to circumcenter
    xi(cn)=x2+rvec(1)*per/pi;                                              % Compute xi coordinate at a radius of the half period distance along the vector pointing to the incenter
    yi(cn)=y2+rvec(2)*per/pi;                                              % Compute yi coordinate at a radius of the half period distance along the vector pointing to the incenter
end