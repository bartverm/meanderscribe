function [rw_out,col_out]=get_zcline(zcr,rw_start,col_start)
% Find a zero crossing line given its starting position (singular point)
%
%   [RW, COL] = get_zcline(ZCR, RW_START, COL_START) track the zero crossing
%   line staring at (ROW, COL) base on the matrix ZCR. This matrix holds a
%   1 for negative to positive zc line, a -1 for a positive to negative
%   zc-line and a zero otherwise.
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

% The procedure is to find in each next period (from large to small) the
% nearest zero-crossing point with the same sign (i.e. -1 for pos-to-neg
% and 1 for neg-to-pos)

rwend=size(zcr,1);                                                         % Row at which tracking stops
sgn=zcr(rw_start,col_start);                                               % Store sign of zero-crossing we a tracking (pos-to-neg or neg-to-pos)
assert(sgn==1 | sgn==-1);                                                  % Make sure the starting point is on a zero-crossing line

n_out=rwend-rw_start+1;                                                    % Length of zero-crossing line (how many periods does it span?)
col_out=deal(nan(n_out,1));                                                % vector to store columns of zero-crossing line
rw_out=(rw_start:rwend)';                                                  % vector to store rows of zero-crossing line

[rw_ss, col_ss]=find(zcr==sgn);                                            % Find rows and columns of zero-crossing points with same sign as the given starting zero-crossing point

col_out(1)=col_start;                                                      % The first point of the line is the given starting point

% Make the line
for cr=2:n_out                                                             % loop over all periods between the starting point and given ending period
    cur_col=col_ss(rw_ss==rw_out(cr));                                     % Store column number of all (equally signed) zero-crossings at current period
    col_diff=abs(cur_col-col_out(cr-1));                                   % compute distance between previous zc-point and zc-points at current period
    col_out(cr)=cur_col(find(col_diff==min(col_diff),1));                  % keep the nearest zc-point at current period
end                                                                        % continue to next period

