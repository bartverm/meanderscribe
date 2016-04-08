function bounds=meander_bounds(poly,meander_id,peak_row)
% Computes the start and end-point of a meander
%
% BOUNDS = meander_bounds( POLY, MEANDER_ID, PEAK_ROW) computes the start 
%   and end-point of a meander based on the bounding zero-crossing lines at
%   the period where the meander wavelet power peak occurs. MEANDER_ID 
%   (from detect_meanders) contains the indices of the tree nodes that have
%   been detected as meanders. PEAK_ROW (from find_peak_in_pol) contains 
%   the row indicating  the period at which the highes power peak occurs in
%   the tree. POLY (from scale_space_tree) is a cell holding the boundaries
%   (including the zero crossing lines) of the nodes in the tree. BOUNDS 
%   has the same number of rows as MID, holding in the first column the 
%   starting indices of the meanders and in the last column the ending 
%   indices of the meanders.
%
%   See also: detect_meanders, meander_shape, scale_space_tree,
%       find_peak_in_pol
%
%   This function implements parts of Section 2.5 in Vermeulen, et al. 2016
%   (sentence starting with: Given the meander node, ....)
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

assert(iscell(poly) && all(cellfun(@(x) ismatrix(x) && size(x,2)==2 &&...
    isnumeric(x), poly)),'Wrong POLY')                                     % Check validity of POLY

assert(isnumeric(meander_id) && all(mod(meander_id,1)==0) &&...
    all(meander_id>0) &&...
    all(meander_id<=numel(poly)),'Wrong MID')                              % Check validity of MEANDER_ID

assert(isnumeric(peak_row) && all(isnan(peak_row) |...
    mod(peak_row,1)==0),'Wrong ROW')                                       % Check validity of PEAK_ROW

m_row=peak_row(meander_id);                                                % Get period index of meanders
m_pol=poly(meander_id);                                                    % Get scale-space region polygons of the meander nodes
bounds=nan(numel(m_row),2);                                                % Initialize output
for c_pol=1:numel(m_row)                                                   % Loop over all meanders
    bounds(c_pol,:)=unique(m_pol{c_pol}(m_pol{c_pol}(:,1)==...
         m_row(c_pol),2));                                                 % Get the polygon coordinate at the meander peak period (these coordinate are on the bounding zero-crossing lines
end

