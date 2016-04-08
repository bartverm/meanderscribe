function [peak_row, peak_col, peak_pwr]=find_peak_in_pol(poly,wave)
% Detect peaks of a 2D function within the given simple polygons
%
%   [ROW, COL, PWR]=find_peak_in_pol(POLY, WAVE) Given a 2D matrix WAVE 
%   (from Torrence and Compo), and polygons in POLY (from scale_space_tree)
%   returns the position and power of the strongest local extreme of WAVE 
%   within the polygons given in POLY. POLY is a cell of Nx2  matrices, 
%   each row defining the points of the polygon. The coordinates are given 
%   in row and column number of the matrix WAVE. PEAK_ROW and PEAK_COL are 
%   arrays with the same size as poly containing the row and column number 
%   of the strongest peak occuring within the corresponsing polygon in 
%   POLY. If no extreme is found, nans  are returned. PEAK_PWR contains the
%   power of the peak.
%
%   This function implements parts of Section 2.4 of Veremulen, et al.
%   (2016)
%
%   See also: scale_space_tree, detect_meanders, sadext, hexl
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

power=wave.^2;
se_wave=sadext(wave);                                                     % Detect saddle points and extremes in spectrum (according to Kuijper) (see function for details)
[rw_xtr, col_xtr]=find(se_wave==0);                                        % Find extremes (id = 0)


        
[peak_pwr,peak_row,peak_col]=deal(nan(size(poly)));
s = warning('off','MATLAB:inpolygon:ModelingWorldLower');                  % Disable warning for too small bounding box
for cp=1:numel(poly)                                                       % For each polygon
    finpol=find(inpolygon(rw_xtr,col_xtr,poly{cp}(:,1),poly{cp}(:,2)));    % Search for extreme withing current polygon
    if isempty(finpol)                                                     % If no extreme lies in polygon
        continue                                                           % Continue to next polygon
    else                                                                   % If extremes were found in polygon
        ppeak=power(sub2ind(size(power),rw_xtr(finpol),col_xtr(finpol)));  % Get the power of the peaks
        peak_pwr(cp)=max(ppeak);                                                % Store strongest peak
        fmax=ppeak==peak_pwr(cp);                                               % find index of strongest peak
        peak_row(cp)=rw_xtr(finpol(fmax));                                      % Store row of peak
        peak_col(cp)=col_xtr(finpol(fmax));                                    % Store column of peak
    end
end
warning(s);                                                                % Re-enable warning
