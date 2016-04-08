function lat=hexl(inmat)
% Constructs a hexagonal lattice for saddle points and extremes detection
%   
%   LAT = hexl( MAT ) given the 2D matrix MAT it returns the hexagonal
%   lattice of the matrix around each point in LAT. The size of LAT is
%   [ size(MAT), 6 ], i.e. it has six elements in the third dimension each
%   containing one of the six surrounding elements in the lattice, in the
%   following order: top left, top center, mid left, mid right, bottom
%   center, bottom right.
%
%   This functions was written to be used by sedext.m which detects saddle
%   points and extremes based on the method of Kuijper (2004).
%   
%  See also: sadext
%
%  References:
%  Kuijper, A. (2004), On detecting all saddle points in 2D images, Pattern
%       Recogn. Lett., 25 (15), 1665-1672, doi:10.1016/j.patrec.2004.06.017

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
assert(ismatrix(inmat), 'Input should be a 2D matrix');                    % INMAT should be a matrix
[nr,nc]=size(inmat);                                                       % Get size of input
assert(nr>2 && nc>2, 'Matrix too small to build a hexagonal lattice');     % Make sure it's large enough to make the hexagonal lattice

%% Make the lattice
[lat1,lat2,lat3,lat4,lat5,lat6]=deal(inmat);                               % initialize variables with the input matrix
lat1([nr-1 nr],:)=[]; lat1(:,[nc-1 nc])=[];                                % top left
lat2([1    nr],:)=[]; lat2(:,[nc-1 nc])=[];                                % top center
lat3([nr-1 nr],:)=[]; lat3(:,[1    nc])=[];                                % mid left
lat4([1     2],:)=[]; lat4(:,[1    nc])=[];                                % mid right
lat5([1    nr],:)=[]; lat5(:,[1     2])=[];                                % bottom center
lat6([1     2],:)=[]; lat6(:,[1     2])=[];                                % bottom right
lat=cat(3,lat1,lat2,lat4,lat6,lat5,lat3);                                  % make lattice
