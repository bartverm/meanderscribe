function id=sadext(inmat)
% Detect saddle points and extremes in a 2D image (Kuijper, 2004)
%   
%  ID = sadext( MAT ) detecte saddle points and extremes in the image given
%  in MAT. ID has the same size as MAT and contains an id which can have
%  one of the following values:
%
%  0: The point is a local extreme
%  2: The point is a regular point
%  4: The point is a saddle point
%  6: The point is a degenerate saddle point
%
%  This function requires the function hexl.m to build the hexagonal
%  lattice.
%
%  See also: hexl
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

assert(ismatrix(inmat) && isnumeric(inmat),...
    'Input should be a numerical 2D matrix');                              % input checking

% create hexagonal lattice
inmat6=hexl(inmat);                                                        % Create hexagonal lattice (see function for details)
ccell=inmat;                                                               % copy input
ccell([1 end],:)=[];                                                       % remove first and last row
ccell(:,[1 end])=[];                                                       % remove first and last column
id=sum(abs(diff(bsxfun(@gt,cat(3,inmat6,inmat6(:,:,1)),ccell),1,3)),3)/2;  % compute id 
id(any(isnan(cat(3,inmat6,ccell)),3))=nan;                                 % if any of the elements and the surrounding lattice is nan, make the id nan
id=id([1 1:end end],:);                                                    % duplicate first and last column
id=id(:,[1 1:end end]);                                                    % duplicate first and last row
