% Example script to use the meanderscribe functions
%
%   This script shows how to use the functions of the meanderscribe package
%   that implements the method described in:
%
%   Vermeulen, B., A. J. F. Hoitink, G. Zolezzi, J. D. Abad, and R. Aalto 
%       (2016), Multi-scale structure of meanders, Geophys. Res. Lett., 43,
%       doi:10.1002/2016GL068238.
%
%   Before you run this script, make sure you download the Torrence and
%   Compo wavelet package, and point the relevant line in the script below
%   to the right location.


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


%% cleanup
restoredefaultpath                                                         % start with a clean path
clearvars                                                                  % start with a clean workspace
close all                                                                  % close all figures

%% load packages
addpath ../../wavelets                                                     % Load Torrence and Compo 
                                                                           % point this to location where you installed the Torrence and Compo package
                                                                           % available at: http://paos.colorado.edu/research/wavelets/)
addpath ../                                                                % Load Meanderscribe
                                                                           % point this to the location where you placed the meanderscribe functions
%% load data
load mahak_example                                                         % Load Mahakam River data

%% Compute wavelet (section 2.2 in paper)
ds=nanmean(diff(S));                                                       % compute ds

% wavelet parameters
pad = 1;                                                                   % pad the time series with zeroes
dj = 5e-2;                                                                 % spacing of the scales
s0 = -1;                                                                   % default 2*ds 
j1 = -1;                                                                   % use default
mother = 'Dog';                                                            % select mother wavelet: derivative of gaussian
m=2;                                                                       % Wavelet parameter (2nd order)
% Compute wavelet spectrum
[wave,period,scale,coi] = wavelet(C,ds,pad,dj,s0,j1,mother,m);

%% plot the power spectrum
fig_spect=figure('Name','Wavelet specturm and trees','NumberTitle','off'); % Make figure
imagesc(log2(wave.^2));                                                    % plot the spectrum
axis xy                                                                    % use normal axis
colormap(gray);                                                            % make it gray
colorbar                                                                   % add colorbar

%% make the scale space tree (Vermeulen et al. 2016, Section 2.3)
[conn, regions, poly, zc_lines, zc_sign]=scale_space_tree(wave);           % construct the scale space tree 

%% Plot the raw ternary scale_space tree as is done by Witkin 1984
fig_witkin=figure('Name','Scale space ternary tree','NumberTitle','off');  % Make new figure
plot_witkin(regions,'k');                                                  % Plot the tree
xreg=nanmean(regions(:,[3,4]),2);                                          % compute the center of the 2D segments
yreg=nanmean(regions(:,[1,2]),2);                                          % compute the center of the 2D segments
plottree(conn,xreg,yreg,'.-','color',[.5 .5 .5])                           % plot the resulting tree

%% Draw regions and tree over the wavelet power spectrum
set(groot,'CurrentFigure', fig_spect)                                      % Activate right figure
hold on                                                                    % Add to the figure what is coming
for czc=1:numel(zc_lines)                                                  % Loop over all zero crossing lines
    if zc_sign(czc)==-1, col='r'; else col='g'; end                        % Define the color based on the direction of sign change
    plot(zc_lines{czc}(:,2),zc_lines{czc}(:,1),col,'linewidth',2)          % Plot the zero crossing line
end
for cpol=1:numel(poly)                                                     % Loop over all polygons
    npts=(size(poly{cpol},1)+1)/2;                                         % Here we plot only half of each polygon, since all polygons overlap
    plot(poly{cpol}(1:npts,2),poly{cpol}(1:npts,1),'k--')                  % Plot the polygons (regions in the spectrum corresponding to tree nodes)
end

%% Detect meanders (Vermeulen et al. 2016, Section 2.4)
[peak_row,peak_col,peak_pwr]=find_peak_in_pol(poly,wave);                  %  detect peaks in spectrum within each 
conn=remove_nodes(conn,find(isnan(peak_pwr)));                             % remove nodes without a spectrum peak
meander_id=detect_meanders(wave, conn, peak_row, peak_col);                % detect meanders (see function for details)

%% make multiple loop tree Section 2.4
set(groot,'CurrentFigure', fig_spect)
ml_tree=clean_mltree(conn, meander_id);                                    % Clean the tree, resultin in the multiple loop tree (Section 

%% Plot the multiple loop tree on the spectrum
plottree(ml_tree,peak_col,peak_row,'.-','color',[.1 .1 .1])                % plot the resulting tree


%% plot the multiple loop tree on the planform (Bottom panel figure 5)
[xi,yi]=get_centers(ml_tree, peak_row,peak_col,period,ds,X,Y);             % Compute centers of curving sections
figure('Name','Multiple loop tree on the planform','NumberTitle','off');   % Make a figure
plottree(ml_tree,xi,yi)                                                    % Plot the tree
axis equal                                                                 % Make axis scaling equal
hold on,plot(X,Y,'k')                                                      % Plot planform

%% Shape of meanders (Section 2.5)
bounds=meander_bounds(poly,meander_id,peak_row);                           % get the meander boundaries on the planform
[sk, fat]=meander_shape(wave,period, peak_row ,peak_col,...
    meander_id,bounds,scale,ds);                                           % compute shape parameters of meanders

%% Plot the shape parameters
figure('Name','Shape parameters','NumberTitle','off');
plot(S(peak_col(meander_id))*W/1000,sk,'s',S(peak_col(meander_id))*W/1000,fat,'^')
xlabel('Along channel distance (km)')
legend('Skewness (-)','Fattening (-)')