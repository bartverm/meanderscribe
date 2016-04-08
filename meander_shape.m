function [skval,fatval]=meander_shape(wave, period, peak_row, peak_col, meander_id, bounds, scale, ds)
% Computes the meander shape parameters fattening and skewing
%
%   [SK, FAT] = meander_shape( WAVE, PERIOD, PEAK_ROW, PEAK_COL, MEANDER_ID
%       , BOUNDS, SCALE, DS) computes the meander shape parameters for all
%       meanders. WAVE is the wavelet spectrum (from Torrence and Compo), 
%       PERIOD contains the periods at which the wavelet spectrum is 
%       computed (from Torrence and Compo), PEAK_ROW and PEAK_COL
%       define the location of the spectral peak corresponding to the
%       meanders (from find_peak_in_pol), MEANDER_ID is the meander index 
%       (from detect_meanders), BOUNDS contains the spatial boundaries of 
%       the meanders (from meander_bounds), SCALE contains the scaling of 
%       the wavelet spectrum for each period (from Torrence and Compo), DS 
%       is the sampling distance of the meandering planform.
%
%   See also: find_peak_in_pol, detect_meanders, meander_bounds
%
%   This function implements parts of Section 2.5 in Vermeulen, et al. 2016
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

assert(isnumeric(wave) && ismatrix(wave), 'Wrong WAVE');
assert(isnumeric(period) && numel(period) == size(wave,1), 'Wrong PERIOD');
assert(isnumeric(peak_row) && isnumeric(peak_col) &&...
    isequal(size(peak_row),size(peak_col)) &&...
    all(isnan(peak_row) | mod(peak_row,1)==0) &&...
    all(isnan(peak_col) | mod(peak_col,1)==0) &&...
    all(isnan(peak_row) | peak_row <= size(wave,1)) &&...
    all(isnan(peak_col) | peak_col <= size(wave,2)),...
    'Invalid PEAK_ROW or PEAK_COL')
assert(isnumeric(meander_id) && all(mod(meander_id,1)==0) &&...
    all(meander_id>0) &&...
    all(meander_id<=numel(peak_row)),'Wrong MID')                          % Check validity of MEANDER_ID
assert(isnumeric(bounds) && size(bounds,1)==numel(meander_id) &&...
    size(bounds,2)==2 && all(mod(bounds(:),1)==0) &&...
    all(bounds(:)>0) && all(bounds(:)<=size(wave,2)),'Wrong BOUNDS')       % Check validity of BOUNDS
assert(isnumeric(scale) && numel(scale)==size(wave,1),'Wrong SCALE')       % Check validity of SCALE
assert(isnumeric(ds) && isscalar(ds) && ds>0,'Wrong DS')                   % Check validity of DS


% find primary and secondary wave
[skval, fatval]=deal(nan(numel(meander_id),1));                            % initialize output variables
for cm=1:numel(meander_id)                                                 % loop over all half meanders
    mfilt=bounds(cm,1):bounds(cm,2);                                       % filter to select current half meander from the whole river
    mpfilt=peak_row(meander_id(cm));                                       % filter to select scale of current meander
    phi=(mfilt-mfilt(1))./(mfilt(end)-mfilt(1))*pi;                        % compute phi
    mper=period(mpfilt);                                                   % Meander scale
    spfilt=diff((period-mper/3)>0)==1;                                     % Filter to select scale for shape parameters (1/3 of meander scale)
    sper=period(spfilt);                                                   % Scale of shape
    if isempty(sper)                                                       % if we don't have the shape-scale in the spectrum (meander too small)
        continue                                                           % skip current meander
    end
    mwave=wave(mpfilt,mfilt)*sign(wave(peak_row(meander_id(cm)),...
        peak_col(meander_id(cm))));                                        % spectrum at meander scale
    swave=wave(spfilt,mfilt)*sign(wave(peak_row(meander_id(cm)),...
        peak_col(meander_id(cm))));                                        % spectrum at shape scale
    
    % detect secondary peaks
    speak_idx=find(diff(diff(swave)>0)==-1)+1;                             % locate were derivative of shape spectrum goes from positive to negative (maxima)
    speak_max=max(swave(speak_idx));                                       % Find strongest maximum in shape spectrum
    if isempty(speak_max), speak_max=0; end                                % If there is no peak in shape spectrum, set peak-strength to zero
    mpeak_max=max(mwave);                                                  % Get the strength of the meander peak (always exist and corresponds with absolute maximum)
    scaling=sqrt(ds./scale(spfilt));                                       % sampling frequency factor
    skv=cos(3*phi);                                                        % compute skewness function
    fatv=sin(3*phi);                                                       % compute fattening function
    skval(cm)=scaling*sum(swave.*skv)./sum(skv.^2)*speak_max./mpeak_max;   % Compute skewness parameter
    fatval(cm)=scaling*sum(swave.*fatv)./sum(fatv.^2)*speak_max./mpeak_max;% Compute fattening parameter
end


