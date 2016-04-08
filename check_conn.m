function check_conn(conn)
% Check the validity of connection vector describing a directed tree
%
% See also: scale_space_tree, detect_meanders, find_peak_in_pol

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

id='check_conn:wrong_conn';

do_throw=false;
cause='';
if ~isvector(conn), cause=[cause, 'CONN shoud be a vector\n'];do_throw=true; end
if ~isnumeric(conn), cause=[cause, 'CONN shoud be numeric\n'];do_throw=true; end
if ~isfinite(conn), cause=[cause, 'CONN shoud be finite\n'];do_throw=true; end
if ~all(conn>=0), cause=[cause, 'CONN shoud hold positive values\n'];do_throw=true; end
if ~all(mod(conn,1)==0), cause=[cause, 'CONN shoud hold integers\n'];do_throw=true; end
if ~all(conn<=numel(conn)), cause=[cause, 'CONN shoud be a vector\n'];do_throw=true; end

if do_throw
    throw(MException(id,cause))
end