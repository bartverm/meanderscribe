function varargout=plottree(conn,x,y,varargin)
% Plot a tree given the locations of the nodes
%
%   plottree(CONN, X, Y) plots the tree specified in CONN, placing the
%   nodes on the locations defined in X and Y. CONN is a vector with
%   pointers to parent nodes, typically generated with scale_space_tree. X
%   and Y have the same size as CONN and specify the coordinates of the
%   nodes of the tree. CONN is typically generated with the function
%   scale_space_tree
%
%   plottree(CONN, X, Y, ...) any additional argument are passed to the
%   plot function
%
%   h=plottree(...) returns the handle to the object returned by the plot
%   function
%
%   See also: scale_space_tree, plot

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

% Check input
check_conn(conn);                                                          % Check connectivity vector
assert(isequal(numel(conn),numel(x),numel(y)),...
    'CONN, X and Y should have the same number of elements');              % Make sure X and Y have enough elements
assert(isnumeric(x) && isnumeric(y),'X and Y should be numeric');          % Make sure X and Y are numeric

% Plot
cur_hold=get(gca,'NextPlot');                                              % Store hold status of current axes
nc=nchild(conn);                                                           % Compute number of children for each node
h=nan(size(conn));                                                         % Initialize variable to hold handles
for cp=1:numel(conn)                                                       % loop over all nodes
    if cp==2, hold on, end                                                 %    Enable hold at second node
    if conn(cp)==0                                                         %    If node is a root note
        if nc(cp)==0                                                       %        If root node has no children
            h(cp)=plot(x(cp),y(cp),varargin{:});                           %            Plot the root node
        else                                                               %        Otherwise (root node has children)
            continue                                                       %            Do nothing
        end                                                                %        End if
    else                                                                   %    Otherwise (node is not a root node)
        h(cp)=plot(x([cp conn(cp)]), y([cp conn(cp)]),varargin{:});        %        Plot connection
    end                                                                    %    End if
end                                                                        % End loop

set(gca,'NextPlot',cur_hold)                                               % Restore hold status of current axes

if nargout>0                                                               % If more than zero output arguments
    varargout{1}=h;                                                        %    return handles
end                                                                        % End if