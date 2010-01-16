function [str] = dparams2str(params)
% Create short(?) description string for div*Param objects
% str = dparams2str(params)
% 
% Input:
%   params: any data type parameter object to be converted (cell not supported)
% Output:
%   str: the short string representation
%
% See also: div*Params
%
% $Id$
% Copyright 2010 iocane project. All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%  - Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  - Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%  - Neither the name of the iocane project nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

str = '';
if isempty(params)
    return;
end

if iscell(params)
    warning('Cell type not supported.');
    return;
end

if isa(params, 'function_handle')
    str = func2str(params);
    return;
end

if isstruct(params)
    if isfield(params, 'desc')
	str = params.desc;
    else
	str = '(';
	fns = fieldnames(params);
	for k = 1:length(fns)
	    val = getfield(params, fns{k});
	    str = [str sprintf('%s=%s', fns{k}, dparams2str(val))];
	    if k ~= length(fns)
		str = [str ','];
	    end
	end
	str = [str ')'];
    end
    return;
end

if isnumeric(params)
    if numel(params) == 1
	str = num2str(params);
	return;
    end
    % Now it's a array of something
    str = '[';
    for k = 1:numel(params)
	str = [str sprintf('%s', dparams2str(params(k)))];
	if k ~= numel(params)
	    str = [str ','];
	end
    end
    str = [str ']'];
    return;
end

if ischar(params) || isstring(params)
    str = params;
    return;
end
