function [divDist] = empDivDist(spikeTrainsArr, divHandle, divParams, indexIterator, verbose)
% Computes the empirical distribution from data collected under a hypothesis
% divDist = empDivDist(spikeTrainsArr, divHandle, divParams, indexIterator)
%
% Spike trains from a hypothesis is supplied along with a divergence measure.
% This function will compute the pairwise divergence values and return
% the structure related to the empirical distribution.
%
% This function is EXPERIMENTAL.
% 
% Input:
%   spikeTrainsArr: (array of spikeTrains structures)
%                   pairs of spikeTrains from this array will be selected via
%                   indexIterator function to compute the divergence dist
%   divHandle: (@) one of the divergence functions for spike trains
%   divParams: (struct) the corresponding parameters for the divHandle
%   indexIterator: (@: (k1, k2) -> (k1, k2)) iterates through the possible
%                  pairs of spike trains to compute the divergence while
%                  maintaining the hypothesis. Given the last indices,
%                  the next indices should be returned. Loop starts with (1,1)
%                  Default computes all possible pairwise interactions.
%   verbose: (default: true) print calculation time estimates
% Output:
%   divDist: (struct) empirical computation results and related function handles
%
% See also: div*
%
% $Id$
% Copyright 2009 iocane project. All rights reserved.

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

M = length(spikeTrainsArr);
nSurro = M * (M - 1) / 2;
d_pairwise = zeros(nSurro, 1);

if nargin > 4
    next = indexIterator;
else
    next = @defaultNext;
end

if nargin < 5
    verbose = true;
end

k = 1; k1 = 1; k2 = 1;
tt = 0;
while true
    if verbose; t = tic; end
    [k1, k2] = next(k1, k2);
    if k1 < 0
	break;
    end

    d_pairwise(k) = divHandle(spikeTrainsArr(k1), spikeTrainsArr(k2),divParams);
    if verbose; tt = tt + toc(t);
	waitMins = (tt / k) * (nSurro-k) / 60;
	if waitMins > 1
	    fprintf('Estimated time remaining (surrogate): %f mins [%d]\r', ...
		waitMins, k);
	end
    end
    k = k + 1;
end

if any(isnan(d_pairwise))
    warning('There are NaN values in the divergence!');
end

if any(isinf(d_pairwise))
    warning('There are Inf values in the divergence!');
end

divDist.unsorted = d_pairwise;
divDist.values = sort(d_pairwise);
divDist.N = length(d_pairwise);
divDist.divHandle = divHandle;
divDist.divParams = divParams;
divDist.method = mfilename;
divDist.pValue = @(alpha)(sum(divDist.values >= alpha) / divDist.N);
divDist.threshold = @(p)(divDist.values(ceil((1 - p) * divDist.N)));

    function [k1, k2] = defaultNext(k1, k2)
    % independent spike train structures
	if k2 == M
	    k1 = k1 + 1;
	    k2 = k1 + 1;
	else
	    k2 = k2 + 1;
	end

	if k1 == M
	    k1 = -1;
	end
    end

end
% vim:ts=8:sts=4:sw=4
