function [div] = divSPD(spikeTrains1, spikeTrains2, params)
% Divergence using strictly positive definite kernels
% [div] = divSPD(spikeTrains1, spikeTrains2, params)
%
% Input:
%   spikeTrains1, spikeTrains2: (struct) sets of spike trains for comparison
%   params: (struct) see divSPDParams
% Output:
%   div: (1) divergence value
%
% See also: divSPDParams, divCount, divL2Poisson
%
% Original idea by Memming and Sohan Seth
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

M1 = cellfun('length', spikeTrains1.data);
M2 = cellfun('length', spikeTrains2.data);
maxM = max(max(M1),max(M2));

% sigma = 1e-3;
% kernel = @(x,y,N)(exp(-sum((x-y).^2)/2/sigma/sqrt(N)));
% kernel = @(x,y,N)(exp(-sum((x-y).^2)));
% kernel = @(x,y,N)(1);
kernel = params.kernel;

% N = 0 case
k1 = (M1 == 0);
k2 = (M2 == 0);

dxx = sum(k1)^2; % TODO: multiply by K(0)?
dyy = sum(k2)^2;
dxy = sum(k1) * sum(k2);

for N = 1:maxM % for each dimension
    k1 = (M1 == N); k2 = (M2 == N); % find all the spike trains with N spks

    NX = sum(k1); NY = sum(k2);

    X = spikeTrains1.data(k1);
    Y = spikeTrains2.data(k2);

    if NX ~= 0
	tdxx = 0;
	for k1 = 1:NX
	    for k2 = (k1+1):NX
		tdxx = tdxx + kernel(X{k1}, X{k2}, N);
	    end
	end
	dxx = dxx + tdxx * 2 + NX;
    end

    if NY ~= 0
	tdyy = 0;
	for k1 = 1:NY
	    for k2 = (k1+1):NY
		tdyy = tdyy + kernel(Y{k1}, Y{k2}, N);
	    end
	end
	dyy = dyy + tdyy * 2 + NY;
    end

    if NX ~= 0 && NY ~= 0
	tdxy = 0;
	for k1 = 1:NX
	    for k2 = 1:NY
		tdxy = tdxy + kernel(X{k1}, Y{k2}, N);
	    end
	end
	dxy = dxy + tdxy;
    end
end

if spikeTrains1.N == 0 || spikeTrains2.N == 0
    div = 0;
else
    div = dxx/spikeTrains1.N^2 + dyy/spikeTrains2.N^2 ...
	- 2 * dxy/spikeTrains1.N/spikeTrains2.N;
end
