function [p, power] = evaluateExperiment(spikeTrains1, spikeTrains2, M, alpha, verbose, divMeasures)
% Compute the statistical power for each divergence method
% [p, power] = evaluateExperiment(spikeTrains1, spikeTrains2, M, alpha, verbose, divMeasures)
%
% Perform hypothesis testing given lots of samples from the null hypothesis
% (spikeTrains1) and alternative hypothesis (spikeTrains2) and return the
% statistical power of the test.
%
% Input:
%   spikeTrains1: independent spike trains from the null hypothesis
%   spikeTrains2: independent spike trains from the alternate hypothesis
%   M: number of spike trains (min spikeTrains1, spikeTrains2)
%   alpha: (default: 0.05) threshold for type-I error
%   verbose: (default: true) print result
%   divMeasures: (default) cell array of div method and parameters
% Output:
%   p: p-value of samples from H1 on empirical distribution of H0
%   power: estimated statistical power for each divergence method
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

if nargin < 4
    alpha = 0.05;
end

if nargin < 5
    verbose = true;
end

% List of divergence measures that will be tested
if nargin < 6
divMeasures = {...
    @divMeanFiringRate, []; ...
    @divFF, []; ...
    @divISI, divISIParams('KS'); ...
    @divCount, divCountParams('KS'); ...
    @divL2Poisson, divL2PoissonParams(1e-3, 'fixed', 100e-3); ...
    @divL2Poisson, divL2PoissonParams(1e-3, 'fixed', 10e-3); ...
    %@divL2Poisson, divL2PoissonParams(1e-3, 'optimal'); ...
    @divTTFS, divTTFSParams('KS'); ...
    @divHilbertian, divHilbertianParams('Hellinger', 'default', 100e-3); ...
    @divHilbertian, divHilbertianParams('Hellinger', 'default', 10e-3); ...
    };
end

for k = 1:size(divMeasures, 1)
    divHandle = divMeasures{k,1};
    divParams = divMeasures{k,2};
    % create null hypothesis
    dist{k} = empDivDist(spikeTrains1, divHandle, divParams);
    kk = 1;
    for k1 = 1:(M-1)
	for k2 = (k1+1):M
	    d12 = divHandle(spikeTrains1(k1), spikeTrains2(k2), divParams);
	    p(k,kk) = dist{k}.pValue(d12);
	    kk = kk + 1;
	end
    end
end

if verbose; fprintf('Statistical power with alpha = %f\n', alpha); end
for k = 1:size(divMeasures, 1)
    power(k) = sum(p(k,:) < alpha) / size(p, 2);
    if verbose
	fprintf('%f - \t %s\n', power(k), func2str(divMeasures{k,1}));
    end
end
% vim:ts=8:sts=4:sw=4
