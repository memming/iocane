function [afppm] = estimateAFPPM(spikeTrains, kNearest)
% Estimate adaptive FPPM via nonparametric estimator using Parzen window
% afppm = estimateAFPPM(spikeTrains, kernelSizeHandle)
%
% Input:
%   spikeTrains: (struct) spikeTrains structure (see README.txt)
%   kNearest: (1) the k for k-NN based adaptive kernel size
% Output:
%   afppm: (struct) Estimated structure
%
% See also: estimateFPPM, likelihoodAFPPM
%
% In addition to FPPM, AFPPM computes the pairwise distance in each dimension
% and save the k-nearest neighbor, so that it can be used for adaptive kernel
% density estimation.
%
% Sohan Seth gave the idea to try this kind of estimator.
%
% $Id$
% Copyright 2010 iocane. All rights reserved.

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

afppm.spikeTrains = spikeTrains;
afppm.duration = spikeTrains.duration;
afppm.M = cellfun('length', spikeTrains.data);
afppm.maxM = max(afppm.M);
afppm.histM = histc(afppm.M, 0:1:afppm.maxM);

afppm.prN = afppm.histM / spikeTrains.N;
afppm.likelihood = @likelihoodAFPPM;

afppm.subSt = cell(afppm.maxM, 1);
for i = 1:afppm.maxM
    N = afppm.histM(i+1);
    subArray = zeros(N, i);
    kList = find(afppm.M == i);
    for kIdx = 1:length(kList)
	k = kList(kIdx);
	subArray(kIdx,:) = spikeTrains.data{k};
    end
    afppm.subSt{i} = subArray;

    if N == 0
	afppm.sigmas{i} = [];
	continue;
    end

    % compute the pairwise distance
    d = zeros(N, N);
    for k = 1:N-1
	for kk = (k+1):N
	    d(k, kk) = sum(subArray(k,:) - subArray(kk,:)).^2;
	end
    end
    d = d + d';
    d = sort(d, 2, 'ascend');
    dkNN = zeros(N, 1);
    if kNearest+1 > N
	dkNN = d(:,end);
    else
	dkNN = d(:,kNearest+1);
    end
    afppm.sigmas{i} = (10 / kNearest * sqrt(dkNN)).^(i+1);
end

% vim:ts=8:sts=4:sw=4
