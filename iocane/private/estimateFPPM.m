function [fppm] = estimateFPPM(spikeTrains, kernelSizeHandle, sigma1)
% Estimate FPPM via nonparametric estimator using Parzen window
% fppm = estimateFPPM(spikeTrains, kernelSizeHandle)
%
% Input:
%   spikeTrains: (struct) spikeTrains structure (see README.txt)
%   kernelSizeHandle: (@(dim,N)) kernel size handle
%   sigma1: (1) kernel size in 1-D
% Output:
%   fppm: (struct) Estimated structure
%
% See also: likelihoodFPPM, generateRealizationsFPPM, divHilbertian
%
% References
% [1] Il Park, Sohan Seth, Jose C. Principe. "Divergence on finite point 
%   processes for multiple trial spike train observations",
%   (submitted to NIPS 2009)
%
% $Id$
% Copyright 2009 Memming. All rights reserved.

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

fppm.kernelSizeHandle = kernelSizeHandle;
fppm.sigma1 = sigma1;
fppm.spikeTrains = spikeTrains;
fppm.duration = spikeTrains.duration;
fppm.M = cellfun('length', spikeTrains.data);
fppm.maxM = max(fppm.M);
fppm.histM = histc(fppm.M, 0:1:fppm.maxM);

fppm.prN = fppm.histM / spikeTrains.N;

fppm.subSt = cell(fppm.maxM, 1);
for i = 1:fppm.maxM
    subArray = zeros(fppm.histM(i+1), i);
    kList = find(fppm.M == i);
    for kIdx = 1:length(kList)
	k = kList(kIdx);
	subArray(kIdx,:) = spikeTrains.data{k};
    end
    fppm.subSt{i} = subArray;
end

% vim:ts=8:sts=4:sw=4
