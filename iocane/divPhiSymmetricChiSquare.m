function [div] = divPhiSymmetricChiSquare(spikeTrains1, spikeTrains2, params)
% Divergence between estimated the finite point processs via symmetric chi^2
% [div] = divPhiSymmetricChiSquare(spikeTrains1, spikeTrains2, params)
% NOT same as divPhi(P, (P+Q)/2), because divPhi symmetrizes
% 
% Input:
%   spikeTrains1, spikeTrains2: (struct) 2 sets of spike trains for comparison
%   params: (struct) see divPhiSymmetricChiSquareParams
% Output:
%   div: (1) divergence value
%
% See also: divPhi, divRatioChiSquare
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

% spikeTrainsPQ = mergeSpikeTrains(spikeTrains1, spikeTrains2);

fppm1 = estimateFPPM(spikeTrains1, params.kernelSizeHandle, params.sigma1);
fppm2 = estimateFPPM(spikeTrains2, params.kernelSizeHandle, params.sigma1);

dTemp1 = zeros(fppm1.spikeTrains.N, 1);
for k = 1:fppm1.spikeTrains.N
    st = fppm1.spikeTrains.data{k};
    j1 = likelihoodFPPM(fppm1, st);
    j2 = likelihoodFPPM(fppm2, st);
    dTemp1(k) = (j2/(j1+j2) - 1)^2;
end

dTemp2 = zeros(fppm2.spikeTrains.N, 1);
for k = 1:fppm2.spikeTrains.N
    st = fppm2.spikeTrains.data{k};
    j1 = likelihoodFPPM(fppm1, st);
    j2 = likelihoodFPPM(fppm2, st);
    dTemp2(k) = (j2/(j1+j2) - 1)^2;
end

div = mean([dTemp1; dTemp2]);

if isnan(div)
    error('NaN in $Id$');
end

if isinf(div)
    error('Inf in $Id$');
end

function spikeTrains3 = mergeSpikeTrains(spikeTrains1, spikeTrains2)
    spikeTrains3 = spikeTrains1;
    spikeTrains3.N = spikeTrains1.N + spikeTrains2.N;
    spikeTrains3.data = cell(spikeTrains3.N, 1);
    for k = 1:spikeTrains1.N
	spikeTrains3.data{k} = spikeTrains1.data{k};
    end
    for k = 1:spikeTrains2.N
	spikeTrains3.data{k+spikeTrains1.N} = spikeTrains2.data{k};
    end
end

end
% vim:ts=8:sts=4:sw=4
