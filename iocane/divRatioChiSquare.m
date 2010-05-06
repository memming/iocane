function [div] = divRatioChiSquare(spikeTrains1, spikeTrains2, params)
% Kernel estimation of the Radon-Nikodym derviative in Chi-square divergence
% [div] = divRatioChiSquare(spikeTrains1, spikeTrains2, params)
% The Radon-Nikodym derivative is estimated using regularized kernel least
% squares using the given kernel.
%
% Input:
%   spikeTrains1, spikeTrains2: (struct) sets of spike trains for comparison
%   params: (struct) see divSPDParams
% Output:
%   div: (1) divergence value
%
% g = argmin_f \int (dP/dQ - f)^2 dQ
% div = \int (g - 1)^2 dQ
%
% See also: 
%
% Original idea by Sohan Seth
% Caution: We haven't showed which SPD kernels are universal. Only universal
%    kernels are guarantteed to work in this approach.
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

lambdaXX = 1 / (spikeTrains1.N);
lambdaYY = 1 / (spikeTrains2.N);

kernel = params.kernel;

if spikeTrains1.N == 0 || spikeTrains2.N == 0
    div = 0;
    return
end

[Kxx, Kxy, Kyy] = kernel(spikeTrains1, 1:spikeTrains1.N, spikeTrains2, 1:spikeTrains2.N);

alphaXX = (Kxx * Kxx + lambdaXX * eye(size(Kxx))) ...
    \ (Kxy * ones(spikeTrains1.N, 1));
divXX = mean((Kxx * alphaXX - 1).^2);

alphaYY = (Kyy * Kyy + lambdaYY * eye(size(Kyy))) ...
    \ (Kxy' * ones(spikeTrains2.N, 1));
divYY = mean((Kyy * alphaYY - 1).^2);

div = (divXX + divYY) / 2;
