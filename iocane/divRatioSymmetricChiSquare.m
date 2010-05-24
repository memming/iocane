function [div] = divRatioSymmetricChiSquare(spikeTrains1, spikeTrains2, params)
% Kernel estimation of R-N derviative in symmetric chi-square divergence
% div = divRatioSymmetricChiSquare(spikeTrains1, spikeTrains2, params)
% The Radon-Nikodym derivative is estimated using regularized kernel least
% squares using the given kernel.
%
% Input:
%   spikeTrains1, spikeTrains2: (struct) sets of spike trains for comparison
%   params: (struct) use divRatioChiSquareParams_I
% Output:
%   div: (1) divergence value
%
% g = argmin_f \int (dP/d(P+Q) - f)^2 dQ
% div = \int (g - 1)^2 d(P+Q)/2
%
% See also: divRatioChiSquare, divRatioChiSquareParams_I
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
lambdaZZ = 1 / (spikeTrains1.N + spikeTrains2.N);

kernel = params.kernel;

if spikeTrains1.N == 0 || spikeTrains2.N == 0
    div = 0;
    return
end

[Kxx, Kxy, Kyy] = kernel(spikeTrains1, 1:spikeTrains1.N, spikeTrains2, 1:spikeTrains2.N);

% When P+Q is used as the basis
%
% Kzy = [Kxy; Kyy];
% Kzx = [Kxx; Kxy'];
% Kzz = [Kxx Kxy; Kxy' Kyy]; % = [Kzx Kzy];
% 
% alphaZY = (Kzz * Kzz + lambdaZZ * eye(size(Kzz))) ...
%     \ (Kzy * ones(spikeTrains2.N, 1));
% divZY = mean((Kzz * alphaZY - 1).^2);
% 
% alphaZX = (Kzz * Kzz + lambdaZZ * eye(size(Kzz))) ...
%     \ (Kzx * ones(spikeTrains1.N, 1));
% divZX = mean((Kzz * alphaZX - 1).^2);
% 
% div = (divZX + divZY) / 2;

% When P is used for the basis and dP/d(P+Q) is evaluated
% ((Kxy * Kyx + Kxx * Kxx) / 2)^-1 * Kxx * 1
% If median kernel size is used, it should be estimated for P
alphaXY = (Kxy * Kxy' + Kxx * Kxx + lambdaXX * eye(size(Kxx)) / 2) ...
    \ Kxx * ones(spikeTrains1.N, 1);
divXY = mean((alphaXY' * Kxy - 1).^2);

% When Q is used for the basis and dQ/d(P+Q) is evaluated
% ((Kyx * Kxy + Kyy * Kyy) / 2)^-1 * Kyy * 1
% If median kernel size is used, it should be estimated for Q
if strcmp(params.sigma, 'median')
    [Kxx, Kxy, Kyy] = kernel(spikeTrains2, 1:spikeTrains2.N, spikeTrains1, 1:spikeTrains1.N); % Now X is Y and Y is X. That's confusing! :P
end
alphaYX = (Kxy * Kxy' + Kxx * Kxx + lambdaYY * eye(size(Kxx)) / 2) ...
    \ Kxx * ones(spikeTrains1.N, 1);
divYX = mean((alphaYX' * Kxy - 1).^2);
% alphaYX = ((Kxy' * Kxy + Kyy * Kyy + lambdaYX * eye(size(Kyy)) / 2) ...
%     \ Kyy * ones(spikeTrains2.N, 1);
% divYX = mean((Kxy' * alphaYX - 1).^2);

div = (divXY + divYX) / 2;
