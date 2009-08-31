function [div] = divHilbertian(spikeTrains1, spikeTrains2, params)
% Divergence between estimated the finite point processs via Hilbertian metric
% div = divHilbertian(spikeTrains1, spikeTrains2, params)
% 
% Input:
%   spikeTrains1, spikeTrains2: (struct) 2 sets of spike trains for comparison
%   params: (struct) see divHilbertianParams
% Output:
%   div: (1) divergence value
%
% See also: divHilbertianParams
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

% FPPM case
fppm1 = estimateFPPM(spikeTrains1, params.kernelSizeHandle, params.sigma1);
fppm2 = estimateFPPM(spikeTrains2, params.kernelSizeHandle, params.sigma1);
if params.isSampleOnly
    div = fppHilbertianMetricSamples(fppm1, fppm2, params.dist2Handle);
else
    div = fppHilbertianMetricMC(fppm1, fppm2, params.dist2Handle, params.NMC);
end
% vim:ts=8:sts=4:sw=4
