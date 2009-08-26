function d = fppHilbertianMetricSamples(fppm1, fppm2, dist2Handle)
% Compute the Hilbertian metric between finite point process models
% d = fppHilbertianMetricSamples(fppm1, fppm2, dist2Handle)
% using Monte Carlo method, no samples are generated. All samples
% in the FPP model are used.
%
% Input:
%   fppm1, fppm2: (struct) FPPM (see estimateFPPM)
%   dist2Handle: (@) 1-homogeneous suqare metric on positive reals
%                 (see dist2HandleFactory)
% Output:
%   d: (double) computed square divergence
%
% References
% [1] Il Park, Sohan Seth, Jose C. Principe. "Divergence on finite point 
%   processes for multiple trial spike train observations",
%   (submitted to NIPS 2009)
%
% See also fppHilbertianMetricMC, estimateFPPM, dist2HandleFactory
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

dTemp1 = zeros(fppm1.spikeTrains.N, 1);
dTemp2 = zeros(fppm2.spikeTrains.N, 1);

for k = 1:fppm1.spikeTrains.N
    % use measure p/2
    st = fppm1.spikeTrains.data{k};
    j1 = likelihoodFPPM(fppm1, st);
    j2 = likelihoodFPPM(fppm2, st);
    dTemp1(k) = dist2Handle(j1, j2) * 2 / (j1 + j2);
end

for k = 1:fppm2.spikeTrains.N
    % use measure q/2
    st = fppm2.spikeTrains.data{k};
    j1 = likelihoodFPPM(fppm1, st);
    j2 = likelihoodFPPM(fppm2, st);
    dTemp2(k) = dist2Handle(j1, j2) * 2 / (j1 + j2);
end

d = mean([dTemp1;dTemp2]);
% vim:ts=8:sts=4:sw=4
