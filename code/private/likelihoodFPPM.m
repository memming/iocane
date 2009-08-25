function [likelihood] = likelihoodFPPM(fppm, spikeTrain)
% likelihood that cares about the order of APs
% likelihood = likelihoodFPPM(fppm, spikeTrain)
%
% Input:
%   fppm: (struct) FPPM struct (see estimateFPPM)
%   spikeTrain: (mx1) a single spike train represented as sorted times
% Output:
%   likelihood: (double) likelihood value
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

m = length(spikeTrain);
if m > fppm.maxM
    likelihood = 0;
    return;
end
sigma = fppm.kernelSizeHandle(m, fppm.histM(m+1));

if m == 0
    likelihood = fppm.prN(m+1);
    return;
end

if fppm.histM(m+1) == 0
    likelihood = 0;
    return;
end

dd = 0;
for k = 1:fppm.histM(m+1)
    d = normpdf(fppm.subSt{m}{k} - spikeTrain, zeros(m, 1), sigma * ones(m, 1));
    dd = dd + prod(d);
end

likelihood = fppm.prN(m+1) * dd / fppm.histM(m+1);
% vim:ts=8:sts=4:sw=4
