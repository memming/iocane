function [spikeTrains] = generateRealizationsFPPM(fppm, N)
% [spikeTrains] = generateRealizationsFPPM(fppm, N)
% Generates N realizations from a FPPM structure
% Input:
%   fppm: (strcut) FPPM structure (See estimateFPPM)
%   N: (natural number) Number of realizations to generate
% Output:
%   spikeTrains: (struct) spikeTrains structure
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

spikeTrains.N = N;
spikeTrains.type = 'FinitePointProcessModel';
spikeTrains.source = '$Id$';
spikeTrains.param = fppm;
spikeTrains.duration = fppm.duration;
spikeTrains.data = cell(N, 1);

for k = 1:N
    m = find(cumsum(fppm.prN) < rand, 1, 'last');
    if isempty(m)
        spikeTrains.data{k} = [];
        continue;
    end
    templateIdx = floor(rand * (fppm.histM(m+1))) + 1;
    st = fppm.subSt{m}{templateIdx};
    sigma = fppm.kernelSizeHandle(m, fppm.histM(m+1));
    st = sort(st + randn(size(st)) * sigma);

    % TODO: what about fppm.duration? just warnings? truncate the PDF?
    %{
    if st(1) < 0
	warning('Negative spike timing');
    elseif st(end) > fppm.duration
	warning('Last spike exceeded the max time [%f > %f] [sigma = %f]', st(end), fppm.duration, sigma);
    end
    %}

    spikeTrains.data{k} = st;
end
% vim:ts=8:sts=4:sw=4
