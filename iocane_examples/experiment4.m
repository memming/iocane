% Hypothesis test example experiment
% - Serialy correlated vs uncorrelated renewal processes
%
% Reference: Maurice J. Chacron, Benjamin Lindner, André Longtin.
%  Noise Shaping by Interval Correlations Increases Information Transfer.
%  Physical Review Letters, Vol. 92, No. 8. (25 Feb 2004), 080601.
%  10.1103/PhysRevLett.92.080601
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

rand('seed', 20091026);
randn('seed', 20091026);

N = 320; % Number of realizations
M = 46; % Number of point processes per class

T = 125e-3; % total duration
mISI = 50e-3; % mean ISI
urISI = 5e-3; % 5 ms uniform distributed jitter to half ISI
mISI = mISI - urISI;

spikeTrains.N = N;
spikeTrains.duration = T;
spikeTrains.source = '$Id$';
spikeTrains.data = cell(N, 1);
spikeTrains.samplingRate = Inf;

for kM = 1:M
    spikeTrains1(kM) = spikeTrains; % PTST
    spikeTrains2(kM) = spikeTrains; % Poisson

    for k = 1:N
	st = [0];
	while st(end) < T
	    isi = mISI + rand * urISI + rand * urISI;
	    st = [st; st(end) + isi];
	end
	st = st(2:end-1);
	spikeTrains1(kM).data{k} = st;

	st = [0];
	lastJitter = rand * urISI;
	while st(end) < T
	    currentJitter = rand * urISI;
	    isi = mISI + currentJitter + lastJitter;
	    lastJitter = currentJitter;
	    st = [st; st(end) + isi];
	end
	st = st(2:end-1);
	spikeTrains2(kM).data{k} = st;
    end
end

evaluateExperiment(spikeTrains1, spikeTrains2, M);
% vim:ts=8:sts=4:sw=4
