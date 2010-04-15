% Hypothesis test example experiment
% - Precisely Timed Spike Train (PTST) vs Poisson approximation
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

N = 40; % Number of realizations
%M = 46; % Number of point processes per class
M = 13; % Number of point processes per class

L = 4; % Number of precisely timed action potentials (or bumps in Poisson case)
T = 1; % total duration

% The mean position of APs
mu = rand(L, 1) * T/2 + T/4;
% The std of each AP
sigma = rand(L, 1) * 0.010 + 0.001;
% Probability of lossing each AP
p = rand(L, 1);
npcum = cumsum(p); % normalized cumulative for randomly choosing one for Poisson
npcum = npcum / npcum(end);

spikeTrains.N = N;
spikeTrains.duration = T;
spikeTrains.source = '$Id$';
spikeTrains.data = cell(N, 1);
spikeTrains.samplingRate = Inf;

for kM = 1:M
    spikeTrains1(kM) = spikeTrains; % PTST
    spikeTrains2(kM) = spikeTrains; % Poisson

    for k = 1:N
	st = [];
	for kk = 1:L
	    if rand < p(kk)
		st = [st; randn * sigma(kk) + mu(kk)];
	    end
	end
	spikeTrains1(kM).data{k} = sort(st);

	st = [];
	nPoiss = poissrnd(sum(p));
	for kk = 1:nPoiss
	    kkk = find(rand < npcum, 1, 'first');
	    st = [st; randn * sigma(kkk) + mu(kkk)];
	end
	spikeTrains2(kM).data{k} = sort(st);
    end
end

divMeasures = {...
    @divL2CuIF, []; ...
    %@divSPD, []; ...
    %@divSPD, divSPDParams_I('int_exp'); ...
    %@divSPD, divSPDParams_I('exp_int'); ...
    %@divCDF, divCDFParams(Inf, 'sum'); ...
    %@divCDF, divCDFParams(2, 'sum'); ...
    %@divCDF, divCDFParams(Inf, 'sup'); ...
%    @divISF, []; ...
};

[p, power, dist, d12] = evaluateExperiment(spikeTrains1, spikeTrains2, M, 0.05, true, divMeasures);

%[p, power, dist, d12] = evaluateExperiment(spikeTrains1, spikeTrains2, M);
% vim:ts=8:sts=4:sw=4
