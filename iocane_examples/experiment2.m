% Hypothesis test example experiment
% - generate two homogeneous Poisson processes
%   One with mean rate 3/trial and the other with 5/trial
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

rand('seed', 20090523);
randn('seed', 20090523);

N = 40; % Number of realizations
M = 46; % Number of point processes per class
%M = 10; % Number of point processes per class

lambda1 = 3;
lambda2 = 5;
tOffset = 0.2;
duration = 0.05;

spikeTrains.N = N;
spikeTrains.duration = 2 * tOffset + duration;
spikeTrains.source = '$Id$';
spikeTrains.data = cell(N, 1);
spikeTrains.samplingRate = Inf;

for kM = 1:M
    spikeTrains1(kM) = spikeTrains;
    spikeTrains2(kM) = spikeTrains;

    for k = 1:N
	spikeTrains1(kM).data{k} = ...
	    tOffset + sort(rand(poissrnd(lambda1), 1)) * duration;
	spikeTrains2(kM).data{k} = ...
	    tOffset + sort(rand(poissrnd(lambda2), 1)) * duration;
    end
end

% divMeasures = { @divCDF, Inf };
divMeasures = { ...
    @divH, []; ...
    %@divL2CuIF, []; ...
    %@divSPD, [] ; ...
    %@divPhi, divPhiParams('Hellinger', 'default', 10e-3); ...
    %@divPhi, divPhiParams('Hellinger', 'silverman', 10e-3); ...
    %@divPhi, divPhiParams('Hellinger', 'modsilverman', 10e-3);
    %@divSPD, divSPDParams_I('int_exp'); ...
    %@divSPD, divSPDParams_I('exp_int'); ...
    %@divCount, divCountParams('CM'); ...
    %@divCDF, Inf ; ...
    %@divICDF, [] ; ...
    %@divCount, divCountParams('KS') ; ...
};
% divMeasures = {...
%     @divCount, divCountParams('KS'); ...
% };
[p, power, dist, d12] = evaluateExperiment(spikeTrains1, spikeTrains2, M, 0.05, true, divMeasures);

%[p, power, dist, d12] = evaluateExperiment(spikeTrains1, spikeTrains2, M);
% vim:ts=8:sts=4:sw=4
