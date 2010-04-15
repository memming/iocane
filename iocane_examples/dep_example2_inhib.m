% Neuron A is inhibiting neuron B
% Neuron A, B are homogeneously Poisson except that
% whenever neuron A fires, spikes of neuron B within 10 ms
%                          are silenced with a given probability.
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

rand('seed', 20100316);
randn('seed', 20100316);

%N = 40; % Number of realizations
N = 100; % Number of realizations
M = 100; % Number of sets of realizations (= number of points for H1)

tOffset = 0.1;
jitterSigma = 0.01;
duration = 1;

lambdaA = 20;
lambdaB = 5;
inhibitionTime = 10e-3;
probInhibition = 0.5;
fprintf('Lambda A %f B %f inhibitionTime %f\n', lambdaA, lambdaB, inhibitionTime);

spikeTrains.N = N;
spikeTrains.duration = duration + tOffset * 2;
spikeTrains.source = '$Id$';
spikeTrains.data = cell(N, 1);
spikeTrains.samplingRate = Inf;

for kM = 1:M
    spikeTrains1(kM) = spikeTrains;
    spikeTrains2(kM) = spikeTrains;
    for k = 1:N
	stA = tOffset + (rand(poissrnd(lambdaA), 1)) * duration;
	stB = tOffset + (rand(poissrnd(lambdaB), 1)) * duration;
	di = (repmat(stB, 1, size(stA,1)) - repmat(stA, 1, size(stB, 1))');
	ij = (di >= 0 & di < inhibitionTime);
	idx = (sum(ij,2) ~= 0) & (rand(size(stB,1),1) >= probInhibition);
	stB(idx) = [];
	%stB = stB(sum(ij,2) == 0);

	spikeTrains1(kM).data{k} = stA;
	spikeTrains2(kM).data{k} = stB;
    end
end

%{
binSize = 1e-3;
maxLag = 100;
for kM = 1:1
    stsA = spikeTrains1(kM);
    stsB = spikeTrains2(kM);
    xc = zeros(maxLag*2+1, 1);
    for k = 1:N
	stA = binSpikeTrain(stsA.data{k}, stsA.duration, binSize);
	stB = binSpikeTrain(stsB.data{k}, stsB.duration, binSize);
	xc = xc + xcorr(stA, stB, maxLag);
    end
    figure(77); hold all; plot(xc);
end
return;
%}

depMeasures = {...
    %@depSPD, divSPDParams_I('int_exp'); ...
    % @depSPD, divSPDParams_I('int_exp', 4); ...
    % @depSPD, divSPDParams_I('int_exp', 8); ...
    %@depSPD, divSPDParams_I('exp_int'); ...
    %@depSPD, divSPDParams_I('exp_int', 4); ...
    %@depSPD, divSPDParams_I('exp_int', 8); ...
    @depSPD, divSPDParams_nci2(10e-3, 10e-3, 'gaussian'); ...
    @depSPD, divSPDParams_nci2(10e-3, 100e-3, 'gaussian'); ...
    @depSPD, divSPDParams_I('exp_int', 16); ...
    @depSPD, divSPDParams_I('identity'); ...
};

nSurrogate = M;
depHandle = depMeasures{1,1};
depParams = depMeasures{1,2};
alpha = 0.05;

for k = 1:size(depMeasures, 1)
    for kk = 1:nSurrogate
	nn = randperm(M);
	n1 = nn(1); n2 = nn(2);
	dSurr(k,kk) = depHandle(spikeTrains1(n1), spikeTrains2(n2), depParams);
    end

    for kk = 1:M
	d(k,kk) = depHandle(spikeTrains1(kk), spikeTrains2(kk), depParams);
    end

    power(k) = (sum(d(k,:) >= quantile(dSurr(k,:), 1 - alpha)) / M);
    fprintf('%f - \t %s %s\n', power(k), func2str(depMeasures{k,1}), dparams2str(depMeasures{k,2}));
end

% Expected output
% 0.840 - 	 depSPD (kernel=divSPDParams_I/k,kappa=identity,sigma=2)
% 0.810 - 	 depSPD (kernel=divSPDParams_I/k,kappa=int_exp,sigma=2)
% 0.840 - 	 depSPD (kernel=divSPDParams_I/k,kappa=int_exp,sigma=3)
% 0.840 - 	 depSPD (kernel=divSPDParams_I/k,kappa=int_exp,sigma=6)
% 0.860 - 	 depSPD (kernel=divSPDParams_I/k,kappa=int_exp,sigma=12)
% 0.800 - 	 depSPD (kernel=divSPDParams_I/k,kappa=exp_int,sigma=2)
% 0.840 - 	 depSPD (kernel=divSPDParams_I/k,kappa=exp_int,sigma=3)
% 0.770 - 	 depSPD (kernel=divSPDParams_I/k,kappa=exp_int,sigma=6)
% 0.760 - 	 depSPD (kernel=divSPDParams_I/k,kappa=exp_int,sigma=12)

% vim:ts=8:sts=4:sw=4
