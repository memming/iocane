% Generate bivariate Poisson process that are dependent via synchrony.
% Two Poisson processes will have their own statistics,
% and a common Poisson process will be superpositioned to both as synchrony.
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

N = 40; % Number of realizations
M = 100; % Number of sets of realizations (= number of points for H1)

tOffset = 0.2;
jitterSigma = 0.01;
duration = 0.5;

% lambda = 4;
% lambdaSync = 16;
% cc = lambdaSync / (lambda + lambdaSync);
cc = 0.5;
totalFiringRate = 15;
lambdaSync = totalFiringRate * cc;
lambda = totalFiringRate - lambdaSync;

spikeTrains.N = N;
spikeTrains.duration = duration + tOffset * 2;
spikeTrains.source = '$Id$';
spikeTrains.data = cell(N, 1);
spikeTrains.samplingRate = Inf;

for kM = 1:M
    spikeTrains1(kM) = spikeTrains;
    spikeTrains2(kM) = spikeTrains;
    for k = 1:N
	common = tOffset + (rand(poissrnd(lambda), 1)) * duration;
	st = tOffset + (rand(poissrnd(lambdaSync), 1)) * duration;
	common2 = tOffset + (rand(poissrnd(lambda), 1)) * duration;
	if jitterSigma ~= 0
	    st1 = sort(unique([st; common + jitterSigma*rand(size(common))]));
	    st3 = sort(unique([st; common2 + jitterSigma*rand(size(common2))]));
	else
	    st1 = sort(unique([st; common]));
	    st3 = sort(unique([st; common2]));
	end

	spikeTrains1(kM).data{k} = st1;
	spikeTrains2(kM).data{k} = st3;
    end
end

depMeasures = {...
    @depSPD, divSPDParams_I('identity'); ...
    @depSPD, divSPDParams_I('int_exp'); ...
    @depSPD, divSPDParams_I('int_exp', 3); ...
    @depSPD, divSPDParams_I('int_exp', 6); ...
    @depSPD, divSPDParams_I('int_exp', 12); ...
    @depSPD, divSPDParams_I('exp_int'); ...
    @depSPD, divSPDParams_I('exp_int', 3); ...
    @depSPD, divSPDParams_I('exp_int', 6); ...
    @depSPD, divSPDParams_I('exp_int', 12); ...
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
    fprintf('.3%f - \t %s %s\n', power(k), func2str(depMeasures{k,1}), dparams2str(depMeasures{k,2}));
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
