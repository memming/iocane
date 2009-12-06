% Consistency for hypothesis test example experiment
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

rand('seed', 20091105);
randn('seed', 20091105);

greyColor = [0.7 0.7 0.7];
blueColor = [0.7 0.7 1];

% m is realization for Monte Carlo integration
mRange = ceil(logspace(1,log10(500),10)); % 1 to 1000
mMax = length(mRange);
% n is realization of multiple Monte Carlo's
nRealizations = 100;
alpha = 0.05;

% divHandle = @divHilbertian;
% divParams = divHilbertianParams('Hellinger', 'default', 10e-3);

divHandle = @divPhi;
divParams = divPhiParams('Hellinger', 'default', 10e-3);
% divParams = divPhiParams('total variation', 'default', 10e-3);

param.T = T;
param.mISI = mISI;
param.urISI = urISI;
param.type = 'correlated';
spikeTrains1 = genSerialCorr(N, M, param);
param.type = 'uncorrelated'
spikeTrains2 = genSerialCorr(N, M, param);

for mIdx = 1:length(mRange)
    m = mRange(mIdx);
    if mIdx > 1
	howManyMore = m - mRange(mIdx-1);
    else
	howManyMore = m;
    end

    for n = 1:nRealizations
	spikeTrains1(n).N = m;
	spikeTrains2(n).N = m;

	for k = 1:howManyMore
	    st = [];
	    for kk = 1:L
		if rand < p(kk)
		    st = [st; randn * sigma(kk) + mu(kk)];
		end
	    end
	    spikeTrains1(n).data{m-k+1} = sort(st);

	    st = [];
	    nPoiss = poissrnd(sum(p));
	    for kk = 1:nPoiss
		kkk = find(rand < npcum, 1, 'first');
		st = [st; randn * sigma(kkk) + mu(kkk)];
	    end
	    spikeTrains2(n).data{m-k+1} = sort(st);
	end
    end

    for k1 = 1:nRealizations
	d2(mIdx,k1) = divHandle(spikeTrains1(k1), spikeTrains2(k1), divParams);
    end
    for k1 = 1:nRealizations-1
	d2a(mIdx,k1) = divHandle(spikeTrains1(k1), spikeTrains1(k1+1), divParams);
	d2b(mIdx,k1) = divHandle(spikeTrains2(k1), spikeTrains2(k1+1), divParams);
    end
    d2a(mIdx,nRealizations) = divHandle(spikeTrains1(k1+1), spikeTrains1(1), divParams);
    d2b(mIdx,nRealizations) = divHandle(spikeTrains2(k1+1), spikeTrains2(1), divParams);

    d2a_sorted = sort(d2a(mIdx,:),'ascend');
    d2b_sorted = sort(d2b(mIdx,:),'ascend');
    thresholdIdx = (ceil(nRealizations * (1-alpha)));
    pa = sum(d2(mIdx,:) > d2a_sorted(thresholdIdx))/nRealizations;
    pb = sum(d2(mIdx,:) > d2b_sorted(thresholdIdx))/nRealizations;
    fprintf('m: %d done...%f\t%f\n', m, pa, pb);
end

figure;
hold on;
for k = 1:20 %nRealizations
    plot(mRange, d2(:,k), 'o-', 'Color', greyColor + [0 0 0.2]);
    plot(mRange, d2a(:,k), 'd-', 'Color', greyColor + [0.1 0 0]);
    plot(mRange, d2b(:,k), 'x-', 'Color', greyColor + [0 0.2 0]);
end
ph1 = errorbar(mRange, mean(d2,2), std(d2,[],2), 'b--o');
ph2 = errorbar(mRange, mean(d2a,2), std(d2a,[],2), 'r--d');
ph3 = errorbar(mRange, mean(d2b,2), std(d2b,[],2), 'g--x');
ylim([0,2]);
xlabel('# of samples');
ylabel('HL10');
set(gca, 'XTick', mRange)
set(gca, 'XScale', 'log')
xlim([mRange(1) mRange(end)])
legend([ph1 ph3 ph2], 'H_0 vs H_1', 'H_0 (Poisson)', 'H_0 (PTST)');
% vim:ts=8:sts=4:sw=4
