function simulate(simCond, simIdx)
% TODO: description
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

tic;

rand('seed', simCond.randseed);
randn('seed', simCond.randseed);

N = simCond.nRange(end);
spikeTrains1all = simCond.dataHandle1(N, simCond.M, simCond.dataHandleParams1);
spikeTrains2all = simCond.dataHandle2(N, simCond.M, simCond.dataHandleParams2);

spikeTrains1 = spikeTrains1all;
spikeTrains2 = spikeTrains2all;

divHandle = simCond.divHandle;
divParams = simCond.divParams;

for nIdx = 1:length(simCond.nRange)
    n = simCond.nRange(nIdx);

    for m = 1:simCond.M
	spikeTrains1(m).N = n;
	spikeTrains2(m).N = n;
	spikeTrains1(m).data = spikeTrains1all(m).data(1:n);
	spikeTrains2(m).data = spikeTrains2all(m).data(1:n);
    end

    empDist1{nIdx} = empDivDist(spikeTrains1, divHandle, divParams);
    empDist2{nIdx} = empDivDist(spikeTrains2, divHandle, divParams);

    for m = 1:simCond.M
	d2(nIdx,m) = divHandle(spikeTrains1(m), spikeTrains2(m), divParams);
	p1(nIdx,m) = empDist1{nIdx}.pValue(d2(nIdx,m));
	p2(nIdx,m) = empDist2{nIdx}.pValue(d2(nIdx,m));
    end
    power1(nIdx) = sum(p1(nIdx,:) < simCond.alpha) / size(p1, 2);
    power2(nIdx) = sum(p2(nIdx,:) < simCond.alpha) / size(p2, 2);
end

tictoc = toc;

save(sprintf('simResults/sim%03d', simIdx), 'power1', 'power2', 'empDist1', 'empDist2', 'd2', 'tictoc');
fprintf('Simulation for [%d] done in [%.2f min]\r', simIdx, tictoc/60);
% vim:ts=8:sts=4:sw=4
