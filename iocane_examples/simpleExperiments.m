function [spikeTrains1, spikeTrains2, p, power, dist, d12, rtime, divMeasures] = simpleExperiments(expNum, isPlotOn, divMeasures, Noveride, M)
% Run a set of simple experiments with fixed parameters
% [spikeTrains1, spikeTrains2, p, power, dist, d12, rtime, divMeasures] = simpleExperiments(expNum, isPlotOn, divMeasures, Noveride);
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

if nargin < 2
    isPlotOn = true;
end

if nargin < 3
divMeasures = {...
    @divL2CuIF, []; ...
    @divMeanFiringRate, []; ...
    @divFF, []; ...
    @divISI, divISIParams('KS'); ...
    %@divCount, divCountParams('KS'); ...
    %%@divH, []; ...
    %%@divISF, []; ...
    %%@divSPD_stratified, divSPDParams_stratified('gaussian', '', 0.1e-3); ...
    %%@divRatioChiSquare, divSPDParams_I('identity'); ...
    %%@divRatioChiSquare, divSPDParams_I('int_exp', 0.067); ...
    %%@divRatioChiSquare, divSPDParams_I('exp_int', 0.067); ...
    %%@divRatioChiSquare, divSPDParams_I('int_exp', 0.67); ...
    %%@divRatioChiSquare, divRatioChiSquareParams_I('exp_int', 0.067); ...
    %@divRatioChiSquare, divRatioChiSquareParams_I('exp_int', 1); ...
    %@divRatioChiSquare, divRatioChiSquareParams_I('exp_int', 0.1); ...
    %@divRatioChiSquare, divRatioChiSquareParams_I('exp_int', 'median'); ...
    %@divRatioChiSquare, divRatioChiSquareParams_I('identity'); ...
    %%@divSPD, divSPDParams_I('int_exp'); ...
    %%@divSPD, divSPDParams_I('exp_int'); ...
    %%@divPD, []; ...
    %%@divCDF, divCDFParams(Inf, 'sum'); ...
    %%@divCDF, divCDFParams(2, 'sum'); ...
    %@divL2Poisson, divL2PoissonParams(1e-3, 'fixed', 10e-3); ...
    %@divL2Poisson, divL2PoissonParams(10e-3, 'fixed', 100e-3); ...
    %@divL2Poisson, divL2PoissonParams(1e-3, 'hist'); ...
    %%@divPhi, divPhiParams('Hellinger', 'default', 10e-3);
    %%@divHilbertian, divHilbertianParams('Hellinger', 'default', 10e-3); ...
    %%@divHilbertian, divHilbertianParams('Hellinger', 'kNN', 10); ...
    %@divHilbertian, divHilbertianParams('chi-square'); ...
    %@divHilbertian, divHilbertianParams('chi-square', 'modsilverman', 1e-3); ...
    %@divHilbertian, divHilbertianParams('chi-square', 'modsilverman', 10e-3); ...
    %@divHilbertian, divHilbertianParams('chi-square', 'modsilverman', 100e-3); ...
    %@divHilbertian, divHilbertianParams('chi-square', 'silverman', 1e-3); ...
    %@divHilbertian, divHilbertianParams('chi-square', 'silverman', 10e-3); ...
    %@divHilbertian, divHilbertianParams('chi-square', 'silverman', 100e-3); ...
    %@divRatioSymmetricChiSquare, divRatioChiSquareParams_I('exp_int', 1); ...
    @divRatioSymmetricChiSquare, divRatioChiSquareParams_I('exp_int', 0.1); ...
    %@divRatioSymmetricChiSquare, divRatioChiSquareParams_I('exp_int', 'median'); ...
};
end

if nargin < 4
    Noveride = 0;
end

if nargin < 5
    M = 46;
end

switch(expNum)
case 1
fprintf('=== Experiment 1: Two correlated vs independent APs\n');
N = 40;
if Noveride ~= 0; N = Noveride; end
spikeTrains1 = genTwoAPex(N, M, struct('jitter',10e-3,'type','correlated'));
spikeTrains2 = genTwoAPex(N, M, struct('jitter',10e-3,'type','uncorrelated'));
case 2
fprintf('=== Experiment 2: Homogeneous Poisson processes\n');
N = 40;
if Noveride ~= 0; N = Noveride; end
homogeneousParam1.tOffset = 0.25;
homogeneousParam1.duration = 0.50;
homogeneousParam1.lambda = 3;
homogeneousParam2 = homogeneousParam1;
homogeneousParam2.lambda = 5;
spikeTrains1 = genHomogeneousPoisson(N, M, homogeneousParam1);
spikeTrains2 = genHomogeneousPoisson(N, M, homogeneousParam2);
case 3
fprintf('=== Experiment 3: PTST vs equi-Poisson process\n');
N = 40;
if Noveride ~= 0; N = Noveride; end
PTSTParams.L = 3;
PTSTParams.T = 1;
PTSTParams.type = 'PTST';
PoissonPTSTParams = PTSTParams;
PoissonPTSTParams.type = 'equPoisson';
spikeTrains1 = genPTST(N, M, PTSTParams);
spikeTrains2 = genPTST(N, M, PoissonPTSTParams);
case 4
fprintf('=== Experiment 4: Renewal vs Serially correlated\n');
N = 300;
if Noveride ~= 0; N = Noveride; end
renewalParam = struct('T',125e-3,'mISI',50e-3,'urISI',5e-3,'type','renewal');
serialCorrParam = struct('T',125e-3,'mISI',50e-3,'urISI',5e-3,'type','correlated');
spikeTrains1 = genSerialCorr(N, M, renewalParam);
spikeTrains2 = genSerialCorr(N, M, serialCorrParam);
case 5
fprintf('=== Experiment 5: Stratified specification (Gaussian)\n');
N = 40;
if Noveride ~= 0; N = Noveride; end
stratParam1.D = 2;
stratParam1.p = [0.5 0.5];
stratParam1.T = 1;
stratParam1.mu{1} = 0.6;
stratParam1.sigma{1} = 20e-3^2;
stratParam1.mu{2} = [0.4 0.6];
stratParam1.sigma{2} = [20e-3^2, 0; 0, 20e-3^2];
stratParam2 = stratParam1;
stratParam2.mu{2} = [0.4 0.5];
spikeTrains1 = genStratifiedGaussian(N, M, stratParam1);
spikeTrains2 = genStratifiedGaussian(N, M, stratParam2);
case 6
fprintf('=== Experiment 6: Step Poisson process with different rate\n');
N = 40;
if Noveride ~= 0; N = Noveride; end
stepPoissonParam1 = struct('sectionLength', 100e-3, 'lambda1', 3, 'lambda2', 2);
stepPoissonParam2 = struct('sectionLength', 100e-3, 'lambda1', 2, 'lambda2', 3);
spikeTrains1 = genStepPoisson(N, M, stepPoissonParam1);
spikeTrains2 = genStepPoisson(N, M, stepPoissonParam2);
case 7
fprintf('=== Experiment 7: Renewal process with different interval structure\n');
N = 40;
if Noveride ~= 0; N = Noveride; end
renewalParam1.T = 1;
renewalParam1.rate = 10;
renewalParam2 = renewalParam1;
renewalParam1.shape = 3;
renewalParam2.shape = 0.5;
spikeTrains1 = genRenewalGamma(N, M, renewalParam1);
spikeTrains2 = genRenewalGamma(N, M, renewalParam2);
case 8
fprintf('=== Experiment 8: Renewal process vs Poisson\n');
N = 25;
if Noveride ~= 0; N = Noveride; end
renewalParam1.T = 1;
renewalParam1.rate = 10;
%renewalParam1.shape = 2;
renewalParam1.shape = 10;
spikeTrains1 = genRenewalGamma(N, M, renewalParam1);
homogeneousParam1.tOffset = 0;
homogeneousParam1.duration = 1;
homogeneousParam1.lambda = 10;
spikeTrains2 = genHomogeneousPoisson(N, M, homogeneousParam1);

otherwise
    error('No such experiment number');
end
fprintf('Data generation complete...[N = %d, M = %d]\n', N, M);

if isPlotOn
    for m = 1:M
	figure(777); clf;
	subplot(2,1,1); plotRaster(spikeTrains1(m));
	xlim([0 spikeTrains1(1).duration])
	subplot(2,1,2); plotRaster(spikeTrains2(m));
	xlim([0 spikeTrains2(1).duration])
	drawnow;
	reply = input('Start computing (enter) next sample (any key): ', 's');
	if isempty(reply)
	    break;
	end
    end
end

alpha = 0.05;
[p, power, dist, d12, rtime] = evaluateExperiment(spikeTrains1, spikeTrains2, M, alpha, true, divMeasures);

if isPlotOn
    cm = colormap('lines');
    figure(778); clf; hold all;
    for k = 1:size(divMeasures,1)
	% ksdensity(dist{k}.values); ksdensity(d12(k,:));
	[h0, x0] = hist(dist{k}.values, 20);
	h0 = h0 / sum(h0 * (x0(2) - x0(1)));
       	[h12, x12] = hist(d12(k,:), 20);
	h12 = h12 / sum(h12 * (x12(2) - x12(1)));
	ph(k) = plot(x0, h0, '.-', 'Color', cm(k,:));
	plot(x12, h12, 'o:', 'Color', cm(k,:));
	line(dist{k}.threshold(0.05) * [1 1], [0 max(h0)], 'Color', cm(k,:));
	lgdStr{k} = func2str(divMeasures{k,1});
    end
    legend(ph, lgdStr, 'Location', 'North');
    set(gca, 'YTick', []);
end

% vim:ts=8:sts=4:sw=4
