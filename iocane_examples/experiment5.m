% 1D, 2D controlled spike trains
%
% Very very artificial example....suggested by Sohan Seth
%
% $Id$

rand('seed', 20090523);
randn('seed', 20090523);

N = 40; % Number of realizations
M = 46; % Number of point processes per class

T = 1; % total duration

mu2_1 = [0.4 0.6]';
mu2_2 = [0.2 0.6]';
%mu2_2 = mu2_1;
sigma2_1 = [0.04 0.04]';
sigma2_2 = [0.04 0.04]';

%mu1_1 = 0.4;
mu1_1 = 0.6;
mu1_2 = 0.6;
sigma1_1 = 0.04;
sigma1_2 = 0.04;

spikeTrains.N = N;
spikeTrains.duration = T;
spikeTrains.source = '$Id$';
spikeTrains.data = cell(N, 1);
spikeTrains.samplingRate = Inf;

for kM = 1:M
    spikeTrains1(kM) = spikeTrains;
    spikeTrains2(kM) = spikeTrains;

    % Identical 2D
    for k = 1:N/2
	st = [randn(2,1) .* sigma2_2 + mu2_1];
	spikeTrains1(kM).data{k} = sort(st);
	st = [randn(2,1) .* sigma2_2 + mu2_2];
	spikeTrains2(kM).data{k} = sort(st);
    end

    % Distinct 1D
    for k = N/2+1:N
	st = [randn(1,1) * sigma1_1 + mu1_1];
	spikeTrains1(kM).data{k} = st;
	st = [randn(1,1) * sigma1_2 + mu1_2];
	spikeTrains2(kM).data{k} = st;
    end
end

divMeasures = { ...
    %@divCDF, 1; ...
    %@divSPD, divSPDParams_nci2(50e-3, 1, 'gaussian'); ...
    %@divRatioChiSquare, divSPDParams_nci2(10e-3, 1e-3, 'gaussian'); ...
    @divRatioChiSquare, divSPDParams_I('int_exp'); ...
    %@divPD, []; ...
    %@divISF, []; ...
    %@divHilbertian, divHilbertianParams('Hellinger', 'default', 10e-3); ...
};

[p, power, dist, d12] = evaluateExperiment(spikeTrains1, spikeTrains2, M, 0.05, true, divMeasures);

% evaluateExperiment(spikeTrains1, spikeTrains2, M);
% vim:ts=8:sts=4:sw=4
