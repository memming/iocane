function simpleScan(expNum, divMeasures, nMCoveride)

NList = round(logspace(log10(10), log10(150), 10));
if nargin > 2
    nMC = nMCoveride;
else
    nMC = 5;
end

%% Fast test
%nMC = 2; NList = [15, 25];

power = zeros(nMC, length(NList), size(divMeasures,1));
rtime = zeros(nMC, length(NList), size(divMeasures,1));

for kMC = 1:nMC
    for NIdx = 1:length(NList)
	N = NList(NIdx);
	[s1, s2, p, power_simple, d, d12, rtime_simple] = ...
	    simpleExperiments(expNum, false, divMeasures, N);
	power(kMC,NIdx,:) = power_simple;
	rtime(kMC,NIdx,:) = rtime_simple;
    end
end

comp_arch = computer;
bench_result = bench(3);
[a, hname] = system('hostname');
save(['simpleScan_' num2str(expNum) '_' datestr(now,30) '.mat'], 'power', 'NList', 'nMC', 'comp_arch', 'rtime', 'bench_result');

% Plot results
color = {'b', 'k', 'r', 'g', 'c', 'y'};
marker = {'-', ':', 'o-', 'x:', '--', 'd:', '+-', '.-'};
figure; hold all;
kkk = 1;
for k = 1:size(divMeasures, 1)
    errorbar(NList, squeeze(mean(power(:, :, k),1)), squeeze(std(power(:, :, k),1)), [color{2} marker{k}]);
    %lgd{kkk} = [func2str(divMeasures{k,1}) ' ' num2str(jitterSigmaList(jitterSigmaIdx))]; kkk = kkk + 1;
end
%legend(lgd(:), 'Location', 'Northwest');
set(gca, 'XScale', 'log');
xlabel('Number of samples');

%for k = 1:7; a = squeeze(h(k, :, :)); meanh(k) = mean(a(:)); stdh(k) = std(mean(a,2)); end
%errorbar(NList, meanh, stdh)

figure; hold all;
title('Running time');
for k = 1:size(divMeasures, 1)
    errorbar(NList, squeeze(mean(rtime(:, :, k),1)), squeeze(std(rtime(:, :, k),1)), [color{2} marker{k}]);
end

% vim:ts=8:sts=4:sw=4
