function spikeTrains = lossySmoothing(spikeTrains, params)

% TODO: current implementation always doubles the number of samples

N = spikeTrains.N;
% lossyAPs = @(st,p)(st(rand(size(st)) >= p));
lossyAPs = @(st)(st(rand(size(st)) >= params.lossyP));

for k = 1:N
    spikeTrains.data{N+k} = lossyAPs(spikeTrains.data{k});
end

spikeTrains.N = 2 * N;
