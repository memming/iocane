function [d2, d2_fg, d2_gf] = fppPhiDivRatio(fppm1, fppm2, phiHandle)
% d2 = fppHellingerRatio(fppm1, fppm2)
% Hellinger distance estimated from ratio
% NOTE: This function is not symmetric on its arguments.
% \int (\sqrt{f2 / f1} - 1)^2 dF1
%
% Input:
%   fppm1, fppm2: (struct) FPPM (see estimateFPPM)
%
% Output:
%   d2: computed squared divergence
%
% See also fppHilbertianMetricSamples
%
% $Id$
% Copyright 2009 Memming. All rights reserved.

dTemp1 = zeros(fppm1.spikeTrains.N, 1);
for k = 1:fppm1.spikeTrains.N
    st = fppm1.spikeTrains.data{k};
    j1 = likelihoodFPPM(fppm1, st);
    j2 = likelihoodFPPM(fppm2, st);
    dTemp1(k) = phiHandle(j2/j1);
end
d2_fg = mean(dTemp1);
clear dTemp1;

dTemp2 = zeros(fppm2.spikeTrains.N, 1);
for k = 1:fppm2.spikeTrains.N
    st = fppm2.spikeTrains.data{k};
    j1 = likelihoodFPPM(fppm1, st);
    j2 = likelihoodFPPM(fppm2, st);
    dTemp2(k) = phiHandle(j1/j2);
end
d2_gf = mean(dTemp2);
clear dTemp2;

d2 = mean([d2_fg, d2_gf]);
% vim:ts=8:sts=4:sw=4

