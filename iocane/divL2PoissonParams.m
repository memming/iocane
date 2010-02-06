function [params] = divL2PoissonParams(dt, mode, kernelSize)
% Generates parameters for the div based on marginal intensity function
% params = divL2PoissonParams(dt, mode, kernelSize)
%
% The intensity function is smoothed with a kernel (Gaussian).
%
% Input:
%   dt: (1) time bin size for the computation of the difference (around 1e-3)
%   mode: (string) type of smoothing kernel size selection method
%         Valid values: fixed, optimal, hist
%         Fixed uses the third argument. Optimal doesn't.
%   kernelSize: (1) fixed kernel size for smoothing
%
% Output:
%   params: (struct) ready to use for divL2Poisson
%
% See also: divL2Poisson, sskernel
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

params.dt = dt;
params.mode = mode;

switch(lower(mode))
case {'hist'}
    %params.kernelSizeHandle = @(x)(kernelSize);
    params.estimateMarginalIntensity = @(tr,spks)(histc(flattenCell(spks.data),tr));
case {'optimal', 'sskernel'}
    % requires sskernel from Shimazaki
    requireThirdParty('sskernel');
    %params.kernelSizeHandle = @(x)(sskernel(x));
    params.estimateMarginalIntensity = @(tr,spks)(smoothedEstimator(tr,spks,sskernel(spks)));
case {'fixed'}
    if nargin > 2
	params.kernelSizeHandle = @(x)(kernelSize);
    else
	error('kernel size is required');
    end
    params.estimateMarginalIntensity = @(tr,spks)(smoothedEstimator(tr,spks,kernelSize));
otherwise
    error('Unknown mode');
end
end % end function

function lambda = smoothedEstimator(tr, spikeTrains, sigma)
    allSpikes = flattenCell(spikeTrains.data);
    if isempty(allSpikes)
	lambda = zeros(size(tr));
	return
    end
    %sigma = params.kernelSizeHandle(allSpikes(:));
    lambda = ksdensity(allSpikes(:), tr, 'width', sigma);
    lambda = lambda * numel(allSpikes) / spikeTrains.N / spikeTrains.duration;
end
% vim:ts=8:sts=4:sw=4
