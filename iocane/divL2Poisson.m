function [div] = divL2Poisson(spikeTrains1, spikeTrains2, params)
% Divergence between the estimated count distributions
% div = divL2Poisson(spikeTrains1, spikeTrains2, params)
% 
% Input:
%   spikeTrains1, spikeTrains2: (struct) 2 sets of spike trains for comparison
%   params: (struct) see divL2PoissonParams
% Output:
%   div: (1) divergence value
%
% REQUIRES: stat toolbox (ksdensity) - not any more (see divL2PoissonParams)
% See also: divL2PoissonParams
%
% This is only a dissimilarity and not a divergence.
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

% params.dt : time step size
tr = 0:params.dt:max(spikeTrains1.duration, spikeTrains2.duration);

% estimate the rate function
lambda1 = params.estimateMarginalIntensity(tr, spikeTrains1);
lambda2 = params.estimateMarginalIntensity(tr, spikeTrains2);

div = sum((lambda1 - lambda2).^2) * params.dt;

end
% vim:ts=8:sts=4:sw=4
