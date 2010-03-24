function [params] = divSPDParams_fgh(f, g, h, dt)
% Generates parameters for the SPD-divergence for point processes.
% [params] = divSPDParams_fgh(f, g, h, dt)
% 
% Input:
%   f: (function_handle) CAUSAL spike train smoothing kernel
%   g: (function_handle) translation invariant positive definite function
%   h: (function_handle) strictly positive weighting function
%   dt: integration (time) step size
%
% Output:
%   params: (struct) ready to use for divSPD
%
% WARNING: This is for prototyping. The SPD kernel provided in this setting
%          is EXTREMELY slow.
%
% See also: divSPD
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

%params.kernel = @(s1,k1,s2,k2)(k(s1,k1,s2,k2,dt));

if ~isa(f, 'function_handle') || ~isa(g, 'function_handle') || ~isa(h, 'function_handle')
    error('f,g,h has to be a function');
end

if ~isnumeric(dt) || (dt <= 0)
    error('dt has to be a positive real number');
end

function [kk] = k(spikeTrains1, k1, spikeTrains2, k2)
    tr = 0:dt:spikeTrains1.duration;
    nNeuron = size(spikeTrains1.data, 2);
    kk = 1;
    for kkk = 1:nNeuron
	st1 = spikeTrains1.data{k1, kkk}; st2 = spikeTrains2.data{k2, kkk};
	x1 = zeros(size(tr)); x2 = zeros(size(tr));
	for k = 1:length(st1)
	    tidx = ceil(st1(k)/dt);
	    x1(tidx:end) = x1(tidx:end) + f(tr(tidx:end) - st1(k));
	end
	for k = 2:length(st2)
	    tidx = ceil(st2(k)/dt);
	    x2(tidx:end) = x2(tidx:end) + f(tr(tidx:end) - st2(k));
	end
	kk = kk * (sum(h(tr) .* g(x1 - x2)) * dt);
    end
end

params.kernel = @k;

end

% vim:ts=8:sts=4:sw=4
