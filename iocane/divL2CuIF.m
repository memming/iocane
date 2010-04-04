function [div] = divL2CuIF(spikeTrains1, spikeTrains2, params)
% L2 distance of (marginal) cumulative intensity functions
% [div] = divL2CuIF(spikeTrains1, spikeTrains2, params)
% 
% Input:
%   spikeTrains1, spikeTrains2: (struct) 2 sets of spike trains for comparison
%   params: not used
% Output:
%   div: (1) divergence value
%
% This is only a dissimilarity and not a divergence.
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

M1 = sum(cellfun('length', spikeTrains1.data))+1;
M2 = sum(cellfun('length', spikeTrains2.data))+1;
st = zeros(M1+M2, 1);
last = 2;
for k = 1:spikeTrains1.N
    s = spikeTrains1.data{k};
    l = length(s);
    st(last:last+l-1) = s;
    last = last + l;
end
last = last + 1;
for k = 1:spikeTrains2.N
    s = spikeTrains2.data{k};
    l = length(s);
    st(last:last+l-1) = s;
    last = last + l;
end

[z, cls] = sort(st);
zz = ones(size(z)) / spikeTrains1.N;
zz(cls > M1) = -1 / spikeTrains2.N;
zz = cumsum(zz);
zz = zz .^2;
div = sum(diff([z; spikeTrains1.duration]) .* zz);
