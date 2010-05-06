function [params] = divSPDParams_I(kappa, sigma)
% Kernel using L2 distance between realizations of counting processes
% [params] = divSPDParams_I
% 
% Input:
%   kappa: (string) kappa type - identity, exp_int, or int_exp
%   sigma: (1) kernel size for kappa. If sigma = 0, median of the data before
%              the exponential operation will be used as sigma^2
%
% Output:
%   params: (struct) ready to use for divSPD
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

% if nargin > 0
%     error('f has to be a function');
% end

if nargin < 2
    sigma = 2;
end
sigma2 = sigma^2;

function [kk] = k(spikeTrains1, k1, spikeTrains2, k2)
    st1 = spikeTrains1.data{k1};
    l1 = length(st1) + 1; % empty spike train
    st2 = spikeTrains2.data{k2};
    [z, cls] = sort([st1(:); 0; st2(:); 0]);
    if z(end) > spikeTrains1.duration
	error('There are more spikes after the end?');
    end
    %
    zz = ones(size(z));
    zz(cls > l1) = -1;
    zz = cumsum(zz);
    % alternatively
    % z1 = cumsum(cls <= l1);
    % z2 = cumsum(cls > l1);
    % zz = (z1 - z2);
    zz = zz .^2;
    switch(kappa)
    case {'identity'}
	kk = sum(diff([z; spikeTrains1.duration]) .* zz);
    case {'exp_int'}
	zzz = sum(diff([z; spikeTrains1.duration]) .* zz);
	if sigma2 == 0
	    median_sigma2 = median(zzz);
	    kk = exp(-zzz/median_sigma2);
	else
	    kk = exp(-zzz/sigma2);
	end
    case {'int_exp'}
	if sigma2 == 0
	    median_sigma2 = median(zz);
	    zz = exp(-zz/median_sigma2);
	else
	    zz = exp(-zz/sigma2);
	end
	kk = sum(diff([z; spikeTrains1.duration]) .* zz);
    otherwise
	error('invalid kernel!');
    end
end

params.kernel = @k;
params.kappa = kappa;
params.sigma = sigma;

end

% vim:ts=8:sts=4:sw=4
