function [params] = divSPDParams_nci2_KM(pwidth, ksize, kern_str)
% params = divSPDParams_nci2_KM(pwidth, ksize, kern_str)
% 
% Returns the nonlinear Cross-Intensity kernel (definition 2 in Paiva 2009).
% \int K(\lambda1(t) - \lambda2(t)) dt
% where lambda is the smoothed spike train with a rectangular function
%
% Input:
%  PWIDTH: Duration of the rectangular pulse smoothing function for
%          intensity estimation (sec).
%   KSIZE: Bandwidth parameter of nonlinear kernel (spk/s).
%  KERNEL: [optional] String describing which kernel to use in the
%          evaluation of the inner product. Default: 'gaussian'.
%          Other known values are: 'laplacian', 'triangular', and 'rectwin'.
%
% Output:
%   params: (struct) ready to use for divSPD
%
% Based on code by Antonio Paiva. (Feb 2008)
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

if (exist('kern_str') ~= 1)
    kern_str = 'gaussian';
end

switch lower(kern_str)
    case 'gaussian'
        K = @(z) exp(-(z.^2) ./ (2*(ksize^2)));
    case 'laplacian'
        K = @(z) exp(-abs(z) ./ ksize);
    case 'triangular'
        K = @(z) ((abs(z) < ksize) .* (1 - abs(z)/ksize));
    case 'rectwin'
        K = @(z) (abs(z) <= ksize);
    case 'cauchy'
        K = @(z) (1 ./ (1 + (z.^2) / ksize^2));
otherwise
	error('Unknown kernel! Try one: ''laplacian'', ''gaussian'', ''triangular'' and ''rectwin''');
end

function [Kxx, Kxy, Kyy] = k(spikeTrainsX, xIdx, spikeTrainsY, yIdx)
    Nx = length(xIdx); Ny = length(yIdx);
    Kxx = zeros(Nx, Nx); Kxy = zeros(Nx, Ny); Kyy = zeros(Ny, Ny);
    T = max(spikeTrainsX.duration, spikeTrainsY.duration);

    function [V] = subKernel(st1, st2, T)
	L1 = length(st1); L2 = length(st2);
	if L1 == 0 && L2 == 0
        V = K(0)*T;
	    return;
	end
	st1 = st1(:)'; st2 = st2(:)';

	% times: when the difference between intensity functions changes
	[times idx] = sort([0, st1-pwidth/2, st1+pwidth/2, ...
						st2-pwidth/2, st2+pwidth/2]);
	% if the difference should increase or decrease
	incr = [0, ones(1,L1), -ones(1,L1), -ones(1,L2), ones(1,L2)];
	incr = incr(idx);
	val = cumsum(incr) / pwidth;
	times(times < 0) = 0;
	times(times > T) = T;
	dtimes = diff([times T]);
	nzidx = (dtimes ~= 0);
	V = sum(dtimes(nzidx) .* K(val(nzidx)));
    end

    for k1 = 1:Nx
	st1 = spikeTrainsX.data{xIdx(k1)};
	for k2 = 1:Nx
	    st2 = spikeTrainsX.data{xIdx(k2)};
	    Kxx(k1, k2) = subKernel(st1, st2, T);
	end
    end

    for k1 = 1:Nx
	st1 = spikeTrainsX.data{xIdx(k1)};
	for k2 = 1:Ny
	    st2 = spikeTrainsY.data{yIdx(k2)};
	    Kxy(k1, k2) = subKernel(st1, st2, T);
	end
    end

    for k1 = 1:Ny
	st1 = spikeTrainsY.data{yIdx(k1)};
	for k2 = 1:Ny
	    st2 = spikeTrainsY.data{yIdx(k2)};
	    Kyy(k1, k2) = subKernel(st1, st2, T);
	end
    end
end

params.kernel = @k;

end

% vim:ts=8:sts=4:sw=4
