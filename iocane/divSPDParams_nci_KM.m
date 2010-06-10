function [params] = divSPDParams_nci_KM(tau, sigma)
% params = divSPDParams_nci_KM(tau, sigma)
% 
% Returns the nonlinear Cross-Intensity kernel (definition 1 in Paiva 2009).
% exp( \int (\lambda1(t) - \lambda2(t))^2 dt )
% where lambda is the smoothed spike train.
% The L2 distance of the smoothed spike trains can be computed using mCI kernel
%
% Input:
%   tau: time constant for the low pass filtering to smooth spike trains
%   sigma: kernel size for the (outter) Gaussian
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

function [Kxx, Kxy, Kyy] = k(spikeTrainsX, xIdx, spikeTrainsY, yIdx)
    Nx = length(xIdx); Ny = length(yIdx);
    Kxx = zeros(Nx, Nx); Kxy = zeros(Nx, Ny); Kyy = zeros(Ny, Ny);

    function [V] = subKernel(st1, st2) % THIS IS mCI (a.k.a. CIP)
	L1 = length(st1); L2 = length(st2);
	V = 0;
	if L1 == 0 || L2 == 0
	    return;
    end

	for kk1 = 1:L1
	    for kk2 = 1:L2
		V = V + exp(-abs(st1(kk1) - st2(kk2))/tau);
	    end
    end
    end

    for k1 = 1:Nx
	st1 = spikeTrainsX.data{xIdx(k1)};
	for k2 = 1:Nx
	    st2 = spikeTrainsX.data{xIdx(k2)};
	    Kxx(k1, k2) = subKernel(st1, st2);
	end
    end

    for k1 = 1:Nx
	st1 = spikeTrainsX.data{xIdx(k1)};
	for k2 = 1:Ny
	    st2 = spikeTrainsY.data{yIdx(k2)};
	    Kxy(k1, k2) = subKernel(st1, st2);
	end
    end

    for k1 = 1:Ny
	st1 = spikeTrainsY.data{yIdx(k1)};
	for k2 = 1:Ny
	    st2 = spikeTrainsY.data{yIdx(k2)};
	    Kyy(k1, k2) = subKernel(st1, st2);
	end
    end

    K = [Kxx, Kxy; Kxy', Kyy];

    % make the pairwise distance matrix from mCI kernel matrix K
    g = diag(K);
    N = Nx + Ny;
    D2 = g(:) * ones(1, N) + ones(N, 1) * g(:)' - 2 * K;

    % transform the distance matrix via Gaussian kernel
    K = exp(-D2/2/sigma^2);

    % decompose back
    Kxx = K(1:Nx, 1:Nx);
    Kyy = K(Nx+1:end, Nx+1:end);
    Kxy = K(Nx+1:end, 1:Nx);
end

params.sigma = sigma;
params.tau = tau;
params.kernel = @k;

end

% vim:ts=8:sts=4:sw=4
