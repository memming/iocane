function [params] = divPhiParams(phiName, kernelSizeName, sigmaOne, lossyP)
% Generates parameters for the Phi-divergence for point processes.
% params = divPhiParams(dist2Name, kernelSizeName, sigmaOne, lossyP)
% 
% Input:
%   phiName: (string) Name of the Phi-function
%   kernelSizeName: (string) Name of the kernel size scaling method
%              as number of samples and dimension grow
%              Valid values: silverman, default
%   sigmaOne: (1) the kernel size for the Parzen estimator at dimension 1
%              Default value is 5 ms.
%   lossyP: (1/optional) smoothing due to lossy APs to double the number of
%	    samples used to estimate the divergence.
% Output:
%   params: (struct) ready to use for divPhi
%
% See also: divPhi
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

if nargin > 1
    params.kernelSizeHandle = kernelSizeScalerFactory(kernelSizeName);
else
    params.kernelSizeHandle = kernelSizeScalerFactory();
end

if nargin > 2
    sigma1 = sigmaOne;
else
    sigma1 = 5e-3;
end
params.sigma1 = sigma1;

if nargin > 0
    params.phiHandle = phiHandleFactory(phiName);
else
    params.phiHandle = phiHandleFactory('Hellinger');
end

% vim:ts=8:sts=4:sw=4
