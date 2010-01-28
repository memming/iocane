function [params] = divHilbertianParams(dist2Name, kernelSizeName, sigmaOne, lossyP)
% Generates parameters for the Hilbertian metric for point processes.
% params = divHilbertianParams(dist2Name, kernelSizeName, sigmaOne, lossyP)
% 
% Input:
%   dist2Name: (string) Name of the 1/2-homogeneous metric (see [1])
%              Valid values: MSC, JS, Total-variation, Hellinger
%              Default value is Hellinger
%   kernelSizeName: (string) Name of the kernel size scaling method
%              as number of samples and dimension grow
%              Valid values: silverman, default, knn
%   sigmaOne: (1) the kernel size for the Parzen estimator at dimension 1
%              Default value is 5 ms. If knn is used, then this is the 'k'.
%   lossyP: (1/optional) smoothing due to lossy APs to double the number of
%	    samples used to estimate the divergence.
%   lossyP: (1/optional) smoothing due to lossy APs to double the number of
%	    samples used to estimate the divergence.
% Output:
%   params: (struct) ready to use for divHilbertian
%
% See also: divHilbertian
%
% References
% [1] Matthias Hein, Olivier Bousquet. "Hilbertian Metrics and Positive Definite
%   Kernels on Probability Measures" In AISTATS (2005)
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

if nargin > 3
    params.lossyP = lossyP;
end

if nargin > 2
    sigma1 = sigmaOne;
else
    sigma1 = 5e-3;
end

if nargin > 1
    if strcmp(lower(kernelSizeName), 'knn')
	if sigmaOne > 0 && (sigmaOne - round(sigmaOne)) == 0
	    params.kNearest = sigmaOne;
	else
	    warning('kNN adaptive estimator requires integer as k');
	    params.kNearest = 1;
	end
    else
	params.kernelSizeHandle = kernelSizeScalerFactory(kernelSizeName);
	params.sigma1 = sigma1;
    end
else
    params.kernelSizeHandle = kernelSizeScalerFactory();
    params.sigma1 = sigma1;
end

if nargin > 0
    params.dist2Handle = dist2HandleFactory(dist2Name);
else
    params.dist2Handle = dist2HandleFactory('Hellinger');
end

params.isSampleOnly = true;

% vim:ts=8:sts=4:sw=4
