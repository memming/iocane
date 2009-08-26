function [params] = divHilbertianParams(dist2Name, kernelSizeName, sigmaOne)
% Generates parameters for the Hilbertian metric for point processes.
% params = divHilbertianParams(dist2Name, kernelSizeName, sigmaOne)
% 
% Input:
%   dist2Name: (string) Name of the 1/2-homogeneous metric (see [1])
%              Valid values: MSC, JS, Total-variation, Hellinger
%              Default value is Hellinger
%   kernelSizeName: (string) Name of the kernel size scaling method
%              as number of samples and dimension grow
%              Valid values: silverman, default
%   sigmaOne: (1) the kernel size for the Parzen estimator at dimension 1
%              Default value is 5 ms.
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

if nargin > 2
    sigma1 = sigmaOne;
else
    sigma1 = 5e-3;
end
params.sigma1 = sigma1;

if nargin > 1
switch lower(kernelSizeName)
case {'silverman'}
    params.kernelSizeHandle = @silvermanKernelSizeScaler;
case {'default'}
    params.kernelSizeHandle = @defaultKernelSizeScaler;
otherwise
    error('Unknown kernel size name (supported: silverman or default');
end
else
    params.kernelSizeHandle = @defaultKernelSizeScaler;
end

if nargin > 0
    params.dist2Handle = dist2HandleFactory(dist2Name);
else
    params.dist2Handle = dist2HandleFactory('Hellinger');
end

params.isSampleOnly = true;

    function [sigma] = defaultKernelSizeScaler(dimension, n)
    % [sigma] = defaultKernelSizeScaler(dimension, n)
    % Computes the kernel size for target dimension given the number of samples
    % Input:
    %   dimension: (natural number) Demension of the data
    %   n: (natural number) Total number of samples
    % External variable scope:
    %   sigma1: (double) Kernel size for one sample in one dimension
    % Output:
    %   sigma: (double) Kernel size
    %
    % Excerpt from [1]
    % we used a spherical normal distribution with
    % $\sigma_n = (\card{X^{(n)}})^{-\frac{1}{5}} \cdot \sqrt{n} \cdot \sigma_1$
    % where $n \neq 1$ is the dimension and $\sigma_1 = 20 \mbox{ms}$
    % is the bandwidth for a single sample in one dimension.

	sigma = sqrt(dimension) * sigma1 * n^-.2;
    end

    function [sigma] = silvermanKernelSizeScaler(dimension, n)
    % [sigma] = silvermanKernelSizeScaler(dimension, n)
    % Uses Silverman's rule 1/N^(-1/(dimension+4))
    % Computes the kernel size for target dimension given the number of samples
    % Input:
    %   dimension: (natural number) Demension of the data
    %   n: (natural number) Total number of samples
    % External variable scope:
    %   sigma1: (double) Kernel size for one sample in one dimension
    % Output:
    %   sigma: (double) Kernel size

	sigma = sigma1 * n^(-1/(4+dimension));
    end
end
% vim:ts=8:sts=4:sw=4
