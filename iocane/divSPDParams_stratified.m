function [params] = divSPDParams_stratified(kernelName, kernelSizeName, sigma);
% Generates parameters for the SPD-divergence for point processes.
% [params] = divSPDParams_stratified(kernelName, kernelSizeName, sigma);
% 
% Input:
%   kernelName: (string) The kernel type. Must be strictly positive definite
%	       to be a divergence. If not only partial statistics is used.
%              Valid values: Guassian, constant
%   kernelSizeName: (string) The kernel size scaler type.
%              Default
%   sigma: (1/optional) the kernel size.
%              Default value is 5 ms.
% Output:
%   params: (struct) ready to use for divSPD
%
% See also: divSPD_stratified
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

if nargin < 3
    sigma = 5e-3;
end
params.sigma = sigma;

if isa(kernelName, 'function_handle')
    params.kernel = kernelName;
else
    switch(lower(kernelName))
    case {'gaussian'}
	params.kernel = @(x,y,N)(exp(-sum((x-y).^2)/2/sigma/sqrt(N)));
	%params.kernel = @(x,y,N)(exp(-sum((x-y).^2)));
    case {'constant'}
	params.kernel = @(x,y,N)(1);
    otherwise
	error('Unknown kernel name [%s]', kernelName);
    end
end

% vim:ts=8:sts=4:sw=4
