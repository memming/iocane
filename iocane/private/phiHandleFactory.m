function [phiHandle] = phiHandleFactory(name)
% Generates the phi function handle (continuous on R+ and convex)
% phiHandle = phiHandleFactory(name)
% Input:
%   name: (string or function handle) Name of the phi-function
%          Valid values: Total-variation, Hellinger
% Output:
%   sigma: (@) function handle that takes the density ratio as argument
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

switch lower(name)
case {'hellinger'}
    phiHandle = @(x)(sqrt(x) - 1)^2;
case {'total variation'}
    phiHandle = @(x)(abs(x - 1));
case {'kl', 'kullbak-leibler'} % KL
    phiHandle = @(x)(x * log(x) - x + 1);
case {'mdi', 'minimum discrimination information'}
    phiHandle = @(x)(-log(x) + x - 1);
case {'chi-square', 'pearson'}
    phiHandle = @(x)(((x-1)^2)/2);
otherwise
    error('Unknown phi-divergence');
end
% vim:ts=8:sts=4:sw=4
