function [dist2Handle] = dist2HandleFactory(name)
% Generates the distance handle object on postivie reals for Hilbertian metric
% dist2Handle = dist2HandleFactory(name)
% Input:
%   name: (string) Name of the 1/2-homogeneous metric
%          Valid values: MSC, JS, Total-variation, Hellinger
% Output:
%   sigma: (@) function handle for the metric
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

switch lower(name)
case {'mean square difference', 'msd', 'msc', 'chi-square'}
    dist2Handle = @(p,q)( (p-q)^2 / (p+q) );
case {'hellinger'} % Hellinger
    dist2Handle = @(p,q)( (sqrt(p) - sqrt(q))^2 );
case {'js', 'jensen-shannon'}
    dist2Handle = @(p,q)( 1/log(2) * ...
	(p * log(2 * p/ (p+q)) + q * log(2 * q/ (p+q))) );
case {'total variation'}
    dist2Handle = @(p,q)( abs(p - q) );
otherwise
    error('Unknown R+ metric [%s]', name);
end
% vim:ts=8:sts=4:sw=4
