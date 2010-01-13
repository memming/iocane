function div = sdivCDF(X, Y, p, p1, p2)
% Divergence based on Lp norm of empirical CDF difference
% div = sdivCDF(X, Y, p)
%
% Input
%   x, y: (NxM) N vectors of M dimension (column vectors)
%   p: (1) p for the Lp norm (0 < p < Inf)
%   p1, p2: (1/optional) probability of X, Y happening...when merging for divCDF
% Output:
%   div^p: \frac{1}{n} \sum_{i=1}^n |F_{X}(x_i) - G_{Y}(x_i)|^p
%
% Original author - Sohan Seth; Date - 12.08.2009
%
% $Id$
%
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

if nargin < 1 | nargin > 5
    error('Invalid number of input arguments')
end

if nargin == 2
    p = 1;
end

if nargin < 4
    p1 = 1;
    p2 = 1;
end

[nx, d1] = size(X); [ny, d2] = size(Y);

if d1 ~= d2
    error('X-Y dimension mismatch')
end

d = d1;
Kxx = ones(nx,nx); Kyx = ones(ny,nx); Kxy = ones(nx,ny); Kyy = ones(ny,ny);

for countDim = 1:d
    Kxx = Kxx .* (repmat(X(:,countDim),1,nx) <= repmat((X(:,countDim))',nx,1));
    Kyx = Kyx .* (repmat(Y(:,countDim),1,nx) <= repmat((X(:,countDim))',ny,1));
    Kxy = Kxy .* (repmat(X(:,countDim),1,ny) <= repmat((Y(:,countDim))',nx,1));
    Kyy = Kyy .* (repmat(Y(:,countDim),1,ny) <= repmat((Y(:,countDim))',ny,1));
end

if isempty(Kxx); Kxx = 0; end
if isempty(Kxy); Kxy = 0; end
if isempty(Kyx); Kyx = 0; end
if isempty(Kyy); Kyy = 0; end

if isinf(p)
    div = 0.5 * max(abs(p1 * mean(Kxx) - p2 * mean(Kyx))) ...
	+ max(abs(p1 * mean(Kxy) - p2 * mean(Kyy)));
else
    div = 0.5 * ((mean((abs(p1 * mean(Kxx) - p2 * mean(Kyx))).^p))^(1/p) ...
	+ (mean((abs(p1 * mean(Kxy) - p2 * mean(Kyy))).^p))^(1/p));
end

% vim:ts=8:sts=4:sw=4
