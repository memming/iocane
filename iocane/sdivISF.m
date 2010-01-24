function div = sdivISF(X, Y, p1, p2)
% Integration of square distance between survival function based divergence.
% div = sdivICDF(X, Y, p1, p2)
%
% Input
%   x, y: (NxM) N vectors of M dimension (column vectors)
%   p1, p2: (1/optional) probability of X, Y happening..when merging for divICDF
% Output:
%   div: \int (F_{X}(x) - G_{Y}(x))^2 dx where F and G are 
%	empirical survival functions.
%
% Original idea f - Sohan Seth
%
% Reference
% [1] Ting Chen, Baba Vemuri, Anand Rangarajan, Stephan Eisenschenk.
%   Group-Wise Point-Set Registration Using a Novel CDF-Based Havrda-Charvát
%   Divergence. International Journal of Computer Vision, Vol. 86, No. 1.
%   (1 January 2010), pp. 111-124.
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

if nargin < 1 | nargin > 4
    error('Invalid number of input arguments')
end

if nargin < 3
    p1 = 1;
    p2 = 1;
end

[nx, d1] = size(X); [ny, d2] = size(Y);

if d1 ~= d2
    error('X-Y dimension mismatch')
end

d = d1;

Kxx = ones(nx,nx); Kyx = ones(ny,nx); Kyy = ones(ny,ny);

for countDim = 1:d
    Kxx = Kxx .* min(repmat(X(:,countDim),1,nx), repmat((X(:,countDim))',nx,1));
    Kyx = Kyx .* min(repmat(Y(:,countDim),1,nx), repmat((X(:,countDim))',ny,1));
    Kyy = Kyy .* min(repmat(Y(:,countDim),1,ny), repmat((Y(:,countDim))',ny,1));
end

if isempty(Kxx); Kxx = 0; end
if isempty(Kyx); Kyx = 0; end
if isempty(Kyy); Kyy = 0; end

ox = ones(1, nx); oy = ones(1, ny);

sKxx = ox * Kxx * ox';
sKyy = oy * Kyy * oy';
sKyx = oy * Kyx * ox';

div = p1^2 * sKxx + p2^2 * sKyy - 2 * p1 * p2 * sKyx;
% vim:ts=8:sts=4:sw=4
