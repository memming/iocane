function div = sdivCM(data1, data2, params)
% Cramer-von-Miese test statistic via empirical cdf estimator
% div = sdivCM(data1, data2, params)
%
% If any of the data is empty, the divergence is set to 1.
% 
% Input:
%   data1, data2: (Nx1, Mx1) 2 sets of real values for comparison
%   params: (struct) not used
% Output:
%   div: (1) divergence value
%
% See also: divCountParams, divISIParams, divTTFSParams
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

if isempty(data1) || isempty(data2)
    div = 1;
    return
end

% estimate empirical distributions
edges = [-inf; sort([data1(:);data2(:)]); inf];
pmf1 = (histc(data1, edges) / length(data1));
pmf2 = (histc(data2, edges) / length(data2));

div = sum((pmf1 - pmf2).^2);

% vim:ts=8:sts=4:sw=4
