% This file contains the mapping from simulation number to the actual
% conditions. Each simulation runs in this directory loads this file first.
% To run a new simulation condition, one should add new numbers in the end of
% this file.
% The result of simulation will be written into a MAT file corresponding to
% the simulation index number. NEVER CHANGE THE NUMBERING only add new ones.
%
% $Id$
% Copyright 2009 Memming. All rights reserved.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each simulation condition contains the following
%   dataHandle: the function handle that generates independent trials
%	If it is a single handle it should return two spiketrain sets,
%	if it is an array of two handles, each will return a spiketrain
%   nRange: range for number of samples
%   M: independent Monte Carlo runs (each run needs independent trial sets)
%   divHandle, divParams: divergence measure and parameters to be tested.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simCond{1}.randseed = 20091206;
simCond{1}.dataHandle1 = @genSerialCorr;
simCond{1}.dataHandleParams1 = struct('T', 100e-3, 'mISI', 25e-3, 'urISI', 5e-3, 'type', 'uncorrelated');
simCond{1}.dataHandle2 = @genSerialCorr;
simCond{1}.dataHandleParams2 = struct('T', 100e-3, 'mISI', 25e-3, 'urISI', 5e-3, 'type', 'correlated');
simCond{1}.nRange = 2.^[3:11];
simCond{1}.M = 54;
simCond{1}.alpha = 0.05;
simCond{1}.divHandle = @divPhi;
simCond{1}.divParams = divPhiParams('Hellinger', 'default', 10e-3);

divMeasures = {...
    @divMeanFiringRate, []; ...
    @divFF, []; ...
    @divISI, divISIParams('KS'); ...
    @divCount, divCountParams('KS'); ...
    @divL2Poisson, divL2PoissonParams(1e-3, 'fixed', 100e-3); ...
    @divL2Poisson, divL2PoissonParams(1e-3, 'fixed', 10e-3); ...
    @divTTFS, divTTFSParams('KS'); ...
    @divHilbertian, divHilbertianParams('Hellinger', 'default', 100e-3); ...
    @divHilbertian, divHilbertianParams('Hellinger', 'default', 10e-3); ...
    @divPhi, divPhiParams('Hellinger', 'default', 100e-3); ...
    @divPhi, divPhiParams('Hellinger', 'default', 10e-3); ...
    };

for k = 1:size(divMeasures, 1)
    simCond{k+1} = simCond{1};
    simCond{k+1}.divHandle = divMeasures{k,1};
    simCond{k+1}.divParams = divMeasures{k,2};
end

%{
simCond{2}.dataHandle2 = @genHomogeneousPoisson;
simCond{2}.dataHandleParams2 = struct('lambda', 1.8, 'tOffset', 25e-3, 'duration', 50e-3);
simCond{3}.dataHandle1 = @genTwoAPex;
simCond{3}.dataHandleParams1 = struct('type', 'correlated');
simCond{3}.dataHandle2 = @genTwoAPex;
simCond{3}.dataHandleParams2 = struct('type', 'uncorrelated');
%}

% vim:ts=8:sts=4:sw=4
