function requireThirdParty(libraries)
% Check third party libraries
% requireThirdParty(libraries)
%
% The statistic is the time of the first spike in each spike train.
% If the spike train has no spikes, the first time is set to twice the duration
% 
% Input:
%   libraries: (cell of strings or string) name of the third party library
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

persistent isAddedPath loadedLibraries

if isempty(isAddedPath)
    pathstr = fileparts(mfilename('fullpath'));
    addpath(genpath(pathstr)); % add the absolute directory
    isAddedPath = true;
end

if isempty(loadedLibraries)
    loadedLibraries = {};
end

if iscell(libraries)
    for k = 1:length(libraries)
	requireThirdPartyLibrary(libraries{k});
    end
    return
else
    requireThirdPartyLibrary(libraries);
end


end

function requireThirdPartyLibrary(libraryName)

persistent loadedLibraries 

for k = 1:length(loadedLibraries)
    if strcmp(loadedLibraries{k}, libraryName)
	return;
    end
end

libfile = which(libraryName);
if ~isempty(libfile)
    loadedLibraries{length(loadedLibraries)+1} = libraryName;
    return;
end

fprintf('library [%s] is missing. Download and put it in the thirdparty directory\n', libraryName);
fprintf('See thirdparty/README.txt for more information\n');

switch(libraryName)
case 'sskernel'
    web('http://www.mathworks.com/matlabcentral/fileexchange/24959');
    error('sskernel is a BSD-copyrighted library, downloadable at http://www.mathworks.com/matlabcentral/fileexchange/24959');
    % urlread
otherwise
    error('Unknown library');
end

end
% vim:ts=8:sts=4:sw=4
