function [params] = divRatioChiSquareParams_I(kernelString, sigma)
% Kernels using L2 distance between realizations of counting processes
% [params] = divSPDParams_I(kernelString, sigma)
% The kernels are symmetric.
%
% Input:
%   kernelString: (string) 'identity', 'exp_int', or 'int_exp'
%   sigma: (postive value or string) kernel size (standard deviation) or
%	'median' (to choose the median of values before exp operation)
%
% 'identity' kernel: squared L2 distance between counting process.
%   sigma is ignored.
% 'exp_int' kernel: exponential of the squared L2 distance
% 'int_exp' kernel: squared L2 distance of the exponential of counting process
%
% Output:
%   params: (struct) ready to use for divRatioChiSquare
%
% See also: divRatioChiSquare

switch(kernelString)
case {'identity'}
    isIntFirst = true;
case {'exp_int'}
    isIntFirst = true;
    isExpInt = true;
case {'int_exp'}
    isIntFirst = false;
otherwise
    error('Unknown kernel type');
end

isMedianSigma = false;
if isnumeric(sigma)
    if sigma <= 0
	error('Kernel size (sigma) has to be strictly positive');
    end
else
    switch(sigma)
    case {'median'}
	if ~isExpInt
	    error('Only exp_int kernel can use median kernel size');
	end
	isMedianSigma = true;
    end
end

function [Kxx, Kxy, Kyy] = k(spikeTrainsX, xIdx, spikeTrainsY, yIdx)
    Nx = length(xIdx); Ny = length(yIdx);
    Kxx = zeros(Nx, Nx); Kxy = zeros(Nx, Ny); Kyy = zeros(Ny, Ny);
    T = max(spikeTrainsX.duration, spikeTrainsY.duration);

    function zzz = subKernel(st1, st2, T)
	% computes the kernel or part of it between pairs of spike trains
	[z, cls] = sort([st1(:); 0; st2(:); 0]);
	if z(end) > T
	    error('There are more spikes after the end?');
	end
	zz = ones(size(z));
	zz(cls > length(st1) + 1) = -1;
	zz = cumsum(zz).^2;
	if isIntFirst % 'identity', 'exp_int'
	    zzz = sum(diff([z; T]) .* zz);
	else % 'int_exp'
	    zz = exp(-zz/sigma^2);
	    zzz = sum(diff([z; T]) .* zz);
	end
    end

    % first pre-compute the counting processes
    for k1 = 1:Nx
	st1 = spikeTrainsX.data{xIdx(k1)};
	for k2 = 1:Nx
	    st2 = spikeTrainsX.data{xIdx(k2)};
	    Kxx(k1, k2) = subKernel(st1, st2, T);
	end
    end

    if isMedianSigma
	% We fix the kernel size from X and apply to Kxy and Kyy
	sigma = median(Kxx(:));
	disp(sigma)
    end

    if isExpInt % 'exp_int'
	Kxx = exp(-Kxx/sigma^2);
    end

    for k1 = 1:Nx
	st1 = spikeTrainsX.data{xIdx(k1)};
	for k2 = 1:Ny
	    st2 = spikeTrainsY.data{yIdx(k2)};
	    Kxy(k1, k2) = subKernel(st1, st2, T);
	end
    end
    
    if isExpInt % 'exp_int'
	Kxy = exp(-Kxy/sigma^2);
    end

    for k1 = 1:Ny
	st1 = spikeTrainsY.data{yIdx(k1)};
	for k2 = 1:Ny
	    st2 = spikeTrainsY.data{yIdx(k2)};
	    Kyy(k1, k2) = subKernel(st1, st2, T);
	end
    end

    if isExpInt % 'exp_int'
	Kyy = exp(-Kyy/sigma^2);
    end
end

params.kernel = @k;
params.kernelString = kernelString;
params.sigma = sigma;

end

% vim:ts=8:sts=4:sw=4
