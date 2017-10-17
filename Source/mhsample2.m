function xsamps = mhsample2(x0, nsamps, logpdf, logproppdf, proprnd, options)
% Input arguments
%   x0
%       Row vector indicating the starting position.
%   nsamps
%       Integer scalar indicating the number of samples to be taken.
%   logpdf
%       Function handle for a function p = logpdf(x), which, for a value of
%       x, returns the natural log of the probability density function. The
%       value can be incorrect by an additive factor, corresponding to an
%       unknown normalization factor.
%   logproppdf
%       Function handle for a function p = logproppdf(y,x), which, for a
%       proposal value y and a current value x, returns the natural log of
%       the probability density function of the proposal distribution. The
%       distribution can be incorrect by an additive factor, corresponding
%       to an unknown normalization factor.
%   proprnd
%       Function handle for a function y = proprnd(x), which generates a
%       sample y from the proposal distribution with current value x.
%   options
%       Struct array with the following fields:
%       .Verbose
%           Set to 1 to display iteration information. Set to 0 to disable.
%       .BurnIn
%           A nonnegative integer value indicating the number of initial
%           iterations to throw away.
%       .Thin
%           A positive integer value. One out of every this many samples
%           will be kept.
%
% Output arguments
%   xsamps
%       nsamps-by-nx array of samples drawn from the distribution defined
%       by logpdf.
%
% Written by David Flowers
% Copyright 2017

if nargin < 6
    options = [];
end

testing = nargin == 0;
if testing
    nvars = 1;
    
    x0 = rand(1,nvars);
    nsamps = 1e4;
    
    mu = rand(1,nvars);
    V = rand(nvars);
    V = V'*V;
    
    logpdf = @(x)log(chi2pdf(x, 3));
    logproppdf = @(y,x)log(normpdf(y, x+0.1, 3));
    proprnd = @(x)x+0.1 + 3*randn;
    
    options.Verbose = 1;
    options.BurnIn = 1e3;
    options.Thin = 2;
    options.Rng = [];
end

if isempty(options)
    options = struct;
end

opts_table = {
    'Verbose'   1
    'BurnIn'    0
    'Thin'      1
    'Rng'       []
    };
for i = 1:size(opts_table,1)
    
    fieldi = opts_table{i,1};
    vali = opts_table{i,2};
    
    % Ensure fields exist
    if ~isfield(options, fieldi)
        options.(fieldi) = [];
    end
    
    % Set defaults
    if isempty(options.(fieldi))
        options.(fieldi) = vali;
    end
end

verbose = options.Verbose;
burnin = options.BurnIn;
thin = options.Thin;

% Set the RNG seed, if one is provided
if ~isempty(options.Rng)
    rng(options.Rng);
end

% Calculate the number of samples we need
nsamps_new = burnin + nsamps*thin;

nx = numel(x0);
xsamps = nan(nsamps_new, nx);
x = x0;

athresh = log(rand(nsamps_new,1));
naccepts = 0;
for i = 1:nsamps_new
    y = proprnd(x);
    
    ly = logpdf(y);
    lx = logpdf(x);
    lpy = logproppdf(y,x);
    lpx = logproppdf(x,y);
    a = ly + lpx - lx - lpy;
    
    accept = a > athresh(i);
    
    if accept
        xsamps(i,:) = y;
        x = y;
        naccepts = naccepts + 1;
    else
        xsamps(i,:) = x;
    end
    
    if verbose > 0
        if mod(i,25) == 1
            fprintf('\n')
            fprintf('%10s%10s%15s%12s%12s%15s%15s%10s%20s\n', 'Iter', 'Accepted', 'Accept rate', 'logp(x)', 'logp(y)', 'logproppdf(x)', 'logproppdf(y)', '|x-x0|', '|log10(x./x0)|')
            fprintf([repmat('-', 1, 119) '\n'])
        end
        fprintf('%10d%10.0g%15.2g%12.3g%12.3g%15.3g%15.3g%10.3g%20.3g\n', i, accept, ...
            naccepts/i, lx, ly, lpx, lpy, norm(x-x0), norm(log10(x)-log10(x0)));
    end
    
end

% Remove burn-in samples
xsamps(1:burnin,:) = [];

% Thin samples
xsamps = xsamps(thin:thin:end, :);

if testing
    if nvars == 1
        figure;
        histogram(xsamps, 'Normalization', 'pdf')
        hold on
        fplot(@(x)exp(logpdf(x)), [min(xsamps) max(xsamps)])
        hold off
    else
    end
end

end