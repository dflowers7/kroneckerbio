function obs = observationLinearWeightedSumOfSquares(outputlist, timelist, sd, name)
%observationLinearWeightedSumOfSquares Create an observation scheme for a
%   list of points, whose uncertainty is a Gaussian distribution with a
%   mean of the true value and standard deviation that is an arbitrary
%   function of the true value.
%
%   obs = observationLinearWeightedSumOfSquares(outputlist, timelist, sd, 
%                                               name)
%
%   Inputs
%   outputlist: [ positive integer vector n ]
%       The indexes of the outputs for each measurement.
%   timelist: [ nonnegative vector n ]
%       The times for each measurement
%   sd: [ handle @(t,i,y) returns positive ]
%       The standard deviation function that takes the time, the index, and
%       the true value and returns the standard deviation of the
%       measurement.
%   name: [ string ]
%       An arbitrary name for the observation scheme
%
%   Outputs
%   obs: [ observation scheme structure ]

% (c) 2015 David R Hagen
% This work is released under the MIT license.

% Clean up inputs
if nargin < 4
    name = '';
end

% Standardize inputs
outputlist = vec(outputlist);
timelist = vec(timelist);

% Find unique timelist
discrete_times = row(unique(timelist));

if isempty(name)
    name = 'WeightedSumOfSquares';
end

n = numel(outputlist);
assert(numel(timelist) == n, 'KroneckerBio:objectiveWeightedSumOfSquares:timelist', 'Input "timelist" must be a vector length of "outputlist".')

obs.Type = 'Observation.Data.LinearWeightedSumOfSquares';
obs.Name = name;
obs.Complex = false;

obs.tF = max([0;timelist]);
obs.DiscreteTimes = discrete_times;

obs.Simulation = @simulation;
obs.Objective  = @objective;

obs.F = @(sol)F(outputlist, timelist, discrete_times, sd, sol);

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        y_all = int.y;
        
        ybar = zeros(n,1);
        yhat = zeros(n,1);
        for i = 1:n
            ind = discrete_times == timelist(i);
            ybar(i) = y_all(outputlist(i),ind);
            yhat(i) = ybar(i) + randn * sd(timelist(i), outputlist(i), ybar(i));
        end
        
        sim.Type = 'Simulation.Data.LinearWeightedSumOfSquares';
        sim.Name = name;
        sim.int = int;
        sim.outputlist = outputlist;
        sim.timelist = timelist;
        sim.true_measurements = ybar;
        sim.measurements = yhat;
        sim.sd = sd;
    end

    function obj = objective(measurements, includeLogSdTerm)
        if nargin < 2
            includeLogSdTerm = [];
        end
        measurements = vec(measurements);
        assert(numel(measurements) == n , 'KroneckerBio:observationWeightedSumOfSquares:measurements', 'Input "measurements" must be a vector length of "outputlist"')
        obj = objectiveLinearWeightedSumOfSquares(outputlist, timelist, measurements, sd, name, includeLogSdTerm);
    end

    function obj = objectiveLinearWeightedSumOfSquares(outputlist, timelist, measurements, sd, name, includeLogSdTerm)
        % Input arguments:
        %   includeLogSdTerm [ boolean {true} ]
        %       Set to true to include a log(sd) term in G that penalizes
        %       larger values for the standard error. If the sd values are
        %       constant with respect to the outputs, the term only
        %       adds a constant value to the objective function and does
        %       not change which parameter set is the optimum. Note that
        %       this option does not change any of the information
        %       theoretical functions (likelihood and Fisher information
        %       functions).
        if nargin < 6
            includeLogSdTerm = [];
        end
            
        if isempty(includeLogSdTerm)
            includeLogSdTerm = true;
        end
        
        % Find unique timelist
        n = numel(outputlist);
        discrete_times = row(unique(timelist));
        
        % Inherit observation
        obj = observationLinearWeightedSumOfSquares(outputlist, timelist, sd, name);
        
        obj.Type = 'Objective.Data.WeightedSumOfSquares';
        obj.Continuous = false;
        
        obj.G = @G;
        obj.dGdy = @dGdy;
        obj.d2Gdy2 = @d2Gdy2;
        
        obj.err = @err;
        obj.derrdT = @derrdT;
        
        obj.dGdT = @dGdT;
        obj.d2GdT2 = @d2GdT2;
        obj.d2GdT2_approximate = @d2GdT2_approximate;
        obj.yhash_delta = @yhash_delta;
        if includeLogSdTerm
            sigma_best = zeros(numel(outputlist),1);
            for yi = 1:numel(outputlist)
                sigma_best(yi) = sd(timelist(yi), outputlist(yi), measurements(yi));
            end
            obj.Gbest = 2*log(sigma_best);
        else
            obj.Gbest = 0;
        end
        
        obj.p = @p;
        obj.logp = @logp;
        obj.dlogpdT = @dlogpdT;
        obj.d2logpdT2_approx = @d2logpdT2_approx;
        obj.F = @(sol)F(outputlist, timelist, discrete_times, sd, sol);
        
        obj = pastestruct(objectiveZero(), obj);
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Parameter fitting functions %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [val, discrete] = G(int)
            % Evaluate solution
            [ybar, sigma] = evaluate_sol(outputlist, timelist, discrete_times, sd, int);
            e = ybar - measurements;
            
            if includeLogSdTerm
                logsdterm = 2*log(sigma);
            else
                logsdterm = 0;
            end
            
            % Goal function
            val = sum(logsdterm + (e./sigma).^2);
            
            % Return discrete times as well
            discrete = discrete_times;
        end
        
        function val = err(int)
            [ybar, sigma] = evaluate_sol(outputlist, timelist, discrete_times, sd, int);
            val = (ybar - measurements)./sigma;
        end
        
        function val = dGdy(t, int)
            % dGdy = 2 * sigma^-1 * dsigmady + 2 * (y-yhat) * (sigma - (y-yhat)*dsigmady) * sigma^-3
            ny = int.ny;
            
            % Find all data points that have a time that matches t
            ind_t = find(t == timelist);
            n_current = numel(ind_t);
            
            if n_current > 0
                % Extract integration for this time point
                yt = int.y(:,int.t == t);
                
                % Extract the data points with time t
                timelist_t = timelist(ind_t);
                outputlist_t = outputlist(ind_t);
                measurements_t = measurements(ind_t);
                
                % Compute e for each datapoint that has a matching time to t
                ybar_t = zeros(n_current,1);
                sigma_t = zeros(n_current,1);
                dsigmady_t = zeros(n_current,1);
                for i = 1:n_current
                    ybar_t(i) = yt(outputlist_t(i));
                    [sigma_t(i), dsigmady_t(i)] = sd(timelist_t(i), outputlist_t(i), ybar_t(i));
                end
                
                % Gradient value
                e = ybar_t - measurements_t; % Y_
                if includeLogSdTerm
                    logsdterm = 2 ./ sigma_t .* dsigmady_t;
                else
                    logsdterm = 0;
                end
                dGdybar = logsdterm + 2 .* e .* (sigma_t - e.*dsigmady_t) ./ sigma_t.^3; % Y_
                val = accumarray(outputlist_t, dGdybar, [ny,1]); % sum the entries associated with the same output
            else
                val = zeros(ny, 1);
            end
        end
        
        function val = dGdT(int)
            
            [dybardT, sigma, dsigmady] = evaluate_grad(outputlist, timelist, discrete_times, sd, int);
            if includeLogSdTerm
                logsdterm = (2./sigma.*dsigmady).'*dybardT;
            else
                logsdterm = 0;
            end
            derrdT_i = derrdT(int);
            val = 2*err(int).'*derrdT_i + logsdterm;
            val = val(:);
            
%             val_test = zeros(size(val));
%             for i = numel(discrete_times):-1:1
%                 val_test = val_test + dGdy(discrete_times(i), int).'*dybardT(timelist == discrete_times(i),:);
%             end
            
        end
        
        function val = d2GdT2_approximate(int)
            val = derrdT(int);
            val = 2*(val.'*val);
        end
        
        function val = derrdT(int)
            
            %ybar = evaluate_sol(outputlist, timelist, discrete_times, sd, int);
            e = err(int);
            [dybardT, sigma, dsigmady] = evaluate_grad(outputlist, timelist, discrete_times, sd, int);
            val = bsxfun(@times, 1./sigma.*(1-e.*dsigmady), dybardT);
            
        end
        
        function val = d2Gdy2(t, int)
            % d2Gdy2dy1 = 2 * sigma^-2 - 4*e*sigma^-3*dsigmady1 - 4*e*sigma^-3*dsigmady2
            %             + (6*e^2*sigma^-4 - 2*sigma^-2)*dsigmady1*dsigmady2
            %             + (2*sigma^-1 - 2*e^2*sigma^-3)*d2sigmady1dy2
            ny = int.ny;
            
            % Find all data points that have a time that matches t
            ind_t = find(t == timelist);
            n_current = numel(ind_t);
            
            if n_current > 0
                % Extract integration for this time point
                yt = int.y(:,int.t == t);
                
                % Extract the data points with time t
                timelist_t = timelist(ind_t);
                outputlist_t = outputlist(ind_t);
                measurements_t = measurements(ind_t);
                
                ybar_t = zeros(n_current,1);
                sigma_t = zeros(n_current,1);
                dsigmady_t = zeros(n_current,1);
                d2sigmady2_t = zeros(n_current,1);
                for i = 1:n_current
                    ybar_t(i) = yt(outputlist_t(i));
                    [sigma_t(i), dsigmady_t(i), d2sigmady2_t(i)] = sd(timelist_t(i), outputlist_t(i), ybar_t(i));
                end
                
                if includeLogSdTerm
                    logsdterm = -2./sigma_t.^2.*dsigmady_t.^2 + 2./sigma_t.*d2sigmady2_t;
                else
                    logsdterm = 0;
                end
                
                % Curvature value
                e = ybar_t - measurements_t;
                d2Gdybar2 = 2 ./ sigma_t.^2 - 8 .* e ./ sigma_t.^3 .* dsigmady_t + ...
                    (6 .* e.^2 ./ sigma_t.^4) .* dsigmady_t.^2 + ...
                    (-2 .* e.^2 ./ sigma_t.^3) .* d2sigmady2_t ...
                    + logsdterm;
                val = accumarray(outputlist_t, d2Gdybar2, [ny,1]); % y_ % sum the entries associated with the same output
                val = spdiags(val, 0, ny, ny); % y_y
            else
                val = zeros(ny, ny);
            end
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Information theory %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Probability density function
        function val = p(sol)
            % p = tau^(-n/2) * prod(sigma)^-1 * exp(-1/2 * sum(((ybar-yhat)/sigma)^2)
            
            % Evaluate solution
            [ybar, sigma] = evaluate_sol(outputlist, timelist, discrete_times, sd, sol);
            e = ybar - measurements;
            
            val = (2*pi)^(-n/2) * prod(sigma).^-1 * exp(-1/2 * sum((e ./ sigma).^2));
        end
        
        %% Log likelihood
        function val = logp(sol)
            % logp = -n/2 * log(tau) + -sum(log(sigma)) + -1/2 * sum(((ybar-yhat)/sigma)^2)
            
            % Evaluate solution
            [ybar, sigma] = evaluate_sol(outputlist, timelist, discrete_times, sd, sol);
            e = ybar - measurements;
            
            val = -n/2 * log(2*pi) + -sum(log(sigma)) + -1/2 * sum((e ./ sigma).^2);
            % logp = -n/2 * log(2*pi) - sum(log(sigma)) - 0.5*sum(err.^2)
        end
        
        function val = dlogpdT(int)
            
            [sigma, dsigmadT] = getsigma(int);
            err_int = err(int);
            derrdT_int = derrdT(int);
            
            % dlogpdTi = -sum(1./sigma * dsigmadT, sigma) ...
            %   - 0.5*sum(2*err*derrdT, err)
            val = -(1./sigma).'*dsigmadT - err_int.'*derrdT_int;
            val = val(:);
            
        end
        
        function val = d2logpdT2_approx(int)
            
            [sigma, dsigmadT] = getsigma(int);
            d2sigmadT2_approx = getsigma2_approx(int);
            
            nT = int.nT;
            
            % d2log(sigmai)dTjdTk ~ 1./sigmai^2 * dsigmaidTj * dsigmaidTk
            %   - 1./sigmai * d2sigmaidTjdTk_approx
            temp1 = bsxfun(@times, (1./sigma), dsigmadT);
            temp1 = temp1.'*temp1;
            temp2 = -bsxfun(@times, 1./sigma, d2sigmadT2_approx);
            temp2 = sum(temp2, 1);
            temp2 = reshape(temp2, nT, nT);
            val_sigma = temp1 + temp2;
            
            derrdT_int = derrdT(int);
            val_err = -derrdT_int.'*derrdT_int;
            
            val = val_sigma + val_err;
            
        end
        
        function [y,dydT,d2ydT2] = gety(int)
            
            ny = int.ny;
            nT = int.nT;
            nt = numel(int.t);
            
            [~,timeindlist] = ismember(timelist, int.t);
            yinds = sub2ind(size(int.y), outputlist, timeindlist);
            y = int.y(yinds);
            
            if nargout > 1
                
                yTinds_start = sub2ind([ny nT nt], outputlist, ones(size(outputlist)), timeindlist);
                yTinds_end = sub2ind([ny nT nt], outputlist, repmat(nT, size(outputlist)), timeindlist);
                yTinds = zeros(n, nT);
                for i = 1:n
                    yTinds(i,:) = yTinds_start(i):ny:yTinds_end(i);
                end
                dydT = reshape(int.dydT(yTinds), n, nT);
                
                if nargout > 2
                    
                    yTTinds_start = sub2ind([ny nT nT nt], outputlist, 1, 1, timeindlist);
                    yTTinds_end = sub2ind([ny nT nT nt], outputlist, nT, nT, timeindlist);
                    yTTinds = zeros(n, nT^2);
                    for i = 1:n
                        yTTinds(i,:) = yTTinds_start(i):yTTinds_end(i);
                    end
                    d2ydT2 = reshape(int.d2ydT2(yTTinds), n, nT^2);
                    
                end
            end
        end
        
        function [sigma, dsigmadT, d2sigmadT2] = getsigma(int)
            
            nT = int.nT;
            
            [yc,sc] = deal(cell(max(nargout,1),1));
            [yc{:}] = gety(int);
            [sc{:}] = arrayfun(sd, timelist, outputlist, yc{1});
            
            if nargout > 0
                sigma = sc{1};
                
                if nargout > 1
                    dsigmady = sc{2};
                    dsigmadT = bsxfun(@times, dsigmady, yc{2});
                    
                    if nargout > 2
                        % d2sigmaidTjdTk = d2sigmaid2yi*dyidTj*dyidTk +
                        % dsigmaidyi*d2yidTjdTk
                        d2sigmady2 = sc{3};
                        dydT_times_dydT = bsxfun(@times, yc{2}, permute(yc{2}, [1 3 2]));
                        dydT_times_dydT = reshape(dydT_times_dydT, [n nT^2]);
                        d2sigmadT2 = bsxfun(@times, dsigmady, yc{3}) ...
                            + bsxfun(@times, d2sigmady2, dydT_times_dydT);
                    end
                end
            end
            
        end
        
        function d2sigmadT2_approx = getsigma2_approx(int)
            
            nT = int.nT;
            
            yc = cell(2,1);
            sc = cell(3,1);
            [yc{:}] = gety(int);
            [sc{:}] = arrayfun(sd, timelist, outputlist, yc{1});
            
            dydT = yc{2};
            dydT_times_dydT = bsxfun(@times, dydT, permute(dydT, [1 3 2]));
            dydT_times_dydT = reshape(dydT_times_dydT, [n nT^2]);
            
            d2sigmady2 = sc{3};
            d2sigmadT2_approx = bsxfun(@times, d2sigmady2, dydT_times_dydT);
                
        end
        
    end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Helper functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ybar, sigma] = evaluate_sol(outputlist, timelist, discrete_times, sd, int)
n = numel(outputlist);
y_all = int.y;

ybar = zeros(n,1);
sigma = zeros(n,1);
for i = 1:n
    ind = discrete_times == timelist(i);
    ybar(i) = y_all(outputlist(i),ind);
    sigma(i) = sd(timelist(i), outputlist(i), ybar(i));
end
end

function [dydT, sigma, dsigmady] = evaluate_grad(outputlist, timelist, discrete_times, sd, int)
n = numel(outputlist);
nT = int.nT;
ny = int.ny;

y_all = int.y;
dydT_all = int.dydT;

% Get dydT for every point
dydT = zeros(n, nT);
sigma = zeros(n,1);
dsigmady = zeros(n,1);
for i = 1:n
    ind = find(discrete_times == timelist(i), 1);
    t_i = timelist(i);
    
    % Compute this y
    y_i = y_all(outputlist(i),ind);
    
    % Compute dy/dT
    try
        dydT_temp = reshape(dydT_all(:,ind), ny, nT);
    catch ME
        rethrow(ME)
    end
    dydT(i,:) = dydT_temp(outputlist(i),:);
    
    % Compute expected V matrix
    [sigma(i), dsigmady(i)] = sd(t_i, outputlist(i), y_i);
end
end

function [sigma, dsigmadT, d2sigmadT2] = evaluate_sigma(outputlist, timelist, int)

end

%% Fisher information matrix
function val = F(outputlist, timelist, discrete_times, sd, sol)
n = numel(outputlist);

% Evaluate solution
[dydT, sigma, dsigmady] = evaluate_grad(outputlist, timelist, discrete_times, sd, sol);
dsigmadT = diag(dsigmady)*dydT;
dCovdsigma = 2*diag(sigma);
dCovdT = dCovdsigma*dsigmadT;
nT = size(dydT,2);

% Assemble variance matrix
V = spdiags(sigma.^2,0,n,n);

% Fisher information matrix
covarTerm = zeros(nT);
diagCovInv = diag(sigma.^-2);
for m = 1:nT
    dCovdT_m = diag(dCovdT(:,m));
    for n = 1:nT
        dCovdT_n = diag(dCovdT(:,n));
        covarTerm(m,n) = 1/2*trace(diagCovInv*dCovdT_m*diagCovInv*dCovdT_n);
    end
end
val = dydT.' * (V \ dydT) + covarTerm;
val = symmat(val);
end
