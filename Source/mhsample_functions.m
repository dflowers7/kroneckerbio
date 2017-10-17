function [logpdf_fun, logproppdf_fun, proprnd_fun, testing_info] = mhsample_functions(m, con, obj, opts)
% Input arguments
%   m
%   con
%   obj
%   opts
%       .MinimumEigenvalueThreshold
% Output arguments
%   logpdf_fun
%   proplogpdf_fun
%   proprnd_fun
% TODO:
%   add nonlinear constraint functionality

if nargin == 0
    m = 'test';
end

testing = ischar(m) && strcmp(m, 'test');

nmemo = 5;
[ints,Ts] = deal(cell(nmemo,1));

% Fix the opts struct
if ~testing
    [m,con,obj,opts] = FixFitObjectiveOpts(m, con, obj, opts);
    
    % Fix extra options
    field_defaults = {
        'MinimumEigenvalueThreshold'    1
        };
    for ii = 1:size(field_defaults,1)
        if ~isfield(opts, field_defaults{ii,1})
            opts.(field_defaults{ii,1}) = [];
        end
        if isempty(opts.(field_defaults{ii,1}))
            opts.(field_defaults{ii,1}) = field_defaults{ii,2};
        end
    end
end

% Set the function to be used to calculate the log PDF and its
% sensitivities
logpdf_sens_fun = @logpdf_sens;

logpdf_fun = @logpdf;
logproppdf_fun = @proplogpdf;
proprnd_fun = @proprnd;

if testing
    % Define a multinormal distribution to test the sampler on
    nvar = 3;
    V = rand(nvar);
    V = V.'*V;
    detV = det(V);
    Vinv = inv(V);
    mu = rand(nvar,1);
    
    opts.Normalized = false;
    opts.LowerBound = abs(mu)/2;
    opts.UpperBound = abs(mu)*2;
    
    logpdf_sens_fun = @logpdf_sens_test;
    
    testing_info.V = V;
    testing_info.mu = mu;
    
    % Test that the logpdf and logproppdf functions return the same value
    ntests = 3;
    test_T = rand(ntests,nvar);
    for ii = 1:ntests
        p1 = logpdf_fun(test_T(ii,:));
        p2 = logproppdf_fun(test_T(ii,:), rand(1,nvar));
        assert(abs(p1 - p2) < 1e-6 && abs(p1-p2)/p1 < 1e-6, 'logpdf and logproppdf test failure!')
    end
    
    % Compare functions against rejection sampling, which works okay for the test
    % problem
    nsamps = 10000;
    reject_samps = bsxfun(@plus, (chol(V,'lower')*randn(nvar, nsamps)).', mu(:).');
    viol_lb = any(bsxfun(@lt, reject_samps, opts.LowerBound(:).'), 2);
    viol_ub = any(bsxfun(@gt, reject_samps, opts.UpperBound(:).'), 2);
    reject_samps(viol_lb | viol_ub, :) = [];
    reject_mu = mean(reject_samps, 1);
    % Bootstrap confidence intervals
    nbs = 1000;
    nsamps_surv = size(reject_samps, 1);
    bsi = randi(nsamps_surv, nsamps_surv, nbs);
    reject_bs = reject_samps(bsi(:),:);
    reject_bs = reshape(reject_bs, [nsamps_surv nbs 3]);
    reject_bs_mus = squeeze(mean(reject_bs, 2));
    reject_bs_stderrs = std(reject_bs_mus, 0, 1);
    
    % Test that the proposal sampling function is working correctly
    proprnd_samps = nan(nsamps, nvar);
    for ii = 1:nsamps
        proprnd_samps(ii,:) = proprnd_fun(rand(1,nvar));
    end
    proprnd_mu = mean(proprnd_samps, 1);
    
    proprnd_mu_err = abs(proprnd_mu - reject_mu)./reject_bs_stderrs;
    assert(max(proprnd_mu_err) < 3, 'proprnd test failure!')
    
%     mh_samples = mhsample(zeros(1,nvar), nsamps, 'logpdf', logpdf_fun, 'logproppdf', logproppdf_fun, 'proprnd', proprnd_fun);
    mh_samples = mhsample2(zeros(1,nvar), nsamps, logpdf_fun, logproppdf_fun, proprnd_fun);
    mh_samples_mu = mean(mh_samples, 1);
    mh_samples_err = abs(mh_samples_mu - reject_mu)./reject_bs_stderrs;
    assert(max(mh_samples_err) < 3, 'mhsample test failure!')
    
    fprintf('All tests passed.\n')
end


%% mhsample functions

    function p = logpdf(T0)
        % Convert to column vector
        T0 = T0(:);
        
        p = logpdf_sens_fun(T0);
    end

    function p = proplogpdf(T1, T0)
        
        % Convert to column vectors
        T1 = T1(:);
        T0 = T0(:);
        
        [p0,dp0dT,d2p0dT2] = logpdf_sens_fun(T0);
        
        if opts.Normalized
            lT1 = log(T1);
            lT0 = log(T0);
        else
            lT1 = T1;
            lT0 = T0;
        end
        
        ldT = lT1 - lT0;
        p = p0 + dp0dT(:).'*ldT + 0.5*ldT.'*d2p0dT2*ldT;
        
    end

    function T1 = proprnd(T0)
        
        % Convert to column vectors
        T0 = T0(:);
        
        [~, dp0dT, d2p0dT2_approx] = logpdf_sens_fun(T0);
        
        if opts.Normalized
            lT0 = log(T0);
            llb = log(opts.LowerBound);
            lub = log(opts.UpperBound);
        else
            lT0 = T0;
            llb = opts.LowerBound;
            lub = opts.UpperBound;
        end
        
        % The proposal distribution is a Gaussian truncated by the bounds,
        % with mean mu = T0 + dT_mu, where
        %   dp0dT = -d2p0dT2*dT_mu --> dT_mu = -d2p0dT2 \ dp0dT
        ldT_mu = -d2p0dT2_approx \ dp0dT;
        lmu = lT0 + ldT_mu;
        
        % and with variance V = -inv(d2p0dT2_approx)
        lV = -inv(d2p0dT2_approx);
                
        lT1 = mvnbndrndgibbs(lmu, lV, llb, lub, [], 1);
        
        if opts.Normalized
            T1 = exp(lT1);
        else
            T1 = lT1;
        end
        
        % Convert back to row vector
        T1 = T1(:).';
        
    end

%% Helper functions

    function [mT,conT] = update(T)
        [mT,conT] = updateAll(m, con, T, opts.UseParams, opts.UseSeeds, ...
            opts.UseInputControls, opts.UseDoseControls); 
    end

    function int = simulate(T)
        % Look for a matching stored evaluation
        for i = 1:nmemo
            if isequal(T, Ts{i})
                int = ints{i};
                return
            end
        end
        
        % If there isn't one, perform a new integration
        [mT, conT] = update(T);
        if opts.ParallelizeExperiments
            sim = ParallelizeExperiments(@SimulateSensitivity, mT, conT, obj, opts);
        else
            sim = SimulateSensitivity(mT, conT, obj, opts);
        end
        int = reshape([sim.int], size(obj));
        
        % Store the new integration, tossing the oldest one
        ints = [ints(2:nmemo); {int}];
        Ts = [Ts(2:nmemo); {T}];
    end

    function [p, dpdT, d2pdT2] = logpdf_sens(T)
        int = simulate(T);
        order = nargout - 1;
        
        [p,dpdT,d2pdT2] = deal([]);
        
        for i = 1:numel(obj)
            obji = obj(i);
            inti = int(i);
            
            if isempty(p)
                p = obji.logp(inti);
            else
                p = p + obji.logp(inti);
            end
            
            if order > 0
                if isempty(dpdT)
                    dpdT = obji.dlogpdT(inti);
                else
                    dpdT = dpdT + obji.dlogpdT(inti);
                end
                
                if order > 1
                    if isempty(d2pdT2)
                        d2pdT2 = obji.d2logpdT2_approx(inti);
                    else
                        d2pdT2 = d2pdT2 + obji.d2logpdT2_approx(inti);
                    end
                end
                
            end
        end
        
        if order > 1
            % Threshold eigenvalues at a minimum value
            [S,d] = eig(-d2pdT2, 'vector'); % Curvature will be negative. Flip to positive
            eigthresh = opts.MinimumEigenvalueThreshold;
            d(d < eigthresh) = eigthresh;
            d2pdT2 = -S*diag(d)*S.'; % Flip back to negative curvature
        end
    end

%% Test functions

    function [p, dpdT, d2pdT2] = logpdf_sens_test(T)
        p = -0.5*log((2*pi)^nvar*detV) - 0.5*(T-mu).'*Vinv*(T-mu);
        dpdT = -Vinv*(T-mu);
        d2pdT2 = -Vinv;
    end

end