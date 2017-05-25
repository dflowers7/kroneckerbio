function [m,con,obj,opts,localOpts,nT,T0,funopts] = FixFitObjectiveOpts(m, con, obj, opts)
% [m,con,obj,opts,localOpts,nT,T0,funopts] = FixFitObjectiveOpts(m, con, obj, opts)
% Output arguments:
%   funopts
%       Options for building the objective and constraint function, for
%       inputting into GenerateObjective().

% Default options
defaultOpts.Verbose          = 1;

defaultOpts.RelTol           = [];
defaultOpts.AbsTol           = [];

defaultOpts.Normalized       = true;
defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

defaultOpts.ObjWeights       = ones(size(obj));

defaultOpts.UseAdjoint       = true;

defaultOpts.LowerBound       = 0;
defaultOpts.UpperBound       = inf;
defaultOpts.Aeq              = [];
defaultOpts.beq              = [];
defaultOpts.TolOptim         = 1e-5;
defaultOpts.Restart          = 0;
defaultOpts.RestartJump      = 0.001;
defaultOpts.TerminalObj      = -inf;

defaultOpts.MaxStepSize      = 1;
defaultOpts.Solver           = 'fmincon';
if ~isfield(opts, 'Solver') || strcmp(opts.Solver, 'fmincon')
    defaultOpts.Algorithm    = 'active-set';
else
    defaultOpts.Algorithm    = 'trust-region-reflective';
end
defaultOpts.MaxIter          = 1000;
defaultOpts.MaxFunEvals      = 5000;

defaultOpts.OutputFcn              = [];
defaultOpts.ParallelizeExperiments = false;
defaultOpts.TimeoutDuration = [];

defaultOpts.ConstraintObj    = {};
defaultOpts.ConstraintVal    = [];

defaultOpts.IntegrateFunction               = [];
defaultOpts.ObjectiveReductionFunction      = [];
defaultOpts.ConstraintIntegrateFunction     = [];
defaultOpts.ConstraintReductionFunction     = [];
defaultOpts.ScaleConstraints                = false;

defaultOpts.UseImprovedHessianApprox = false;
defaultOpts.HessianApproxMaximumConditionNumber = 1000;
defaultOpts.ApproximateSecondOrderHessianTerm = true;
defaultOpts.HessianApproximation = 'bfgs';
defaultOpts.SubproblemAlgorithm = 'factorization';

defaultOpts.GlobalOptimization = false;
defaultOpts.GlobalOpts         = [];

% Global fit default options
defaultGlobalOpts.Algorithm = 'globalsearch';
defaultGlobalOpts.StartPointsToRun = 'bounds-ineqs';
defaultGlobalOpts.nStartPoints = 10;
defaultGlobalOpts.UseParallel = false;
defaultGlobalOpts.MaxIter = 1000; % for pattern search; fix to better default

% Assign default values and make final options structs
opts = mergestruct(defaultOpts, opts);
opts.GlobalOpts = mergestruct(defaultGlobalOpts, opts.GlobalOpts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Check for global optimization toolbox only if global optimization is specified
%   Note: this isn't necesssary for all global optimization methods, but we depend
%   on this functionality for the current implementation
if opts.GlobalOptimization
    assert(logical(license('test','gads_toolbox')), 'KroneckerBio:FitObjective:GlobalOptimizationToolboxMissing', 'Global optimization requires the global optimization (gads) toolbox.')
end


% Constants
nx = m.nx;
nk = m.nk;
ns = m.ns;

% Ensure structures are proper sizes
[con, n_con] = fixCondition(con);
if ~isempty(obj)
    [obj, n_obj] = fixObjective(obj, n_con);
end

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, n_con);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, n_con, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, n_con, cat(1,con.nh));

nT = nTk + nTs + nTq + nTh;

% Ensure Restart is a positive integer
if ~(opts.Restart >= 0)
    opts.Restart = 0;
    warning('KroneckerBio:FitObjective:NegativeRestart', 'opts.Restart was not nonegative. It has been set to 0.')
end

if ~(opts.Restart == floor(opts.Restart))
    opts.Restart = floor(opts.Restart);
    warning('KroneckerBio:FitObjective:NonintegerRestart', 'opts.Restart was not a whole number. It has been floored.')
end

% Ensure RestartJump is a function handle
if isnumeric(opts.RestartJump)
    opts.RestartJump = @(iter,G)(opts.RestartJump);
end

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
if opts.UseAdjoint
    warning('UseAdjoint is currently not implemented in this version of FitObjective. Switching to false...')
    opts.UseAdjoint = false;
end
if strcmp(opts.Solver, 'lsqnonlin') && opts.UseAdjoint
    warning('KroneckerBio:FitObjective:AdjointNotSupportedForLsqnonlin', ...
        'Adjoint method currently not supported for lsqnonlin. Switching to forward method...')
    opts.UseAdjoint = false;
end
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, n_con, opts.UseAdjoint, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

% Bounds
opts.LowerBound = fixBounds(opts.LowerBound, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
opts.UpperBound = fixBounds(opts.UpperBound, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

obj_constraint = opts.ConstraintObj;
if isstruct(obj_constraint)
    obj_constraint = {obj_constraint};
end
isNonlinearConstraint = ~isempty(obj_constraint);

localOpts = optimoptions(opts.Solver);
localOpts.TolFun                  = opts.TolOptim;
localOpts.TolX                    = 0;
localOpts.MaxFunEvals             = opts.MaxFunEvals;
localOpts.MaxIter                 = opts.MaxIter;

switch opts.Solver
    case 'fmincon'
        localOpts.GradObj                 = 'on';
        localOpts.Hessian                 = 'off'; % unused
        localOpts.RelLineSrchBnd          = opts.MaxStepSize;
        localOpts.RelLineSrchBndDuration  = Inf;
        localOpts.TolCon                  = 1e-6;
        if isNonlinearConstraint
            localOpts.SpecifyConstraintGradient = true;
        end
        localOpts.HessianApproximation    = opts.HessianApproximation; % only used for 'interior-point' algorithm
        localOpts.SubproblemAlgorithm     = opts.SubproblemAlgorithm;
    case 'lsqnonlin'
        localOpts.Jacobian                  = 'on'; % For older versions of MATLAB
        % The following options are for newer versions of MATLAB that
        % change the options' names
        localOpts.FunctionTolerance         = opts.TolOptim; % Similar to TolOptim
        localOpts.StepTolerance             = 0; % Similar to TolX
        localOpts.SpecifyObjectiveGradient  = true; % For newer versions of MATLAB
        localOpts.MaxFunctionEvaluations    = opts.MaxFunEvals;
        localOpts.MaxIterations             = opts.MaxIter;        
        
        assert(~isNonlinearConstraint, 'KroneckerBio:FitObjective:lsqnonlinNonlinearConstraintNotSupported', ...
            'The lsqnonlin solver does not support nonlinear constraints. Use the fmincon solver instead.')
end
localOpts.Algorithm           = opts.Algorithm;
providedOutputFcn             = opts.OutputFcn;

if strcmp(opts.Algorithm, 'sqp')
    % Prevent sqp algorithm from automatically normalizing parameters
    localOpts.ScaleProblem = 'none';
    
    % Set output function that updates what the last parameter set used was
    localOpts.OutputFcn = @outfun;
else
    localOpts.OutputFcn = @isTerminalObj;
end

    function stop = outfun(x, optimValues, state)
        lastT = x;
        stop = isTerminalObj(x, optimValues, state);
    end

if verbose
    localOpts.Display = 'iter';
else
    localOpts.Display = 'off';
end

if opts.Normalized
    % Change relative line search bound to an absolute scale in log space
    % Because fmincon lacks an absolute option, this hack circumvents that
    if strcmp(opts.Solver, 'fmincon') && strcmp(opts.Algorithm, 'active-set')
        % This only works for active-set. Other algorithms don't have the
        % RelLineSrchBnd option, and setting TypicalX won't help on its own.
        localOpts.TypicalX = zeros(nT,1) + log(1 + opts.MaxStepSize)*log(realmax);
        localOpts.RelLineSrchBnd = 1 / log(realmax);
    end
end

% Construct starting variable parameter set
T0 = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

% Apply bounds to starting parameters before optimizing
% fmincon will choose a wierd value if a starting parameter is outside the bounds
belowLowerBounds = T0 < opts.LowerBound;
if any(belowLowerBounds)
    warning('Parameter %g was below its lower bound by %g. Resetting it to its lower bound...', [find(belowLowerBounds).'; (opts.LowerBound(belowLowerBounds)-T0(belowLowerBounds))'])
end
T0(belowLowerBounds) = opts.LowerBound(belowLowerBounds);

aboveUpperBounds = T0 > opts.UpperBound;
if any(aboveUpperBounds)
    warning('Parameter %g is above its upper bound by %g. Resetting it to its upper bound...', [find(aboveUpperBounds).'; (T0(aboveUpperBounds)-opts.UpperBound(aboveUpperBounds)).'])
end
T0(aboveUpperBounds) = opts.UpperBound(aboveUpperBounds);

funopts.ScaleConstraints = opts.ScaleConstraints;
if isempty(opts.HessianApproxMaximumConditionNumber)
    funopts.HessianMaximumConditionNumber = defaultOpts.HessianApproxMaximumConditionNumber;
else
    funopts.HessianMaximumConditionNumber = opts.HessianApproxMaximumConditionNumber;
end
funopts.ApproximateSecondOrderHessianTerm = opts.ApproximateSecondOrderHessianTerm;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Halt optimization on terminal goal %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function stop = isTerminalObj(x, optimValues, state)
        switch opts.Solver
            case 'fmincon'
                fval = optimValues.fval;
            case 'lsqnonlin'
                fval = optimValues.resnorm.^2;
        end
        if strcmp(state,'iter') && fval  < opts.TerminalObj % Only stop at the end of iterations to avoid parameters not being updated
            aborted = true;
            Tabort = x;
            stop = true;
        else
            stop = false;
        end
        
        if ~isempty(opts.OutputFcn)
            stop = stop || providedOutputFcn(x, optimValues, state);
        end
    end

end

