function [m,con,obj,opts,localOpts,nT,T0,funopts] = FixFitObjectiveOpts(m, con, obj, opts)
% [m,con,obj,opts,localOpts,nT,T0,funopts] = FixFitObjectiveOpts(m, con, obj, opts)
% Output arguments:
%   funopts
%       Options for building the objective and constraint function, for
%       inputting into GenerateObjective().

% Default options
defaultOpts.UseDoseControls  = [];

defaultOpts.ObjWeights       = ones(size(obj));

defaultOpts.LowerBound       = 0;
defaultOpts.UpperBound       = inf;
defaultOpts.Aeq              = [];
defaultOpts.beq              = [];
defaultOpts.TolOptim         = 1e-5;
defaultOpts.TolX             = 0;
defaultOpts.ConstraintTolerance = 1e-6;
defaultOpts.Restart          = 0;
defaultOpts.RestartJump      = 0.001;
defaultOpts.TerminalObj      = -inf;

defaultOpts.MaxStepSize      = 1;
defaultOpts.Solver           = 'fmincon';
if ~isfield(opts, 'Solver') || strcmp(opts.Solver, 'fmincon') || strcmp(opts.Solver, 'sqp')
    defaultOpts.Algorithm    = 'active-set';
else
    defaultOpts.Algorithm    = 'trust-region-reflective';
end
defaultOpts.MaxIter          = 1000;
defaultOpts.MaxFunEvals      = 5000;

defaultOpts.OutputFcn              = [];

defaultOpts.ConstraintObj    = {};
defaultOpts.ConstraintVal    = [];

defaultOpts.IntegrateFunction               = [];
defaultOpts.ObjectiveReductionFunction      = [];
defaultOpts.ConstraintIntegrateFunction     = [];
defaultOpts.ConstraintReductionFunction     = [];
defaultOpts.ScaleConstraints                = false;

defaultOpts.UseImprovedHessianApprox = true;
defaultOpts.HessianApproxMaximumConditionNumber = Inf;
defaultOpts.ApproximateSecondOrderHessianTerm = false;
defaultOpts.HessianApproximation = 'bfgs';
defaultOpts.SubproblemAlgorithm = 'factorization';
defaultOpts.HessianGuess = 'I';

defaultOpts.GoodObjectiveFunctionValue = [];
defaultOpts.GoodConstraintFunctionValues = [];

defaultOpts.UseAdjoint = false;

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

% Fix simulation options after doing the above to account for any changes
% to UseAdjoint and its effects on AbsTol
derorder = 1;
opts = FixSimulationOpts(m, con, obj, opts, derorder);

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
ny = m.ny;

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

% Bounds
opts.LowerBound = fixBounds(opts.LowerBound, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
opts.UpperBound = fixBounds(opts.UpperBound, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

obj_constraint = opts.ConstraintObj;
if isstruct(obj_constraint)
    obj_constraint = {obj_constraint};
end
isNonlinearConstraint = ~isempty(obj_constraint);

switch opts.Solver
    case 'sqp'
        % Use optimset because I need a struct, not an object, when using
        % the sqp solver
        localOpts = optimset('fmincon');
    otherwise
        localOpts = optimoptions(opts.Solver);
end
localOpts.TolFun                  = opts.TolOptim;
localOpts.TolX                    = opts.TolX;
localOpts.MaxFunEvals             = opts.MaxFunEvals;
localOpts.MaxIter                 = opts.MaxIter;
localOpts.ConstraintTolerance     = opts.ConstraintTolerance;

switch opts.Solver
    case 'fmincon'
        localOpts.GradObj                 = 'on';
        localOpts.Hessian                 = 'off'; % unused
        localOpts.RelLineSrchBnd          = opts.MaxStepSize;
        localOpts.RelLineSrchBndDuration  = Inf;
        localOpts.TolCon                  = opts.ConstraintTolerance;
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
    case 'sqp'
        localOpts.GradObj     = 'on';
        localOpts.GradConstr  = 'on';
        localOpts.MaxStepSize = opts.MaxStepSize;
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

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1, 0);

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

% Determine the number of data points to estimate a "good" objective function
% value and constraint values, if one wasn't provided
if isempty(opts.GoodObjectiveFunctionValue)
    errs = ObjectiveErrors(m, con, obj);
    nDataPoints = numel(errs);
    opts.GoodObjectiveFunctionValue = chi2inv(0.95, nDataPoints);
    if opts.GoodObjectiveFunctionValue == 0
        warning('No objective function value was provided. Assuming objective function value should not factor into whether the objective function is allowed to increase.')
        opts.GoodObjectiveFunctionValue = Inf;
    end
end
if isempty(opts.GoodConstraintFunctionValues)
    errs_cons = cell(size(opts.ConstraintObj));
    for i = 1:numel(errs_cons)
        errs_cons{i} = ObjectiveErrors(m, con, opts.ConstraintObj{i});
    end
    nDataPoints = cellfun(@numel, errs_cons);
    opts.GoodConstraintFunctionValues = chi2inv(0.95, nDataPoints);
    noLeastSquares = opts.GoodConstraintFunctionValues == 0;
    if any(noLeastSquares)
        warning('Some constraint functions were not provided examples of good values. Assuming constraint values should not factor into whether the objective function should be allowed to increase.')
        opts.GoodConstraintFunctionValues(noLeastSquares) = Inf;
    end
end

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

