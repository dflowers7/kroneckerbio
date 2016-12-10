function [m, con, G, D] = FitObjective(m, con, obj, opts)
%FitObjective Optimize the parameters of a model to minimize a set of
%   objective functions
%
%   Mathematically: T = argmin(G(T))
%
%   [m, con, G, D] = FitObjective(m, con, obj, opts)
%
%   FitObjective uses the derivatives in the Kronecker model and in the
%   objective functions to build a function that can not only evaluate the
%   objective functions at particular parameter sets, but also evaluate the
%   gradient at those parameter sets, thereby pointing in the direction of
%   a more optimum parameter set. This function is built around Matlab's
%   fmincon, which is a gradient descent minimizer. It varies the
%   parameters attempting to find the parameter set that will minimize the
%   objective function while keeping with the bounds.
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector n_con ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix n_obj by n_con ]
%       The objective structures defining the objective functions to be
%       evaluated. Note that this matrix must have a number of columns
%       equal to numel(con) (e.g. one objective for each experimental
%       condition is a row vector and multiple objective structures for
%       a single experimental conditions is a column vector).
%   opts: [ options struct scalar {} ]
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters that will be allowed to vary
%           during the optimization
%       .UseSeeds [ logical matrix ns by nCon | logical vector ns |
%                   positive integer vector {[]} ]
%           Indicates the seeds that will be allowed to vary during the
%           optimzation. If UseModelSeeds is true then UseSeeds can be a
%           vector of linear indexes or a vector of logicals length of ns.
%           If UseModelSeeds is false then UseSeeds can be a matrix of
%           logicals size ns by nCon. It can also be a vector of length ns,
%           and every experiment will be considered to have the same active
%           seed parameters. It can also be a vector of linear indexes into
%           the ns vector and assumed the same for all conditions.
%       .UseInputControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization
%       .UseDoseControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization
%       .LowerBound [ nonegative vector {0} ]
%           The lower bound on the fitted parameters. It can be length
%           nk+n_con*nx, nk+nx, nT, just nk if nTx = 0, or a scalar. The
%           bounds will be interpreted in that order if the length matches
%           multiple orders.
%       .UpperBound [ nonegative vector {0} ]
%           The upper bound for the fitted parameters. It must be the same
%           length as LowerBound.
%     	.ObjWeights [ real matrix nObj by nCon {ones(nObj,nCon)} ]
%           Applies a post evaluation weight on each objective function
%           in terms of how much it will contribute to the final objective
%           function value.
%       .Normalized [ logical scalar {true} ]
%           Indicates if the optimization should be done in log parameters
%           space
%    	.UseAdjoint [ logical scalar {false} ]
%           Indicates whether the gradient should be calculated via the
%           adjoint method or the forward method
%     	.TolOptim [ positive scalar {1e-5} ]
%           The objective tolerance. The optimization stops when it is
%           predicted that the objective function cannot be improved more
%           than this in the next iteration.
%     	.Restart [ nonnegative integer scalar {0} ]
%           A scalar integer determining how many times the optimzation
%           should restart once optimization has stopped.
%     	.RestartJump [ handle @(iter,G) returns nonnegative vector nT or
%                      scalar | nonnegative vector nT or scalar {0.001} ]
%           This function handle controls the schedule for the noise that
%           will be added to the parameters before each restart. The
%           parameters for the next iteration will be normally distributed
%           in log space with a mean equal to the previous iteration and a
%           standard deviation equal to the value returned by this
%           function. The value returned should usually be a scalar, but it
%           can also be a vector with length equal to the number of active
%           parameters. It can also be numeric, and the noise will be
%           treated as this constant value.
%      	.TerminalObj [ real scalar {-inf} ]
%           Optimization is halted when this objective function value is
%           reached
%       .MaxStepSize [ nonegative scalar {1} ]
%           Scalar fraction indicator of the maximum relative step size
%           that any parameter can take in a single interation
%     	.Algorithm [ string {active-set} ]
%           Option for fmincon. Which optimization algorithm to use
%     	.MaxIter [ postive scalar integer {1000} ]
%           Option for fmincon. Maximum number of iterations allowed before
%           optimization will be terminated.
%     	.MaxFunEvals [ postive scalar integer {5000} ]
%           Option for fmincon. Maximum number of objective function
%           evaluations allowed before optimization will be terminated.
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%       .ParallelizeExperiments [ logical scalar {false} ]
%           Set to true to parallelize local optimization over the
%           experiments. The optimization will use the current parallel
%           pool, if one is available, or initialize a parallel pool in the
%           default cluster, if Parallel Toolbox preferences are set to
%           allow initialization of parallel pools on calling of parallel
%           keywords. If no pool is available and cannot be initialized,
%           the optimization will run serially. This option has no effect
%           on global optimization; set opts.GlobalOpts.UseParallel to true
%           to parallelize global optimization.
%       .OutputFcn [ function handle {[]} ]
%           Set to a function handle of the form stop = outfun(x,
%           optimValues, state). The function is run after each iteration
%           of the optimization. x are the parameters after the iteration,
%           optimValues is a struct containing the current iteration's
%           data, and state is a string describing what the algorithm is
%           currently doing. stop is a flag that, when set to true, stops
%           optimization. See the documentation for optimization options
%           for more details.
%       .ConstraintObj [ cell vector of objective struct matrices n_obj_constraints by n_con]
%           If set, the provided objective structs are used as a nonlinear
%           constraint on the fit. The fit will be constrained to keeping
%           the objective value less than or equal to opts.ConstraintVal.
%           This option is only supported for the fmincon solver.
%       .ConstraintVal [ vector of doubles ]
%           Set each element of the vector to the upper bound on the
%           corresponding constraint objective function's value in
%           opts.ConstraintObj.
%       .GlobalOptimization [ logical scalar {false} ]
%           Use global optimization in addition to fmincon
%       .GlobalOpts [ options struct scalar {} ]
%           TODO: API in progress
%
%   Outputs
%       m: [ model scalar ]
%           The model with all the optimum kinetic parameters applied, as
%           well as the IC parameters if UseModelSeeds = true.
%       con: [ experiment vector ]
%           The experimental conditions which will have the optimum IC
%           parameters applied if UseModelSeeds = false.
%       G: [ real scalar ]
%           The optimum objective function value
%       D: [ real vector nT ]
%           The objective gradient at the optimum parameter set

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Check for optimization toolbox
% Note: this isn't necessary for simple local optimzation methods like
%   fminsearch but we don't use that
assert(logical(license('test','optimization_toolbox')), 'KroneckerBio:FitObjective:OptimizationToolboxMissing', 'FitObjective requires the optimization toolbox for fmincon.')

% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:FitObjective:TooFewInputs', 'FitObjective requires at least 3 input arguments')
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:FitObjective:MoreThanOneModel', 'The model structure must be scalar')

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

defaultOpts.ConstraintObj    = [];
defaultOpts.ConstraintVal    = [];

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
[obj, n_obj] = fixObjective(obj, n_con);

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

% Construct starting variable parameter set
T0 = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

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

%% Options structure for integration
%intOpts = opts;

%% Local optimization options

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
        localOpts.SpecifyConstraintGradient = isNonlinearConstraint;
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

%% Global optimization options
% TODO: make sure options are relevant for solver
globalOpts = opts.GlobalOpts;

% Sanity checking
if strcmpi(globalOpts.Algorithm, 'multistart') && globalOpts.UseParallel
    warning('KroneckerBio:FitObjective:InvalidMultistartOpts', 'Using multistart with UseParallel is not supported at this time (due to global variable in obj fun usage).')
end

%% Normalize parameters

if opts.Normalized
    % Normalize starting parameters and bounds
    T0 = log(T0);
    opts.LowerBound = log(opts.LowerBound);
    opts.UpperBound = log(opts.UpperBound);
    
    % Change relative line search bound to an absolute scale in log space
    % Because fmincon lacks an absolute option, this hack circumvents that
    if strcmp(opts.Solver, 'fmincon') && strcmp(opts.Algorithm, 'active-set')
        % This only works for active-set. Other algorithms don't have the
        % RelLineSrchBnd option, and setting TypicalX won't help on its own.
        localOpts.TypicalX = zeros(nT,1) + log(1 + opts.MaxStepSize)*log(realmax);
        localOpts.RelLineSrchBnd = 1 / log(realmax);
    end
end

%% Apply bounds to starting parameters before optimizing
% fmincon will choose a wierd value if a starting parameter is outside the bounds
T0(T0 < opts.LowerBound) = opts.LowerBound(T0 < opts.LowerBound);
T0(T0 > opts.UpperBound) = opts.UpperBound(T0 > opts.UpperBound);

%% Abort in rare case of no optimization
if numel(T0) == 0
    [G, D] = objective(T0);
    return
end

%% Run optimization
% Initialize loop
That = T0;
Gbest = inf;
Tbest = T0;

% Disable experiment parallelization if global optimization is desired
if opts.GlobalOptimization && opts.ParallelizeExperiments
    warning('KroneckerBio:FitObjective:GlobalOptimizationParallelExperimentsNotSupported', ...
        ['Both opts.GlobalOptimization and opts.ParallelizeExperiments ' ...
        'were set to true. Disabling parallelized experiments.'])
    opts.ParallelizeExperiments = false;
end

% Generate the objective and constraint functions
% TODO: set up adjoint integration function
% Fix global optimization (currently untested)
[objective,constraint] = GenerateObjective(m, con, obj, opts);

for iRestart = 1:opts.Restart+1
    % Init abort parameters
    aborted = false;
    Tabort = That;
    lastT = That;
    
    % Run specified optimization
    if opts.GlobalOptimization
        
        % Create local optimization problem
        %   Used as a subset/refinement of global optimization
        localProblem = createOptimProblem('fmincon', 'objective', objective, ...
            'x0', That, 'Aeq', opts.Aeq, 'beq', opts.beq, ...
            'lb', opts.LowerBound, 'ub', opts.UpperBound, ...
            'options', localOpts);
        
        if opts.Verbose
            fprintf('Beginning global optimization with %s...\n', globalOpts.Algorithm)
        end
        
        switch globalOpts.Algorithm
            case 'globalsearch'
                gs = GlobalSearch('StartPointsToRun', globalOpts.StartPointsToRun);
                [That, G, exitflag] = run(gs, localProblem);
            case 'multistart'
                ms = MultiStart('StartPointsToRun', globalOpts.StartPointsToRun, ...
                    'UseParallel', globalOpts.UseParallel);
                [That, G, exitflag] = run(ms, localProblem, globalOpts.nStartPoints);
            case 'patternsearch'
                psOpts = psoptimset('MaxIter', globalOpts.MaxIter, 'UseParallel', globalOpts.UseParallel);
                [That, G, exitflag] = patternsearch(objective, That, [], [], ...
                    opts.Aeq, opts.beq, opts.LowerBound, opts.UpperBound, [], psOpts);
            otherwise
                error('KroneckerBio:FitObjective:InvalidGlobalOptAlgorithm', 'Global optimization algorithm %s not recognized.', globalOpts.Algorithm)
        end
        
        [~, D] = objective(That); % since global solvers don't return gradient at endpoint
        
    else
        if opts.Verbose; fprintf('Beginning gradient descent...\n'); end
        switch opts.Solver
            case 'fmincon'
                [That, G, exitflag, ~, ~, D] = fmincon(objective, That, [], [], opts.Aeq, opts.beq, opts.LowerBound, opts.UpperBound, constraint, localOpts);
            case 'lsqnonlin'
                [That,G,~,exitflag,~,~,D] = lsqnonlin(objective, That, opts.LowerBound, opts.UpperBound, localOpts);
            otherwise
                error('KroneckerBio:FitObjective:UnrecognizedSolver', ...
                    'Unrecognized solver %s', opts.Solver);
        end
    end
    
    % Check abortion status
    % Abortion values are not returned by fmincon and must be retrieved
    if aborted
        That = Tabort;
        [G, D] = objective(That);
    end
    
    % Re-apply stiff bounds
    That(That < opts.LowerBound) = opts.LowerBound(That < opts.LowerBound);
    That(That > opts.UpperBound) = opts.UpperBound(That > opts.UpperBound);
    
    % See if these parameters are better than previous iterations
    if G < Gbest
        Gbest = G;
        Tbest = That;
    end

    % Terminate if goal has been met
    if G <= opts.TerminalObj
        break
    end
    
    % Jump parameters before restarting
    % Retain the deterministic nature of fitting by fixing the stream
    if iRestart < opts.Restart + 1
        rng_state = rng;
        prodThat = prod(That); % model fingerprint
        rng(mod(prodThat/eps(prodThat)*iRestart,2^32));
        if opts.Normalized
            That = Tbest + randn(nT,1) .* vec(opts.RestartJump(iRestart,G));
        else
            That = exp(log(Tbest) + randn(nT,1) .* vec(opts.RestartJump(iRestart,G)));
        end
        rng(rng_state);
        
        % Prevent jumps from leaving bounds
        while any(That < opts.LowerBound) || any(That > opts.UpperBound)
            That(That < opts.LowerBound) = 2*(opts.LowerBound(That < opts.LowerBound)) - That(That < opts.LowerBound);
            That(That > opts.UpperBound) = 2*(opts.UpperBound(That > opts.UpperBound)) - That(That > opts.UpperBound);
        end
    end
end

% Unnormalize
if opts.Normalized
    Tbest = exp(Tbest);
end

% Update parameter sets
[m, con] = updateAll(m, con, Tbest, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        if fval  < opts.TerminalObj
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
