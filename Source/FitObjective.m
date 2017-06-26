function [m, con, G, D, exitflag, output] = FitObjective(m, con, obj, opts)
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
%    	.UseAdjoint [ logical scalar {true} ]
%           Indicates whether the gradient should be calculated via the
%           adjoint method or the forward method
%       .TolX [ positive scalar {0} ]
%           Tolerance on change in parameter values. If the step size is less
%           than this value, optimization stops.
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
%       .IntegrateFunction
%       .ObjectiveReductionFunction [ function handle ]
%           The function, fun(obj,int), has two input arguments: (1) obj,
%           the objective struct, and (2) int, the integration struct
%           returned by integrating obj over the model and experiments. It
%           returns up to two output arguments: (1) the overall objective
%           value to be minimized, and (2) the gradient of the objective
%           value with respect to the fit parameters, expressed as a column
%           vector. The default function sums obj(i).G(int(i)) over i to
%           calculate the objective value.
%       .ConstraintIntegrateFunction
%       .ConstraintReductionFunction
%       .UseImprovedHessianApprox [ false ]
%           If set to true, and opts.Algorith is set to 'interior-point',
%           an improved approximation of the Hessian is used instead of the
%           default bfgs update. The improved Hessian uses the Fisher
%           Information Matrix in addition to other optional terms to
%           approximate the Hessian for least-squares objective functions.
%       .HessianApproxMaximumConditionNumber [ 1000 ]
%           If left empty, the default value will be used. Is only used if
%           .UseImprovedHessianApprox is set to true.
%       .ApproximateSecondOrderHessianTerm [ true ]
%           If set to true, a structured BFGS update is used to approximate
%           the Hessian term corresponding to second derivatives of the
%           error function of least squares objective functions. This is
%           only used if .UseImprovedHessianApprox is set to true.
%       .HessianGuess
%           'FIM' or 'I' (the default). Only has an effect when the
%           .Solver option is set to 'sqp'.
%       .SubproblemAlgorithm
%           'factorization' or 'cg'
%       .TimeoutDuration [ nonnegative scalar {[]} ]
%           Sets an upper limit to the amount of time an integration may
%           take. Any integration taking longer than this throws an error.
%           If empty (the default), no upper limit is set.
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
%       exitflag [ numeric scalar ]
%           Indicates why the fit terminated. Values' meanings differ
%           depending on the fit function used. See documentation for
%           fmincon and other fitting functions for details.

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

[m,con,obj,opts,localOpts,nT,T0,funopts] = FixFitObjectiveOpts(m, con, obj, opts);

% Normalize parameters and bounds
if opts.Normalized
    opts.LowerBound = log(opts.LowerBound);
    opts.UpperBound = log(opts.UpperBound);
    T0 = log(T0);
end

%% Global optimization options
% TODO: make sure options are relevant for solver
globalOpts = opts.GlobalOpts;

% Sanity checking
if strcmpi(globalOpts.Algorithm, 'multistart') && globalOpts.UseParallel
    warning('KroneckerBio:FitObjective:InvalidMultistartOpts', 'Using multistart with UseParallel is not supported at this time (due to global variable in obj fun usage).')
end

%% Generate objective and constraint functions

% Generate the objective and constraint functions
% TODO: set up adjoint integration function
% Fix global optimization (currently untested)
[objective,constraint,hessian,outfun] = GenerateObjective(m, con, obj, opts, ...
    opts.IntegrateFunction, opts.ObjectiveReductionFunction, ...
    opts.ConstraintIntegrateFunction, opts.ConstraintReductionFunction, ...
    funopts, localOpts.OutputFcn);

switch opts.HessianGuess
    case 'I'
        localOpts = struct(localOpts);
        localOpts.HessMatrix = eye(nT);
    case 'FIM'
        localOpts = struct(localOpts);
        if isempty(constraint)
            temp1 = [];
        else
            [temp1,~,temp2,~] = constraint(T0);
        end
        FIM = hessian(T0, zeros(numel(temp1),1));
        localOpts.HessMatrix = FIM;
end

if opts.UseImprovedHessianApprox
    switch opts.Algorithm
        case 'sqp'
            localOpts.HessFun = hessian;
        otherwise
            localOpts.HessianFcn = hessian;
    end
    localOpts.OutputFcn = outfun;
end

%% Abort in rare case of no optimization
if numel(T0) == 0
    [G, D] = fminconObjective(T0);
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

output = struct;

for iRestart = 1:opts.Restart+1
    % Init abort parameters
    aborted = false;
    Tabort = That;
    lastT = That;
    
    % Create local optimization problem
    localProblem = createOptimProblem('fmincon', 'objective', objective, ...
        'x0', That, 'Aeq', opts.Aeq, 'beq', opts.beq, ...
        'lb', opts.LowerBound, 'ub', opts.UpperBound, ...
        'nonlcon', constraint, 'options', localOpts);
    
    % Run specified optimization
    if opts.GlobalOptimization
        
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
        
        [~, D] = fminconObjective(That); % since global solvers don't return gradient at endpoint
        
    else
        if opts.Verbose; fprintf('Beginning gradient descent...\n'); end
        switch opts.Solver
            case 'fmincon'
                %[That, G, exitflag, output, ~, D] = fmincon(objective, That, [], [], opts.Aeq, opts.beq, opts.LowerBound, opts.UpperBound, constraint, localOpts);
                [That, G, exitflag, output, ~, D] = fmincon(localProblem);
            case 'lsqnonlin'
                [That,G,~,exitflag, output, ~, D] = lsqnonlin(objective, That, opts.LowerBound, opts.UpperBound, localOpts);
            case 'sqp'
                localProblem.options = localOpts;
                [That, output, ~, ~, exitflag] = sqp(localProblem);
                % output(8)  = value of the function at the solution
                % output(10) = number of function evaluations
                % output(11) = number of gradient evaluations
                % output(15) = number of iterations
                [G,D] = objective(That);
            otherwise
                error('KroneckerBio:FitObjective:UnrecognizedSolver', ...
                    'Unrecognized solver %s', opts.Solver);
        end
    end
    
    % Check abortion status
    % Abortion values are not returned by fmincon and must be retrieved
    if aborted
        That = Tabort;
        [G, D] = fminconObjective(That);
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



end
