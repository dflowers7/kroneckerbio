function sim = SimulateSensitivity(m, con, obs, opts)
%SimulateSensitivity Integrate the sensitivities of every species with
%   respect to every parameter over all time
%
%   Mathematically: dx/dT = Integral(df/dx * dx/dT + df/dT, t=0:tF)
%   
%   sim = SimulateSensitivity(m, con, opts)
%   
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   opts: [ options struct scalar {} ]
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters whose sensitivities are
%           desired
%       .UseSeeds [ logical matrix ns by nCon | logical vector ns |
%                   positive integer vector {[]} ]
%           Indicates the seed parameters whose sensitivities are desired
%       .UseInputControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the input control parameters whose sensitivites are
%           desired
%       .UseDoseControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the dose control parameters whose sensitivites are
%           desired
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%       .AdjointOutputSensitivities [ logical matrix ny by nCon | logical
%                   vector ny | positive integer vector ([]) ]
%           If .Integrator is set to 'sundials', this option indicates to
%           use the adjoint method to calculate sensitivities for the
%           indicated outputs. If left empty (the default) or if it is a
%           logical matrix of false entries, the forward method will be
%           used to calculate the sensitivities of all states and outputs.
%           If any outputs are indicated, the adjoint method will be used
%           to calculate the sensitivities of only the indicated outputs,
%           bypassing the need to calculate sensitivities for all the
%           states. If .Integrator is anything other than 'sundials', this
%           option is ignored, and the forward method is used.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%       .TimeoutDuration [ nonnegative scalar {[]} ]
%           Sets an upper limit to the amount of time an integration may
%           take. Any integration taking longer than this throws an error.
%           If empty (the default), no upper limit is set.
%
%   Outputs
%   SimulateSensitivity(m, con, opts)
%   	Plots the sensitivities under each condition
%
%   sim = SimulateSensitivity(m, con, opts)
%   	A vector of structures with each entry being the simulation
%       under one of the conditions.
%       .t [ sorted nonnegative row vector ]
%           Timepoints chosen by the ode solver
%       .y [ handle @(t,y) returns matrix numel(y) by numel(t) ]
%           This function handle evaluates some outputs y of the system at
%           some particular time points t. The user may exclude y, in which
%           case all outputs are returned.
%       .x [ handle @(t,x) returns matrix numel(x) by numel(t) ]
%           This function handle evaluates some states x of the system at
%           some particular time points t. The user may exclude x, in which
%           case all states are returned.
%       .dydT [ handle @(t,y) returns matrix numel(y)*nT by numel(t) ]
%           This function handle evaluates the sensitivity of some outputs
%           y to the active parameters of the system at some particular
%           time points t. The user may exclude y, in which case all
%           outputs are returned.
%       .dxdT [ handle @(t,x) returns matrix numel(x) by numel(t) ]
%           This function handle evaluates the sensitivity of some states
%           x to the active parameters of the system at some particular
%           time points t. The user may exclude x, in which case all
%           states are returned.
%       .sol [ odesolver struct scalar ]
%           The integrator solution to the system

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

assert(nargin >= 3, 'KroneckerBio:SimulateSensitivity:TooFewInputs', 'SimulateSensitivity requires at least 2 input arguments')
assert(isscalar(m), 'KroneckerBio:SimulateSensitivity:MoreThanOneModel', 'The model structure must be scalar')

% Ensure structures are proper sizes
[con, n_con] = fixCondition(con);
[obs, n_obs] = fixObservation(obs, n_con);

% Default options
derorder = 1;
[opts,nTstruct] = FixSimulationOpts(m, con, obs, opts, derorder);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Determine which parameters are fit by which experiments
T_experiment = zeros(nTstruct.nT, 1);
T_experiment(1:nTstruct.nTk) = 0;
nTs_con = sum(opts.UseSeeds,1);
nTq_con = cellfun(@sum, opts.UseInputControls(:)');
nTh_con = cellfun(@sum, opts.UseDoseControls(:)');
nT_con = {nTs_con, nTq_con, nTh_con};
endi = nTstruct.nTk;
for i_type = 1:3
    for j_con = 1:n_con
        starti = endi + 1;
        endi = endi + nT_con{i_type}(j_con);
        T_experiment(starti:endi) = j_con;
    end
end
% Determine which T's are fit by which experiments
TisExperiment = bsxfun(@(T_experiment,i_con)T_experiment == 0 | T_experiment == i_con, T_experiment, 1:n_con);

%% Run integration for each experiment
sim = emptystruct([n_obs,n_con]);

for i_con = 1:n_con
    % Modify opts structure
    opts_i = opts;
    opts_i.AbsTol = opts.AbsTol{i_con};
    opts_i.UseSeeds = opts.UseSeeds(:,i_con);
    opts_i.UseInputControls = opts.UseInputControls{i_con};
    opts_i.UseDoseControls = opts.UseDoseControls{i_con};
    opts_i.AdjointOutputSensitivities = opts.AdjointOutputSensitivities(:,i_con);
    
    % Integrate [x; dx/dT] over time
    if verbose; fprintf(['Integrating sensitivities for ' con(i_con).Name '...']); end
    ints = integrateAllSens(m, con(i_con), obs(:,i_con), opts_i);
    [ints.TisExperiment] = deal(TisExperiment(:,i_con));
    if verbose; fprintf('done.\n'); end
    
    for i_obs = 1:n_obs
        sim = insertstruct(sim, obs(i_obs,i_con).Sensitivity(ints(i_obs)), i_obs,i_con);
    end
end
