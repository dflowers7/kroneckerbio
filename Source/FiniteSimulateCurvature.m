function sim = FiniteSimulateCurvature(m, con, obs, opts)
%FiniteSimulateSensitivitySelect Approximate the sensitivities of every 
%   species with respect to every parameter over all time returns the
%   values at select time points
%
%   Mathematically: dx/dT = (x(T1) - x(T2)) / (T1 - T2)
%   
%   sim = FiniteSimulateSensitivitySelect(m, con, tGet, opts)
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
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%       .TimeoutDuration [ nonnegative scalar {[]} ]
%           Sets an upper limit to the amount of time an integration may
%           take. Any integration taking longer than this throws an error.
%           If empty (the default), no upper limit is set.
%
%   Outputs
%   sim = FiniteSimulateSensitivitySelect(m, con, tGet, opts)
%   	A vector of structures with each entry being the simulation
%       under one of the conditions.
%       .t tGet
%       .y [ matrix ny by numel(tGet) ]
%           The value of the outputs at each selected time point
%       .x [ matrix nx by numel(tGet) ]
%           The value of the states at each selected time point
%       .dydT [ matrix ny*nT by numel(tGet) ]
%           The value of the sensitivites of the outputs at each selected
%           time point
%       .dxdT [ matrix nx by numel(tGet) ]
%           The value of the sensitivities of the states at each selected
%           time point

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

% Ensure structures are proper sizes
[con, n_con] = fixCondition(con);
[obs, n_obs] = fixObservation(obs, n_con);

derorder = 2;
opts = FixSimulationOpts(m, con, obs, opts, derorder);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

%% Run integration for each experiment
sim = emptystruct(n_con, 'Type', 'Name', 't', 'y', 'x', 'dxdT', 'dudT', 'dydT', 'd2xdT2', 'd2udT2', 'd2ydT2', 'ie', 'te', 'xe', 'ue', 'ye', 'dxedT', 'duedT', 'dyedT', 'd2xedT2', 'd2uedT2', 'd2yedT2', 'int');

for i_con = 1:n_con
    % Modify opts structure
    opts_i = opts;
    opts_i.AbsTol = opts.AbsTol{i_con};
    opts_i.UseSeeds = opts.UseSeeds(:,i_con);
    opts_i.UseInputControls = opts.UseInputControls{i_con};
    opts_i.UseDoseControls = opts.UseDoseControls{i_con};
    
    % Integrate [x; dx/dT] for each finitely perturbed parameter
    if verbose; fprintf(['Integrating curvature for ' con(i_con).Name '...']); end
    ints = integrateAllCurv(m, con(i_con), obs(:,i_con), opts_i, true);
    if verbose; fprintf('done.\n'); end
    
    for i_obs = 1:n_obs
        sim(i_obs,i_con) = pastestruct(sim(i_obs), obs(i_obs).Curvature(ints(i_obs)));
    end
end
