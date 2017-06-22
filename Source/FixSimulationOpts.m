function opts = FixSimulationOpts(m, con, obs, opts, derorder)
% opts = FixSimulationOpts(m, con, obs, opts, derorder)
% Input arguments:
%   derorder
%       0 for system, 1 for sensitivities, 2 for curvature

%% Constants

% Constants
nx = m.nx;
nk = m.nk;
ns = m.ns;
ny = m.ny;

% Standardize con
[con, n_con] = fixCondition(con);

%% Set default option fields

defaultOpts.Verbose          = 1;

defaultOpts.RelTol           = [];
defaultOpts.AbsTol           = [];
defaultOpts.AbsTolY          = [];
defaultOpts.RelTolY          = [];

defaultOpts.Normalized       = true;
defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

defaultOpts.UseAdjoint       = false;

% May need a more detailed treatment here later
defaultOpts.AdjointOutputSensitivities = false(m.ny,n_con);

defaultOpts.ParallelizeExperiments = false;
defaultOpts.TimeoutDuration = [];

defaultOpts.Integrator = '';

opts = mergestruct(defaultOpts, opts);

%% Process options

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, n_con);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, n_con, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, n_con, cat(1,con.nh));

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obs);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);
opts.RelTolY = fixRelTolY(opts.RelTolY);

opts.AbsTol = fixAbsTol(opts.AbsTol, derorder+1, opts.continuous, nx, n_con, opts.UseAdjoint, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
opts.AbsTolY = fixAbsTolY(opts.AbsTolY, ny, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls, derorder);

% Fix adjoint output sensitivity specification
opts.AdjointOutputSensitivities = fixAdjointOutputSensitivities(opts.AdjointOutputSensitivities, ny, n_con);

end