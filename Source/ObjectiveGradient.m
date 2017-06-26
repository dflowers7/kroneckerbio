function D = ObjectiveGradient(m, con, obj, opts)
%ObjectiveGradient Evaluate the gradient of a set of objective functions
%
%   D = ObjectiveGradient(m, con, obj, opts)
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix n_obj by n_con ]
%       The objective structures defining the objective functions to be
%       evaluated. Note that this matrix must have a number of columns
%       equal to numel(con) (e.g. one objective for each experimental
%       condition is a row vector and multiple objective structures for
%       a single experimental conditions is a column vector).
%   opts: [ options struct scalar {} ]
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Which kinetic parameters the gradient will be calculated on
%       .UseSeeds [ logical matrix nx by nCon | logical vector nx |
%                   positive integer vector {[]} ]
%           Which seed parameters the gradient will be calculated on
%       .UseInputControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the input control parameters whose gradient will be
%           calculated
%       .UseDoseControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the dose control parameters whose gradient will be
%           calculated
%     	.ObjWeights [ real matrix nObj by nCon {ones(nObj,nCon)} ]
%           Applies a post evaluation weight on each objective function
%           in terms of how much it will contribute to the final objective
%           function value
%       .Normalized [ logical scalar {true} ]
%           Indicates if the gradient should be computed in log parameters
%           space
%    	.UseAdjoint [ logical scalar {true} ]
%           Indicates whether the gradient should be calculated via the
%           adjoint method or the forward method.
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
%       D: [ real vector nT ]
%           The sum of all objective function gradients

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

assert(nargin >= 3, 'KroneckerBio:ObjectiveGradient:TooFewInputs', 'ObjectiveGradient requires at least 3 input arguments')
assert(isscalar(m), 'KroneckerBio:ObjectiveGradient:MoreThanOneModel', 'The model structure must be scalar')

% Fix options
derorder = 1;
opts = FixObjectiveOpts(m, con, obj, opts, derorder);

%% Run main calculation
[~, D] = computeObjGrad(m, con, obj, opts);
