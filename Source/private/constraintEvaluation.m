function varargout = constraintEvaluation(m,con,opts)

%% Clean up inputs

assert(nargin >= 3, 'KroneckerBio:ConstraintEvaluation:TooFewInputs', 'ConstraintEvaluation requires at least 3 input arguments')

assert(isscalar(m), 'KroneckerBio:ConstraintEvaluation:MoreThanOneModel', 'The model structure must be scalar')

% Create a dummy obj
n_con = numel(con);
obj = objectiveZero([1 n_con]);

[m,con,~,opts,localOpts,nT,T0,funopts] = FixFitObjectiveOpts(m, con, obj, opts);

order = nargout-1;

%% Generate constraint function

[~,constraint] = GenerateObjective(m, con, obj, opts, ...
    opts.IntegrateFunction, opts.ObjectiveReductionFunction, ...
    opts.ConstraintIntegrateFunction, opts.ConstraintReductionFunction, ...
    funopts);

if order == 0
    varargout = cell(2,1);
elseif order == 1
    varargout = cell(4,1);
else
    error('KroneckerBio:ConstraintEvaluation:UnsupportedOrder', ...
        'constraintEvaluation currently only supports calculating derivatives of order 0 or 1.')
end

[varargout{:}] = constraint(T0);

if order == 0
    varargout = varargout(1);
elseif order == 1
    varargout = varargout([1 3]);
end

end