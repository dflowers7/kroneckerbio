function m = InitializeModelAnalytic(name)
% Create an analytic model. Call FinalizeModel to validate components and
% calculate derivatives for use with the rest of kroneckerbio.
% Inputs:
%   name [ string {''} ]
% Outputs:
%   m [ Model.Analytic struct ]

if nargin < 1
    name = [];
end
if isempty(name)
    name = '';
end

m.Type = 'Model.Analytic';
m.Name = fixModelName(name);

m.Compartments = growCompartments;
m.Parameters   = growParameters;
m.Seeds        = growSeeds;
m.Inputs       = growInputs;
m.States       = growStates;
m.Reactions    = growReactionsAnalytic;
m.Outputs      = growOutputsAnalytic;
m.Rules        = growRules;

m = initializeModelBase(m);

m.add.Compartments = growCompartments;
m.add.Parameters   = growParameters;
m.add.Seeds        = growSeeds;
m.add.Inputs       = growInputs;
m.add.States       = growStates;
m.add.Reactions    = growReactionsAnalytic;
m.add.Outputs      = growOutputsAnalytic;
m.add.Rules        = growRules;
