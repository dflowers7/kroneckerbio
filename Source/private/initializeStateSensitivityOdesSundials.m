function options = initializeStateSensitivityOdesSundials(der, t0, x0, dx0dT, RelTol, AbsTol, delta)
% Initializes SUNDIALS to integrate the state sensitivity system using the
% forward method. To calculate the state sensitivity values, call
% accumulateOdeFwdSundials following this function.

% Clean up inputs
if nargin < 7
    delta = [];
end

nx = numel(x0);
nT = numel(dx0dT)/nx;

% Set integration options
options = CVodeSensSetOptions(...
    'RelTol', RelTol,...
    'AbsTol', reshape(AbsTol, nx, nT)...
    );

if ~isempty(delta)
    dx0dT = dx0dT + delta(t0, x0); % add in the delta at the first timepoint
end

% Initialize solver
CVodeSensInit(nT, der, dx0dT, options); 

end