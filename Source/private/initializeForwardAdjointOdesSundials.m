function options = initializeForwardAdjointOdesSundials(der, jac, g, t0, x0, g0, RelTol, AbsTol, delta, events)
% Initializes SUNDIALS to integrate the state system. To calculate the
% state values, call accumulateOdeFwdSundials following this function.

% Initialize state system
options = initializeStateOdesSundials(der, jac, t0, x0, RelTol, AbsTol, delta, events);

% Initialize quadrature system, if supplied
if ~isempty(g)
    CVodeQuadInit(g, g0);
end

% Initialize adjoint system
maxStepsBetweenCheckpoints = 1500; % I don't know what this value should be in general, so I am using the value from the example code times ten
CVodeAdjInit(maxStepsBetweenCheckpoints, 'Polynomial');