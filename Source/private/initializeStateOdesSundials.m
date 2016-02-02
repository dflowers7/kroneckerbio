function options = initializeStateOdesSundials(der, jac, t0, x0, RelTol, AbsTol, delta, events)
% Initializes SUNDIALS to integrate the state system. To calculate the
% state values, call accumulateOdeFwdSundials following this function.

% Clean up inputs
if nargin < 8
    events = [];
    if nargin < 7
        delta = [];
    end
end

if isempty(events)
    events = @(t,x) [];
end

% Get number of events being tracked
NumRoots = length(events(t0, x0));

% Set integration options
options = CVodeSetOptions(...
    'RelTol', RelTol,...
    'AbsTol', AbsTol,...
    'LinearSolver','Dense',...
    'JacobianFn', jac,...
    'RootsFn', @eventsInterface,...
    'NumRoots', NumRoots...
    );
%     'MonitorFn', 'CVodeMonitor',...
%     'MonitorData', struct('sol', true, 'stats', false)...
%    );

if ~isempty(delta)
    x0 = x0 + delta(t0, x0); % add in the delta at the first timepoint
end

% Initialize solver
CVodeInit(der, 'BDF', 'Newton', t0, x0, options); 

    % Function for setting flag in events. Note that isterminal and
    % direction are ignored. These are checked later using the
    % functions below.
    function [G, flag] = eventsInterface(t, x)
        % TODO: This won't work if the event is expecting more than just
        % the states, i.e., in sensitivity or curvature calculations.
        G = events(t, x);
        flag = 0;
    end

end