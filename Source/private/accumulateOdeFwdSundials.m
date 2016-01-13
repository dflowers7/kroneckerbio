function sol = accumulateOdeFwdSundials(der, jac, t0, tF, ic, discontinuities, t_get, RelTol, AbsTol, delta, events, is_finished, nx)
% Outputs
%   sol: struct with the following fields:
%       .t [ 1-by-nt double ]
%           Times at which an output value was desired. This is equal to
%           t_get for successful integrations unless a terminal event
%           occurs before tF, in which case all 
%       .x [ nx-by-nt double ]
%       .y [ ny-by-nt double ]
%       .ie [ 1-by-nevents double ]
%           Index indicating which event occurred. Indices are determined
%           according to the order of events specified in the events input
%           argument.
%       .te [ 1-by-nevents double ]
%           Times at which events occurred.
%       .xe [ nx-by-nevents double ]
%           Values of states at each event time.
%       .ye [ ny-by-nevents double ]
%           Values of outputs at each event time.

% Clean up inputs
if nargin < 12
    is_finished = [];
    if nargin < 11
        events = [];
        if nargin < 10
            delta = [];
        end
    end
end

if isempty(events)
    events = @(t,x) [];
end

if isempty(is_finished)
    is_finished = @(sol) true; % Set to true so that the isterminal value set in the event function determines whether to halt integration
end

% Get number of events being tracked
NumRoots = length(events(0, ic));

% discontinuities will be a column vector sorted ascending. These are
% discontinuities in the output, not in the RHS.
discontinuities = unique([vec(discontinuities); t0; tF]);
N = numel(discontinuities);

% Set integration options
options = CVodeSetOptions(...
    'RelTol', RelTol,...
    'AbsTol', AbsTol,...
    'LinearSolver','Dense',...
    'JacobianFn', jac,...
    'RootsFn', @eventsInterface,...
    'NumRoots', NumRoots...
    );

% Determine stopping times, which are a combination of t_get and the
% discontinuities
t_stop = unique([vec(discontinuities); vec(t_get)]);

% Initialize solution vectors t and x
sol.t = zeros(1,N);
sol.x = zeros(nx,N);

% Initialize event solution vectors te, xe, and ye
eventArraySize = 10; % maximum number of events that can be stored in the initialized array
sol.ie = zeros(1,eventArraySize);
sol.te = zeros(1,eventArraySize);
sol.xe = zeros(nx,eventArraySize);
nevents = 0; % number of events that have occurred

done = false;
ii = 1;
ti = t0;
xi = ic + delta(t0, ic); % add in any delta at the first timepoint

sol.t(1) = t0;
sol.x(:,1) = xi;

% Iterate over the stopping times, integrating to each one from the last
while ~done
    
    ii = ii + 1;
    
    tstopi = t_stop(ii);
    
    % Initialize solver
    options = CVodeSetOptions(options, 'StopTime', tstopi);
    if ii == 2 % First initialization
        CVodeInit(der, 'BDF', 'Newton', ti, xi, options);
    else
        % Reinitialization updates the states with their post-discontinuity
        % values
        CVodeReInit(ti, xi, options);
    end
    
    % Integrate to first discontinuity or event
    [status, ti, xi] = CVode(tstopi, 'Normal');
    
    % Determine why integration stopped and handle accordingly
    switch status
        case -1 % error
            CVodeFree;
            error('KroneckerBio:accumulateOde:IntegrationFailure', 'Did not integrate through entire interval!');
        case 2 % event(s)
            % Determine which events occurred
            stats = CVodeGetStats;
            eventOccurred = stats.RootInfo.roots; % eventOccurred is a vector with 1's at positions of events that triggered
            
            % Check directions (currently not implemented. Does anyone use this?)
            
            % Record each event
            for ie_i = find(eventOccurred(:).' == 1)
                if is_finished(sol) && isTerminal(ti, xi)
                    done = true;
                end
                nevents = nevents + 1;
                % If more events have occurred than can be stored, double
                % array sizes
                if nevents > eventArraySize
                    sol.ie = [sol.ie zeros(1,eventArraySize)];
                    sol.te = [sol.te zeros(1,eventArraySize)];
                    sol.xe = [sol.xe zeros(nx, eventArraySize)];
                    eventArraySize = 2*eventArraySize;
                end
                sol.ie(nevents) = ie_i;
                sol.te(nevents) = ti;
                sol.xe(:,nevents) = xi;
            end
        otherwise % discontinuity and/or requested time
            if ismember(ti, discontinuities) % Discontinuity
                % Apply delta
                if ~isempty(delta)
                    xi = xi + delta(tstopi, xi);
                end
            end
            if ismember(ti, t_get) % Requested time
                % Record time and states
                sol.t(ii) = ti;
                sol.x(:,ii) = xi;
            end
            if ti == tF % Check if reached tF
                done = true;
            end
    end
    
end

% Filter out reserved but unused values
sol.ie = sol.ie(1:nevents);
sol.te = sol.te(1:nevents);
sol.xe = sol.xe(:,1:nevents);

% Free memory
CVodeFree;

    % Function for setting flag in events. Note that isterminal and
    % direction are ignored. These are checked later using the
    % functions below.
    function [G, flag] = eventsInterface(t, x)
        % TODO: This won't work if the event is expecting more than just
        % the states, i.e., in sensitivity or curvature calculations.
        G = events(t, x);
        flag = 0;
    end

    function val = isTerminal(t, x)
        [~,val] = events(t, x);
        val = val == 1;
    end

    % Currently unused. Might be useful in the future.
    function val = directionIsSatisfied(t, x, xprev)
        % requiredDirections: 0 means all directions trigger, 1 means
        % increasing function value triggers, -1 means decreasing function
        % value triggers
        [~,~,requiredDirections] = events(t,x);
        actualDirections = sign(x-xprev);
        val = requiredDirections(:) == 0 | actualDirections(:) == requiredDirections(:);
    end
end