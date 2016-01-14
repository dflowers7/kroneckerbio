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
N = numel(t_get);

% Set integration options
options = CVodeSetOptions(...
    'RelTol', RelTol,...
    'AbsTol', AbsTol,...
    'LinearSolver','Dense',...
    'JacobianFn', jac,...
    'RootsFn', @eventsInterface,...
    'NumRoots', NumRoots...
    );

% Determine stopping times, which are the discontinuities
%t_stop = unique([vec(discontinuities); vec(t_get)]);
t_stop = discontinuities;

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
tstop_i = t0;
x_stop = ic + delta(t0, ic); % add in any delta at the first timepoint
isDiscontinuity = true;

sol.t(1) = t0;
sol.x(:,1) = x_stop;

% Iterate over the discontinuities, integrating to each one from the last
while ~done
    
    ii = ii + 1;
    
    tstart_i = tstop_i;
    tstop_i = t_stop(find(t_stop > tstart_i, 1, 'first'));
    
    % Initialize options and solver
    if ii == 2 % First initialization
        %options = CVodeSetOptions(options, 'StopTime', tstopi);
        CVodeInit(der, 'BDF', 'Newton', tstart_i, x_stop, options);
    elseif isDiscontinuity
        % Reinitialization updates the states with their post-discontinuity
        % values and moves the final time to the next discontinuity
        %options = CVodeSetOptions(options, 'StopTime', tstopi);
        CVodeReInit(tstart_i, x_stop, options);
    end
    
    % Determine which t_get are in this interval of time
    tget_is_i = t_get > tstart_i & t_get <= tstop_i;
    
    % Set up list of time points to record during this interval of time
    % (includes t_get points determined above and next discontinuity)
    trecord_i = unique([vec(t_get(tget_is_i)); tstop_i]);
    trecord_is_tget = ismember(trecord_i, t_get(tget_is_i));
    trecord_is_tstop = trecord_i == tstop_i;
    
    % Integrate to next discontinuity or event
    [status, t_i, x_i] = CVode(trecord_i, 'Normal');
    
    isDiscontinuity = false;
    
    % Determine why integration stopped and handle accordingly
    switch status
        case -1 % error
            CVodeFree;
            error('KroneckerBio:accumulateOde:IntegrationFailure', 'Did not integrate through entire interval!');
        case 2 % event(s)
            
            % Get time and states of the last time step
            stats = CVodeGetStats;
            tcur = stats.tcur;
            xcur = CVodeGet('DerivSolution', tcur, 0);
    
            % Determine which events occurred
            eventOccurred = stats.RootInfo.roots; % eventOccurred is a vector with 1's at positions of events that triggered
            
            % Check directions (currently not implemented. Does anyone/anything use this?)
            
            % Record each event
            for ie_cur = find(eventOccurred(:).' == 1)
                if is_finished(sol) && isTerminal(tcur, xcur)
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
                sol.ie(nevents) = ie_cur;
                sol.te(nevents) = tcur;
                sol.xe(:,nevents) = xcur;
            end
        otherwise % discontinuity
            %if ismember(ti, discontinuities) % Discontinuity
            
            % Apply delta
            if ~isempty(delta)
                x_i(:,trecord_is_tstop) = x_i(:,trecord_is_tstop) + delta(tstop_i, x_stop);
            end
            x_stop = x_i(:,trecord_is_tstop);
            
            x_get = x_i(:,trecord_is_tget);
            sol.t(tget_is_i) = t_i(trecord_is_tget);
            sol.x(:,tget_is_i) = x_get;
            isDiscontinuity = true;
            %end
%             if ismember(ti, t_get) % Requested time
%                 % Record time and states
%                 sol.t(ii) = ti;
%                 sol.x(:,ii) = xi;
%             end
            if t_i(trecord_is_tstop) == tF % Check if reached tF
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