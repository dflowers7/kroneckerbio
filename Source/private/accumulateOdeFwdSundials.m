function sol = accumulateOdeFwdSundials(options, tF, t_get, discontinuities, x0, dx0dT, g0, nx, nT, ng, delta, delta_sens, events, is_finished, freeMemoryOnFinish)
% Inputs
%   nT: [ integer scalar ]
%       Number of parameters for which sensitivities are to be calculated
%       via the forward method. Set to 0 if forward method sensitivities
%       are not to be calculated.
%   ng: [ integer scalar ]
%       Number of pure quadratures that are to be integrated. Set to 0 if
%       no quadratures are to be calculated.
% Outputs
%   sol: struct with the following fields:
%       .t [ 1-by-nt double ]
%           Times at which an output value was desired. This is equal to
%           t_get for successful integrations unless a terminal event
%           occurs before tF, in which case all 
%       .x [ nx-by-nt double ]
%       .ie [ 1-by-nevents double ]
%           Index indicating which event occurred. Indices are determined
%           according to the order of events specified in the events input
%           argument.
%       .te [ 1-by-nevents double ]
%           Times at which events occurred.
%       .xe [ nx-by-nevents double ]
%           Values of states at each event time.

% Clean up inputs
if nargin < 11
    freeMemoryOnFinish = true;
    if nargin < 10
        is_finished = [];
        if nargin < 9
            events = [];
            if nargin < 8
                delta_sens = [];
                if nargin < 7
                    delta = [];
                end
            end
        end
    end
end

% Determine whether pure quadratures (g) or sensitivities are to be calculated
if isempty(ng); ng = 0; end; % Set empty values to 0
if isempty(nT); nT = 0; end;
isQuad = ng > 0;
isSens = nT > 0;

% Set default events functions
if isempty(events)
    events = @(t,x) [];
end
if isempty(is_finished)
    is_finished = @(sol) true; % Set to true so that the isterminal value set in the event function determines whether to halt integration
end

stats = CVodeGetStats;
t0 = stats.tcur;

% discontinuities will be a column vector sorted ascending. These are
% discontinuities in the output, not in the RHS. Note I don't have to
% include the initial time point here because I already added in delta to
% x0 in the initialize* function.
discontinuities = unique([vec(discontinuities); tF]);
N = numel(t_get);

% Determine stopping times, which are the discontinuities
%t_stop = unique([vec(discontinuities); vec(t_get)]);
t_stop = discontinuities;

% Initialize solution vectors t and x
sol.t = zeros(1,N);
sol.x = zeros(nx,N);
if isQuad
    sol.g = zeros(ng, N);
end
if isSens
    sol.dxdT = zeros(nx, nT, N);
end

% Initialize event solution vectors te, xe, and ye
eventArraySize = 10; % maximum number of events that can be stored in the initialized array
sol.ie = zeros(1,eventArraySize);
sol.te = zeros(1,eventArraySize);
sol.xe = zeros(nx,eventArraySize);
nevents = 0; % number of events that have occurred

% Handle special case where t_get includes the initial time. CVODES errors
% if the first requested time is too close to the initial time.
tget_is_t0 = t_get == t0;
sol.t(tget_is_t0) = t0;
sol.x(:,tget_is_t0) = x0 + delta(t0, x0);
if isSens
    sol.dxdT(:,:,tget_is_t0) = dx0dT + delta_sens(t0, x0);
end
if isQuad
    sol.g(:, tget_is_t0) = g0;
end
t_get(tget_is_t0) = NaN; % Placeholder to ensure sizes of arrays match up, while removing t0 from list

% Iterate over the discontinuities, integrating to each one from the last
ii = 0;
tstop_i = t0;
isDiscontinuity = false;
done = false;
while ~done
    
    ii = ii + 1;
    
    tstart_i = tstop_i;
    tstop_i = t_stop(find(t_stop > tstart_i, 1, 'first'));
    
    % Determine which t_get are in this interval of time
    tget_is_i = t_get > tstart_i & t_get <= tstop_i;
    
    % Set up list of time points to record during this interval of time
    % (includes t_get points determined above and next discontinuity)
    trecord_i = unique([vec(t_get(tget_is_i)); tstop_i]);
    trecord_is_tget = ismember(trecord_i, t_get(tget_is_i));
    trecord_is_tstop = trecord_i == tstop_i;
    
    if isDiscontinuity
        % Reinitialization updates the states and sensitivities with their
        % post-discontinuity values and moves the final time to the next
        % discontinuity
        %options = CVodeSetOptions(options, 'StopTime', tstopi);
        CVodeReInit(tstart_i, x_stop, options);
        if isSens && ~isempty(delta_sens)
            CVodeSensReInit(x_s_stop)
        end
    end
    
    % Integrate to next discontinuity or event
    if ~isQuad && ~isSens
        [status, t_i, x_i] = CVode(trecord_i, 'Normal');
    elseif isQuad && isSens
        [status, t_i, x_i, x_q, x_s] = CVode(trecord_i, 'Normal');
    elseif isQuad
        [status, t_i, x_i, x_q] = CVode(trecord_i, 'Normal');
    elseif isSens
        [status, t_i, x_i, x_s] = CVode(trecord_i, 'Normal');
    end
    
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
            
            % Check for discontinuity, in case we ended exactly on a time
            % point where a discontinuity occurs (to be implemented)
            
        otherwise % discontinuity
            
            % Apply delta
            if ~isempty(delta)
                x_i(:,trecord_is_tstop) = x_i(:,trecord_is_tstop) + delta(tstop_i, x_i(:,trecord_is_tstop));
            end
            if isSens && ~isempty(delta_sens)
                x_s(:,:,trecord_is_tstop) = x_s(:,:,trecord_is_tstop) + delta_sens(tstop_i, x_i(:,trecord_is_tstop));
            end
            
            % Update values of states and sensitivities at tstop
            x_stop = x_i(:,trecord_is_tstop);
            if isSens
                x_s_stop = x_s(:,:,trecord_is_tstop);
            end
            
            % Store requested values
            sol.t(tget_is_i) = t_i(trecord_is_tget);
            x_get = x_i(:,trecord_is_tget);
            sol.x(:,tget_is_i) = x_get;
            if isQuad
                x_q_get = x_q(:,trecord_is_tget);
                sol.x_q(:,tget_is_i) = x_q_get;
            end
            if isSens
                x_s_get = x_s(:,:,trecord_is_tget);
                sol.dxdT(:, :, tget_is_i) = x_s_get;
            end
            
            isDiscontinuity = true;
            
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
if freeMemoryOnFinish
    CVodeFree;
end

    function val = isTerminal(t, x)
        [~,val] = events(t, x);
        val = val == 1;
    end

    % Currently unused. Might be useful if/when event directions are
    % implemented.
    function val = directionIsSatisfied(t, x, xprev)
        % requiredDirections: 0 means all directions trigger, 1 means
        % increasing function value triggers, -1 means decreasing function
        % value triggers
        [~,~,requiredDirections] = events(t,x);
        actualDirections = sign(x-xprev);
        val = requiredDirections(:) == 0 | actualDirections(:) == requiredDirections(:);
    end
end