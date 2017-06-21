function sol = accumulateOdeRevSundials(options_x, index_b, t0, tf, t_get, discontinuities, xf, gf, nx, ng, delta_x, delta_g, freeMemoryOnFinish)
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
if nargin < 12
    freeMemoryOnFinish = true;
    if nargin < 11
        delta_g = [];
        if nargin < 10
            delta_x = [];
        end
    end
end

% Determine whether pure quadratures (g) are to be calculated
if isempty(ng); ng = 0; end; % Set empty values to 0
isQuad = ng > 0;

% Reverse sort t_get
t_get = unique(t_get);
t_get = t_get(end:-1:1);

% discontinuities will be a column vector sorted ascending. These are
% discontinuities in the output, not in the RHS. Note I don't have to
% include the final time point here because I already added in delta to
% xf in the initialize* function.
discontinuities = unique([vec(discontinuities); t0]);
discontinuities = discontinuities(end:-1:1); % Sort from last time to first
N = numel(t_get);

% Determine stopping times, which are the discontinuities and measured
% times
t_stop = unique([vec(discontinuities); vec(t_get)]);
t_stop = t_stop(end:-1:1);
%t_stop = discontinuities;

% Initialize solution vectors t and x
sol.t = zeros(1,N);
sol.x = zeros(nx,N);
if isQuad
    sol.g = zeros(ng, N);
end

% Handle special case where t_get includes the final time. CVODES errors
% if the first requested time is too close to the final time.
tget_is_tf = t_get == tf;
sol.t(tget_is_tf) = tf;
sol.x(:,tget_is_tf) = repmat(xf, 1, sum(tget_is_tf));
if isQuad
    sol.g(:, tget_is_tf) = repmat(gf, 1, sum(tget_is_tf));
end
t_get(tget_is_tf) = NaN; % Placeholder to ensure sizes of arrays match up, while also removing tF from list

% Iterate over the discontinuities, integrating to each one from the last
ii = 0;
tstop_i = tf; % Will be copied to tstart_i
done = false;
while ~done

    ii = ii + 1;
    tstart_i = tstop_i;
    tstop_i = t_stop(find(t_stop < tstart_i, 1, 'first'));
    
    %CVodeSetB(index_b, 'StopTime', tstop_i);
    
    % Reinitialization updates the states and quadratures with their
    % pre-discontinuity values and moves the final time to the previous
    % discontinuity
    if ii > 1
        % Calculate quadrature delta first, because when going backward in time,
        % the reported state value for a time t is post-delta. Therefore,
        % we need the post-delta value for x here.
        if isQuad
            g_delta = delta_g(tstart_i, [], xstop_i);
            if any(g_delta ~= 0)
                gstop_i = gstop_i - g_delta;
                CVodeQuadReInitB(index_b, gstop_i);
            end
        end
        
        x_delta = delta_x(tstart_i, [], xstop_i);
        %if any(x_delta ~= 0)
            xstop_i = xstop_i - x_delta;
            CVodeReInitB(index_b, tstart_i, xstop_i, options_x);
        %end
    end
    
%     % Determine which t_get are in this interval of time
%     tget_is_i = t_get < tstart_i & t_get >= tstop_i;
%     
%     % Set up list of time points to record during this interval of time
%     % (includes t_get points determined above and next discontinuity)
%     trecord_i = unique([vec(t_get(tget_is_i)); tstop_i]);
%     trecord_i = trecord_i(end:-1:1);
%     trecord_is_tget = ismember(trecord_i, t_get(tget_is_i));
%     trecord_is_tstop = trecord_i == tstop_i;
    
    % Integrate to previous discontinuity or event
    if isQuad
        [status, t_i, x_i, x_q] = CVodeB(tstop_i, 'Normal');
    else
        [status, t_i, x_i] = CVodeB(tstop_i, 'Normal');
    end
    
    % Determine why integration stopped and handle accordingly
    switch status
        case -1 % error
            
            CVodeFree;
            error('KroneckerBio:accumulateOde:IntegrationFailure', 'Did not integrate through entire interval!');
            
        otherwise % discontinuity or measured time
            
            % Update values of states and sensitivities at tstop_i
            xstop_i = x_i;
            gstop_i = x_q;
            
            % Store value if requested
            tget_is_tstop = tstop_i == t_get;
            if any(tget_is_tstop)
                sol.t(tget_is_tstop) = t_i;
                sol.x(:,tget_is_tstop) = x_i;
                if isQuad
                    sol.g(:,tget_is_tstop) = gstop_i;
                end
            end
            
            if t_i == t0 % Check if reached t0
                done = true;
            end
    end
    
end

% Free memory
if freeMemoryOnFinish
    CVodeFree;
end

end