function [output_fun, event_fun] = timeoutOutputFunction(secondsToTimeout, otherOutputFun, otherEventFun)

if nargin < 3
    otherEventFun = [];
end
if nargin < 2
    otherOutputFun = [];
end

thistic = [];
output_fun = @outputfun;
event_fun = @eventfun;

    function stop = outputfun(t, y, flag)
        
        thisstop = false;
        
        % Run other output function, if provided
        if isempty(otherOutputFun)
            stop = false;
        else
            stop = otherOutputFcn(t,y,flag);
        end
        
        switch flag
            case 'init'
                % Start timer at beginning of integration
                thistic = tic;
            case ''
                % Stop if time elapsed exceeds the timeout duration
                thisstop = toc(thistic) > secondsToTimeout;
            case 'iter'
                % Stop if time elapsed exceeds the timeout duration
                thisstop = toc(thistic) > secondsToTimeout;
        end
    
        % Use an assertion to issue an error here, since setting the stop
        % flag only stops integration as an event. If the integration would
        % have stopped anyway, don't issue an error.
        assert(~thisstop || stop, 'KroneckerBio:timeoutOutputFunction:IntegrationTimeout', ...
            'A single ODE solver call took longer than %g seconds', secondsToTimeout)

        % Stop if timeout applies or if other output function issues a stop
        % condition
        stop = thisstop || stop;
    end

    function [value, isterminal, direction] = eventfun(t, y)
        
        % Stop if time elapsed exceeds the timeout duration. Skip if
        % thistic isn't initialized yet.
        thisstop = ~isempty(thistic) && toc(thistic) > secondsToTimeout;

        % Use an assertion to issue an error here, since setting the stop
        % flag only stops integration as an event
        assert(~thisstop, 'KroneckerBio:timeoutOutputFunction:IntegrationTimeout', ...
            'A single ODE solver call took longer than %g seconds', secondsToTimeout)
        
        % Run other event function, if provided
        if isempty(otherEventFun)
            % Put dummy values in that never throw an event to prevent errors
            value = 1;
            isterminal = 0;
            direction = 0;
        else
            [value, isterminal, direction] = otherEventFun(t,y);
        end
        
    end

end