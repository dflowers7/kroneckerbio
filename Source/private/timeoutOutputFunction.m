function fun = timeoutOutputFunction(secondsToTimeout, otherOutputFun)

if nargin < 2
    otherOutputFun = [];
end

thistic = [];
fun = @outputfun;

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

end