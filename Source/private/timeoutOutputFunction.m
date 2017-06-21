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
            case 'iter'
                % Stop if time elapsed exceeds the timeout duration
                thisstop = toc(thistic) > secondsToTimeout;
        end
        
        % Stop if timeout applies or if other output function issues a stop
        % condition
        stop = thisstop || stop;
    end

end