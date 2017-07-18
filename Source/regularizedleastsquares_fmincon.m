function [G,D,err,derrdT,Happrox] = regularizedleastsquares_fmincon(obj, ints, minimumSingularValue, Tcenter, goodsumsqerrs)
% [G,D,err,derrdT,Happrox] = regularizedleastsquares_fmincon(obj, ints, minimumSingularValue)
% Use this objective function reduction function by creating an anonymous
% function setting minimumSingularValue, Tcenter, and goodsumsqerrs.

if nargin < 5
    goodsumsqerrs = [];
    if nargin < 4
        Tcenter = [];
    end
end

if isempty(Tcenter)
    Tcenter = 1;
end

if isempty(ints)
    G = 0;
    D = [];
    return
end

if nargout <= 1
    G = sumobjectives_fmincon(obj, ints);
else
    [G,D,err,derrdT,Happrox] = sumobjectives_fmincon(obj, ints);
end

if nargout > 1
    
    nT = numel(ints(1).TisExperiment);
    T = nan(nT,1);
    for i = 1:numel(ints)
        T(ints(i).TisExperiment) = ints(i).T;
    end
    
    if isscalar(Tcenter)
        Tcenter = repmat(Tcenter,nT,1);
    end
    
    isnormalized = ints(1).Normalized;
    if isnormalized
        logTcenter = log(Tcenter);
        logT = log(T);
    else
        logTcenter = Tcenter;
        logT = T;
    end
    
    % Set derrdT singular values to minimum, if they are below the minimum
    derrdT_orig = derrdT;
    [u,s,v] = svd(derrdT);
    sumsqerrs = sum(err.^2);
    nDataPoints = numel(err);
    if isempty(goodsumsqerrs)
        goodsumsqerrs = chi2inv(0.95, nDataPoints);
    end
    if sumsqerrs > goodsumsqerrs
        % Change the singular values of derrdT only if the fit isn't good enough, i.e.,
        % the sum of squares is above the 95 percent upper confidence
        % limit on the chi-squared statistic
        for i = 1:min(size(derrdT))
            s(i,i) = max(minimumSingularValue, s(i,i));
%             if s(i,i) <= minimumSingularValue
%                 % Choose signs on small singular vectors such that the gradient points more towards the
%                 % center of the parameter space
%                 ssign = sign((logT-logTcenter).'*v(:,i)*s(i,i)*u(:,i).'*err);
%                 s(i,i) = ssign*s(i,i);
%             end
        end
        
        derrdT = u*s*v.';
    end
    
    % Recalculate the other values by subtracting off old derrdT contribution
    % and adding new derrdT contribution
    D = D - 2*(derrdT_orig - derrdT).'*err; % Subtract off old derrdT contribution and add new one
    Happrox = Happrox + 2*(-(derrdT_orig.'*derrdT_orig) + (derrdT.'*derrdT));
end

end