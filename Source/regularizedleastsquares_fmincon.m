function [G,D,err,derrdT,Happrox] = regularizedleastsquares_fmincon(obj, ints, minimumSingularValue)
% [G,D,err,derrdT,Happrox] = regularizedleastsquares_fmincon(obj, ints, minimumSingularValue)
% Use this objective function reduction function by creating an anonymous
% function setting minimumSingularValue

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
    % Set derrdT singular values to minimum, if they are below the minimum
    derrdT_orig = derrdT;
    [u,s,v] = svd(derrdT);
    sumerrs = sum(abs(err));
    for i = 1:min(size(derrdT))
        s(i,i) = max(minimumSingularValue, s(i,i));
    end
    derrdT = u*s*v.';
    
    % Recalculate the other values by subtracting off old derrdT contribution
    % and adding new derrdT contribution
    D = D - 2*(derrdT_orig - derrdT).'*err; % Subtract off old derrdT contribution and add new one
    Happrox = Happrox + 2*(-(derrdT_orig.'*derrdT_orig) + (derrdT.'*derrdT));
end

end