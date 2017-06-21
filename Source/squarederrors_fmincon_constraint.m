function [c,ceq,GC,GCeq,derrdT] = squarederrors_fmincon_constraint(obj, int)

% nconstraints = numel(obj);
% for i = nconstraints:-1:1
if nargout <= 2
    c = squarederrors_fmincon(obj, int);
elseif nargout <= 4
    [c,GC] = squarederrors_fmincon(obj, int);
else
    [c,GC,derrdT] = squarederrors_fmincon(obj,int);
end
%end

ceq = zeros(0,1);
if nargout > 2
    GCeq = zeros(size(GC,1),0);
end

end