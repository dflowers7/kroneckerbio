function [c,ceq,GC,GCeq] = sumobjectives_fmincon_constraint(obj, int)

% nconstraints = numel(obj);
% for i = nconstraints:-1:1
if nargout <= 2
    c = sumobjectives_fmincon(obj, int);
else
    [c,GC] = sumobjectives_fmincon(obj, int);
end
%end

ceq = zeros(0,1);
if nargout > 2
    GCeq = zeros(size(GC,1),0);
end

end