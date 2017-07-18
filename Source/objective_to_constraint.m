function [c,ceq,GC,GCeq,err,derrdT,Happrox] = objective_to_constraint(obj_red_fun, obj, int)

% nconstraints = numel(obj);
% for i = nconstraints:-1:1
if nargout <= 2
    c = obj_red_fun(obj, int);
elseif nargout <= 4
    [c,GC] = obj_red_fun(obj, int);
else
    [c,GC,err,derrdT,Happrox] = obj_red_fun(obj, int);
end
%end

ceq = zeros(0,1);
if nargout > 2
    GCeq = zeros(size(GC,1),0);
end

end