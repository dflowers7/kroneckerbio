function [c,ceq,GC,GCeq] = sumobjectives_fmincon_constraint(obj, int)

nconstraints = numel(obj);
for i = nconstraints:-1:1
    if nargout <= 2
        [c(i,1),~] = sumobjectives_fmincon(obj{i}, int);
    else
        [c(i,1),~,GC(:,i),~] = sumobjectives_fmincon(obj{i}, int);
    end
end

ceq = zeros(0,1);
GCeq = zeros(size(GC,1),0);

end