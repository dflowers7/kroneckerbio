function [c,ceq,GC,GCeq,err,derrdT,Happrox]= concatenateerrors_fmincon_constraint(obj, ints, flipsign)

if isempty(ints)
    c = 0;
    GC = [];
    ceq = [];
    GCeq = [];
    return
end

[n_obj,n_con] = size(ints);

nT = numel(ints(1).TisExperiment);

c = cell(size(ints));
GC = cell(size(ints));
Happrox = [];
isObjectiveZero = strcmp({obj.Type}, 'Objective.Data.Zero');
isObjectiveZero = reshape(isObjectiveZero, size(obj));

for i_con = 1:n_con
    
    TisExperiment = ints(1,i_con).TisExperiment;
    
    for i_obj = find(~isObjectiveZero(:,i_con).')
        
        weight = sqrt(ints(i_obj,i_con).ObjWeight);
        
        err_temp = weight*obj(i_obj,i_con).err(ints(i_obj,i_con));
        nDataPoints_i = numel(err_temp);
        if nDataPoints_i > 0
            c{i_obj,i_con} = err_temp;
        end
        
        if nargout > 2
            derrdT_temp = weight*obj(i_obj,i_con).derrdT(ints(i_obj,i_con));
            
            if nDataPoints_i > 0
                GC{i_obj,i_con} = zeros(nT, nDataPoints_i);
                GC{i_obj,i_con}(TisExperiment,:) = derrdT_temp.';
            end
            
        end
        
    end
    
end

c = vertcat(c{:});
if nargout > 2
    GC = [GC{:}];
end

if flipsign
    c = -c;
    if nargout > 2
        GC = -GC;
    end
end

ceq = zeros(0,1);
if nargout > 2
    GCeq = zeros(size(GC,1),0);
end
if nargout > 4
    err = [];
    derrdT = [];
    Happrox = repmat({zeros(nT)}, size(c));
end

end
