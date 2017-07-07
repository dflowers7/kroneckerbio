function [G,D,err,derrdT,Happrox] = sumobjectives_fmincon(obj, ints)

if isempty(ints)
    G = 0;
    D = [];
    return
end

[n_obj,n_con] = size(ints);

nT = numel(ints(1).TisExperiment);

G = 0;
D = zeros(nT,1);
err = cell(size(ints));
derrdT = cell(size(ints));
Happrox = zeros(nT);
isObjectiveZero = strcmp({obj.Type}, 'Objective.Data.Zero');
isObjectiveZero = reshape(isObjectiveZero, size(obj));
isPureLeastSquares = true;

for i_con = 1:n_con
    
    TisExperiment = ints(1,i_con).TisExperiment;
    
    for i_obj = 1:n_obj
        
        ObjWeight = ints(i_obj,i_con).ObjWeight;
        
        
        G_i = obj(i_obj,i_con).G(ints(i_obj,i_con));
        G = G + ObjWeight * G_i;
        
        if nargout > 1
            D_i = obj(i_obj,i_con).dGdT(ints(i_obj,i_con));
            % Don't need to normalize because integration already
            % normalized the values of dydT
%             if ints(i_obj,i_con).Normalized
%                 D_i = D_i.*ints(i_obj,i_con).T;
%             end

            D(TisExperiment) = D(TisExperiment) + ints(i_obj,i_con).ObjWeight * D_i;
            
            if nargout > 2 && isPureLeastSquares
                ObjWeight_sqrt = sqrt(ObjWeight);
                err_temp = ObjWeight_sqrt*obj(i_obj,i_con).err(ints(i_obj,i_con));
                derrdT_temp = ObjWeight_sqrt*obj(i_obj,i_con).derrdT(ints(i_obj,i_con));
                nDataPoints_i = numel(err_temp);
                if nDataPoints_i > 0
                    err{i_obj,i_con} = err_temp;
                    derrdT{i_obj,i_con} = zeros(nDataPoints_i, nT);
                    derrdT{i_obj,i_con}(:,TisExperiment) = derrdT_temp;
                end
            end
            if nargout > 4
                Happrox(TisExperiment,TisExperiment) = Happrox(TisExperiment,TisExperiment) ...
                    + ObjWeight*obj(i_obj,i_con).d2GdT2_approximate(ints(i_obj,i_con));
            end
        end
        
    end
    
end

if nargout > 2 && isPureLeastSquares
    err = vertcat(err{:});
    derrdT = vertcat(derrdT{:});
end

end
