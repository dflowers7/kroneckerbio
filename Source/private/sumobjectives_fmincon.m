function [G,D] = sumobjectives_fmincon(obj, ints)

if isempty(ints)
    G = 0;
    D = [];
    return
end

[n_obj,n_con] = size(ints);

nT = numel(ints(1).TisExperiment);

G = 0;
D = zeros(nT,1);

for i_con = 1:n_con
    
    TisExperiment = ints(1,i_con).TisExperiment;
    
    for i_obj = 1:n_obj
        
        G_i = obj(i_obj,i_con).G(ints(i_obj,i_con));
        G = G + ints(i_obj,i_con).ObjWeight * G_i;
        
        if nargout > 1
            D_i = obj(i_obj,i_con).dGdT(ints(i_obj,i_con));
            % Don't need to normalize because integration already
            % normalized the values of dydT
%             if ints(i_obj,i_con).Normalized
%                 D_i = D_i.*ints(i_obj,i_con).T;
%             end

            D(TisExperiment) = D(TisExperiment) + ints(i_obj,i_con).ObjWeight * D_i;
        end
        
    end
    
end

end