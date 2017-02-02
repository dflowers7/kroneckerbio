function [err,derrdT] = concatobjectives_lsqnonlin(obj, ints)

if isempty(ints)
    err = [];
    derrdT = [];
    return
end

[n_obj,n_con] = size(ints);

nT = numel(ints(1).TisExperiment);

err_cell = cell(n_obj,n_con);
derrdT_cell = cell(n_obj,n_con);

for i_con = 1:n_con
    
    TisExperiment = ints(1,i_con).TisExperiment;
    
    for i_obj = 1:n_obj
        
        err_cell{i_obj,i_con} = obj(i_obj,i_con).err(ints(i_obj,i_con));
        
        if nargout > 1
            derrdT_i = obj(i_obj,i_con).derrdT(ints(i_obj,i_con));
            if ints(i_obj,i_con).Normalized
                derrdT_i = bsxfun(@times, derrdT_i, ints(i_obj,i_con).T(:).');
            end

            derrdT_cell{i_obj,i_con} = zeros(size(derrdT_i,1), nT);
            derrdT_cell{i_obj,i_con}(:,TisExperiment) = derrdT_i;
        end
        
    end
    
end

err = vertcat(err_cell{:});
derrdT = vertcat(derrdT_cell{:});

end