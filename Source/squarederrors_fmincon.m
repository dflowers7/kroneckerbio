function [err2,derr2dT_transpose] = squarederrors_fmincon(obj, ints)

if isempty(ints)
    err2 = [];
    derr2dT_transpose = [];
    return
end

[n_obj,n_con] = size(ints);

nT = numel(ints(1).TisExperiment);

err = nan(100,1);
derrdT_transpose = nan(nT,100);
nGfilled = 0;

for i_con = 1:n_con
    
    TisExperiment = ints(1,i_con).TisExperiment;
    
    for i_obj = 1:n_obj
        
        G_i = obj(i_obj,i_con).err(ints(i_obj,i_con));
        nGadding = numel(G_i);
        
        if nGfilled + nGadding > numel(err)
            nG_resize = max(numel(err), nGadding);
            err = [err; nan(nG_resize,1)];
            derrdT_transpose = [derrdT_transpose nan(nT, nG_resize)];
        end
        
        err(nGfilled+1:nGfilled+nGadding) = G_i;
        
        if nargout > 1
            
            nerr_i = numel(G_i);
            D_i = zeros(nT, nerr_i);
            D_i(TisExperiment,:) = obj(i_obj,i_con).derrdT(ints(i_obj,i_con)).';
            % Don't need to normalize because integration already
            % normalized the values of dydT
%             if ints(i_obj,i_con).Normalized
%                 D_i = D_i.*ints(i_obj,i_con).T;
%             end

            derrdT_transpose(:,nGfilled+1:nGfilled+nGadding) = D_i;
        end
        
        nGfilled = nGfilled + nGadding;
        
    end
    
end

err = err(~isnan(err));
derrdT_transpose = derrdT_transpose(:,~isnan(derrdT_transpose(1,:)));

err2 = err.^2;
if nargout > 1
    derr2dT_transpose = 2*bsxfun(@times, err(:).', derrdT_transpose);
end

end