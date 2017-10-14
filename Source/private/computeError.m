function [G,D,Hv] = computeError(m, con, obj, opts, v)
% [G,D,Hv] = computeError(m, con, obj, opts, v) 
% This function computes the vector of errors G for an objective function
% obj using the simulation defined by the model m and the experimental
% conditions con. It can also calculate the gradient of the errors D,
% returned as an nErrors-by-nT matrix, and the directional second
% derivative of the errors Hv in the direction of nT-by-1 vector v,
% returned as an nErrors-by-nT matrix.
verbose = logical(opts.Verbose);
verbose_all = max(verbose-1,0);
calcJac = nargout > 1;
calcHv = nargout > 2;
assert(nargin > 4 || ~calcHv, 'KroneckerBio:computeError:HvRequiresv', ...
    'Computing a directional second derivative requires the direction v be provided.')

if calcHv
    normv = norm(v);
    dirv = v/normv;
    stepsize = 1e-8;
end

% Constants
n_con = numel(con);
n_obj = size(obj,1);

G = zeros(100,1);
if calcJac
    nTk = nnz(opts.UseParams);
    nTs = nnz(opts.UseSeeds);
    nTq = nnz(cat(1,opts.UseInputControls{:}));
    nTh = nnz(cat(1,opts.UseDoseControls{:}));
    nT  = nTk + nTs + nTq + nTh;
    D = zeros(100,nT);
end

% Need to collect fit parameters for all experiments for priors
T_allcon = collectFitParameters(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

if verbose; disp('Integrating forward...'); end
count = 1;
for i_con = 1:n_con
    if verbose_all; tic; end
    
    % Modify opts structure
    opts_i = opts;
    
    UseSeeds_i = opts.UseSeeds(:,i_con);
    opts_i.UseSeeds = UseSeeds_i;
    
    UseInputControls_i = opts.UseInputControls{i_con};
    opts_i.UseInputControls = UseInputControls_i;
    
    UseDoseControls_i = opts.UseDoseControls{i_con};
    opts_i.UseDoseControls = UseDoseControls_i;
    
    opts_i.AbsTol = opts.AbsTol{i_con};
    opts_i.ObjWeights = opts.ObjWeights(:,i_con);
    
    % Integrate
    if nargout == 1
        ints = integrateAllSys(m, con(i_con), obj(:,i_con), opts_i);
    elseif nargout == 2
        ints = integrateAllSens(m, con(i_con), obj(:,i_con), opts_i);
    end
    
    if calcHv
        inds = ParameterMappings(opts);
        
        % Update model parameters
        dk = zeros(m.nk,1);
        dk(inds.k) = dk(inds.k) + 1i*stepsize*dirv(inds.Tk);
        mv = m.Update(m.k + dk);
        
        % Update experiment parameters
        ds = zeros(m.ns,1);
        ds(inds.s{i_con}) = 1i*stepsize*dirv(inds.Ts{i_con});
        dq = zeros(con(icon).nq,1);
        dq(inds.q{i_con}) = 1i*stepsize*dirv(inds.Tq{i_con});
        dh = zeros(con(icon).nh,1);
        dh(inds.h{i_con}) = 1i*stepsize*dirv(inds.Th{i_con});
        conv = con(i_con).Update(con(i_con).s+ds, con(i_con).q+dq, con(i_con).h+dh);
        
        intsv = integrateAllSens(mv, conv, obj(:,i_con), opts_i);
    end
    
    % Add fields for prior objectives
    [ints.UseParams] = deal(opts.UseParams);
    [ints.UseSeeds] = deal(UseSeeds_i);
    [ints.UseInputControls] = deal(UseInputControls_i);
    [ints.UseDoseControls] = deal(UseDoseControls_i);
    [ints.T] = deal(T_allcon);
    
    if calcJac
        [ints.nT] = deal(nT);
    end
    
    count_start_icon = count;
    for i_obj = 1:n_obj
        G_i = obj(i_obj,i_con).err(ints(i_obj));
        findex = count+numel(G_i)-1;
        addRows = findex > numel(G);
        if addRows
            factor2 = ceil(log2(findex/numel(G)));
            addsize = numel(G)*(2^factor2-1);
            G = [G; zeros(addsize,1)];
        end
        G(count:findex) = G(count:findex) + opts.ObjWeights(i_obj,i_con) .* G_i;
        
        if calcJac
            D_i = obj(i_obj,i_con).derrdT(ints(i_obj));
            if addRows
                D = [D; zeros(addsize,nT)];
            end
            D(count:findex,:) = D(count:findex,:) + opts.ObjWeights(i_obj,i_con) .* D_i;
        end
        
        if calcHv
            Hv_i = obj(i_obj,i_con).derrdT(intsv(i_obj));
        end
        
        count = count + numel(G_i);
    end
    
    if verbose_all; fprintf('iCon = %d\tResidual = %g\tTime = %0.2f\n', i_con, sum(G(count_start_icon:findex).^2), toc); end
end

G = G(1:findex);
if calcJac
    D = D(1:findex,:);
end

if verbose; fprintf('Summary: Residual = %g\n', sum(G.^2)); end
