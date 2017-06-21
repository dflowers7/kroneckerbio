function [m,con,rank,eigval,eigvals,F] = MaximizeFIMRank(m,con,obj,opts,eigvalthreshold_terminalobj,eigvalthreshold_stopmaximizing,scalebyresidual)

if nargin < 7
    scalebyresidual = [];
    if nargin < 6
        eigvalthreshold_stopmaximizing = [];
        if nargin < 5
            eigvalthreshold_terminalobj = [];
        end
    end
end

if isempty(scalebyresidual)
    scalebyresidual = true;
end
if isempty(eigvalthreshold_terminalobj)
    eigvalthreshold_terminalobj = 1e-6;
end
if isempty(eigvalthreshold_stopmaximizing)
    eigvalthreshold_stopmaximizing = eigvalthreshold_terminalobj/1e3;
end

assert(eigvalthreshold_terminalobj >= eigvalthreshold_stopmaximizing, ...
    'KroneckerBio:MaximizeFIMRank:Thresholds', ...
    'The terminal objective value threshold should be greater than or equal to the threshold for ceasing maximization.')

% Calculate the FIM and the residual
[F,res] = evaluateObjective(m,con,obj,opts);

% Calculate the FIM eigenvalues and sort them
eigF = eig(F);
eigF = sort(eigF, 1, 'descend');
if scalebyresidual
    eigF = eigF./res;
end

% Find the index of the first eigenvalue that is less than the threshold
% value
eigindex_start = find(eigF < eigvalthreshold_terminalobj, 1);
rank = eigindex_start-1;
if eigindex_start == 1
    eigval = eigF(eigindex_start);
else
    eigval = eigF(eigindex_start-1);
end

nT = size(F,1);
neigstoincrease = nT - eigindex_start + 1;

% Remove terminal objective function value and replace with threshold value
% if isfield(opts, 'TerminalObj')
%     opts = rmfield(opts, 'TerminalObj');
% end
opts.TerminalObj = -log(eigvalthreshold_terminalobj);

% Use the 'interior-point' algorithm to avoid bfgs approximation to
% Hessian. bfgs doesn't work well here when the order of the eigenvalues
% shifts due to parameter changes, causing discontinuities in the Hessian
% opts.Algorithm = 'interior-point';
% opts.HessianApproximation = 'fin-diff-grads';
% opts.SubproblemAlgorithm = 'cg';

for i = 1:neigstoincrease
    eigindex_i = eigindex_start + i - 1;
    [intfun,objfun] = GenerateFIMEigenvalueFunction(eigindex_i,false,scalebyresidual);
    opts_i = opts;
    opts_i.IntegrateFunction = intfun;
    opts_i.ObjectiveReductionFunction = objfun;
    fprintf('Maximizing eigenvalue %g of %g...\n', eigindex_i, nT)
    [m_i,con_i,neglogeigval_i] = FitObjective(m,con,obj,opts_i);
    eigval_i = exp(-neglogeigval_i);
    if eigval_i < eigvalthreshold_stopmaximizing
        fprintf('Eigenvalue index %g maximized at a value of %g, short of the stop threshold %g. Stopping maximization.\n', eigindex_i, eigval_i, eigvalthreshold_stopmaximizing)
        m = m_i;
        con = con_i;
        break
    else
        fprintf('Eigenvalue index %g maximized to a value of %g. Continuing maximization of next eigenvalue.\n', eigindex_i, eigval_i)
        m = m_i;
        con = con_i;
        rank = eigindex_i;
        eigval = eigval_i;
    end
end

if nargout > 4
    [F,res] = evaluateObjective(m,con,obj,opts);
    eigvals = eig(F);
    if scalebyresidual
        eigvals = eigvals./res;
    end
end

end

function [F,res] = evaluateObjective(m,con,obj,opts)

sim = SimulateSensitivity(m, con, obj, opts);
int = cell(size(sim,2),1);
res = 0;
for i = size(sim,2):-1:1
    int{i} = vertcat(sim(:,i).int);
    for j = size(sim,1):-1:1
        int_ij = int{i}(j);
        err = obj(j,i).err(int_ij);
        res = res + err.'*err;
    end
end
F = ObjectiveInformation(m,con,obj,opts,int);

end