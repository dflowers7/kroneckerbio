function [objfun, constrfun] = GenerateObjective(m, con, obj, opts, integrateFunctions_objective, objectiveFunction, integrateFunctions_constraint, constraintFunctions)

if nargin < 8
    constraintFunctions = [];
    if nargin < 7
        integrateFunctions_constraint = [];
        if nargin < 6
            objectiveFunction = [];
            if nargin < 5
                integrateFunctions_objective = [];
            end
        end
    end
end

% Set default integrate, objective, and constraint functions if no custom
% functions are provided
if isempty(integrateFunctions_objective)
    integrateFunctions_objective = {@integrateAllSys; @integrateAllSens};
end

if isempty(objectiveFunction)
    switch opts.Solver
        case 'fmincon'
            objectiveFunction = @sumobjectives_fmincon;
        case 'lsqnonlin'
            objectiveFunction = @concatobjectives_lsqnonlin;
    end
end

nconstraints = numel(opts.ConstraintVal);

if isempty(integrateFunctions_constraint)
    integrateFunctions_constraint = repmat({{@integrateAllSys;@integrateAllSys;@integrateAllSens;@integrateAllSens}}, nconstraints, 1);
end

if isempty(constraintFunctions)
    switch opts.Solver
        case 'fmincon'
            constraintFunctions = repmat({@sumobjectives_fmincon_constraint}, nconstraints, 1);
        case 'lsqnonlin'
            constraintFunctions = cell(0,1);
    end
end

assert(~strcmp(opts.Solver,'lsqnonlin') || nconstraints == 0, ...
    'KroneckerBio:GenerateObjective:NonlinearConstraintWithLsqnonlin', ...
    'lsqnonlin does not support nonlinear constraint functions. Use fmincon instead.')

nk = m.nk;
ns = m.ns;

% Ensure structures are proper sizes
[con, n_con] = fixCondition(con);
[obj, n_obj] = fixObjective(obj, n_con);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, n_con);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, n_con, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, n_con, cat(1,con.nh));

nT = nTk + nTs + nTq + nTh;

% Extract constraint objectives
constr_obj = opts.ConstraintObj;
constr_vals = opts.ConstraintVal;
hasconstraint = ~isempty(constr_vals);

% Determine which parameters are fit by which experiments
T_experiment = zeros(nT, 1);
T_experiment(1:nTk) = 0;
nTs_con = sum(opts.UseSeeds,1);
nTq_con = cellfun(@sum, opts.UseInputControls(:)');
nTh_con = cellfun(@sum, opts.UseDoseControls(:)');
nT_con = {nTs_con, nTq_con, nTh_con};
endi = nTk;
for i_type = 1:3
    for j_con = 1:n_con
        starti = endi + 1;
        endi = endi + nT_con{i_type}(j_con);
        T_experiment(starti:endi) = j_con;
    end
end
% Determine which T's are fit by which experiments
TisExperiment = bsxfun(@(T_experiment,i_con)T_experiment == 0 | T_experiment == i_con, T_experiment, 1:n_con);

% Check for parallel toolbox if parallel optimization is specified
if opts.ParallelizeExperiments && isempty(ver('distcomp'))
    warning('KroneckerBio:FitObjective:ParallelToolboxMissing', ...
        ['opts.ParallelizeExperiments requires the Parallel '...
        'Computing toolbox. Reverting to serial evaluation.'])
    opts.ParallelizeExperiments = false;
end

% Determine number of workers in the parallel pool
if opts.ParallelizeExperiments
    p = gcp();
    if isempty(p)
        warning('KroneckerBio:FitObjective:NoParallelPool', ...
            ['opts.ParallelizeExperiments was set to true, ' ...
            'but no parallel pool could be initialized. '...
            'Reverting to serial optimization.']);
        opts.ParallelizeExperiments = false;
    else
        NumWorkers = p.NumWorkers;
    end
else
    NumWorkers = 1;
end

% Split experiment indices into worker groups
baseNumPerWorker = floor(n_con/NumWorkers);
remainder = n_con - baseNumPerWorker*NumWorkers;
numPerWorker = repmat(baseNumPerWorker, NumWorkers, 1);
numPerWorker(1:remainder) = numPerWorker(1:remainder) + 1;
[icon_worker,cons,objs,optss,TisWorker] = deal(cell(NumWorkers, 1));
endi = 0;
for ii = 1:NumWorkers
    starti = endi+1;
    endi = sum(numPerWorker(1:ii));
    icon_worker{ii} = (starti:endi)';
    [cons{ii}, objs{ii}, optss{ii}, TisWorker{ii}] = splitExperiments(con, obj, opts, icon_worker{ii}, T_experiment);
end

% Distribute experiments to workers
if opts.ParallelizeExperiments
    spmd
        cons = codistributed(cons);
        objs = codistributed(objs);
        optss = codistributed(optss);
        TisWorker = codistributed(TisWorker);
    end
end

% Set up variables for memoization
[lastT,lastint] = deal(struct('integrateAllSys', {[]}, 'integrateAllSens', {[]}, 'integrateAllCurvVecProd', {[]}));

displayedMemoWarning = false;

objfun = @objective;
if hasconstraint
    constrfun = @constraint;
else
    % Set to empty so that fmincon knows there is no constraint
    constrfun = [];
end

    function int = integr(T, integrateFunctions)
        
        nIntegrates = sum(~cellfun(@isempty,integrateFunctions));
        
        prevint = [];
        
        for j = 1:nIntegrates
            
            int_temp = get_memoized_int(T, integrateFunctions{j});
            if ~isempty(int_temp)
                int = int_temp;
                prevint = int;
                continue
            end
            
            if opts.ParallelizeExperiments
                spmd
                    % Get local experiments
                    con_i = getLocalPart(cons);
                    con_i = con_i{1};
                    obj_i = getLocalPart(objs);
                    obj_i = obj_i{1};
                    opts_i = getLocalPart(optss);
                    opts_i = opts_i{1};
                    TisWorker_i = getLocalPart(TisWorker);
                    TisWorker_i = TisWorker_i{1};
                    
                    if isempty(con_i)
                        ints = [];
                    else
                        % Update model and experiments with provided parameters
                        [m_i,con_i] = updateAll(m, con_i, T(TisWorker_i), opts_i.UseParams, opts_i.UseSeeds, opts_i.UseInputControls, opts_i.UseDoseControls);
                        
                        % Integrate the system
                        ncon_i = numel(con_i);
                        for i = ncon_i:-1:1
                            opts_ii = opts_i;
                            opts_ii.AbsTol = opts_ii.AbsTol{i};
                            opts_ii.UseSeeds = opts_ii.UseSeeds(:,i);
                            opts_ii.UseInputControls = opts_ii.UseInputControls{i};
                            opts_ii.UseDoseControls = opts_ii.UseDoseControls{i};
                            if j == 1
                                ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_i(:,i), opts_ii);
                            else
                                ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_i(:,i), opts_ii, obj, prevint);
                            end
                        end
                    end
                end
            else
                % Get local experiments
                con_i = cons{1};
                obj_i = objs{1};
                opts_i = optss{1};
                TisWorker_i = TisWorker{1};
                
                if isempty(con_i)
                    ints = [];
                else
                    % Update model and experiments with provided parameters
                    [m_i,con_i] = updateAll(m, con_i, T(TisWorker_i), opts_i.UseParams, opts_i.UseSeeds, opts_i.UseInputControls, opts_i.UseDoseControls);
                    
                    % Integrate the system
                    ncon_i = numel(con_i);
                    for i = ncon_i:-1:1
                        opts_ii = opts_i;
                        opts_ii.AbsTol = opts_ii.AbsTol{i};
                        opts_ii.UseSeeds = opts_ii.UseSeeds(:,i);
                        opts_ii.UseInputControls = opts_ii.UseInputControls{i};
                        opts_ii.UseDoseControls = opts_ii.UseDoseControls{i};
                        if j == 1
                            ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_i(:,i), opts_ii);
                        else
                            ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_i(:,i), opts_ii, obj, prevint);
                        end
                    end
                end
                
                ints = {ints};
            end
            
            clear int; % Prevents "subscripted assignment between dissimilar structures" error
            
            % Combine the workers' int structs by concatenating horizontally
            for i = NumWorkers:-1:1
                ints_temp = ints{i};
                if ~isempty(ints_temp)
                    int(:,icon_worker{i}) = ints_temp;
                end
            end
            
            clear ints; % Prevents "subscripted assignment between dissimilar structures" error
            
            % Add other fields needed
            ObjWeights = num2cell(opts.ObjWeights);
            [int.ObjWeight] = deal(ObjWeights{:});
            for i = 1:n_con
                [int(:,i).TisExperiment] = deal(TisExperiment(:,i));
                [int(:,i).T] = deal(T(TisExperiment(:,i)));
                [int(:,i).nT] = deal(sum(TisExperiment(:,i)));
            end
            
            prevint = int;
            
        end
        
        set_memoized_int(int, T, integrateFunctions{nIntegrates});
        
    end

    function varargout = objective(T)
        
        % Unnormalize
        if opts.Normalized
            T = exp(T);
        else
            % If fmincon chooses values that violate the lower bounds, force them to be equal to the lower bounds
            T(T < opts.LowerBound) = opts.LowerBound(T < opts.LowerBound);
        end
        
        fun_i = nargout;
        if fun_i == 0
            fun_i = 1;
        end
        int = integr(T, integrateFunctions_objective(fun_i,:));
        
        varargout = cell(fun_i,1);
        [varargout{:}] = objectiveFunction(obj, int);
        
    end

    function varargout = constraint(T)
        % Outputs: c,ceq,GC,GCeq
        
        % Unnormalize
        if opts.Normalized
            T = exp(T);
        else
            % If fmincon chooses values that violate the lower bounds, force them to be equal to the lower bounds
            T(T < opts.LowerBound) = opts.LowerBound(T < opts.LowerBound);
        end
        
        fun_i = nargout;
        if fun_i == 0
            fun_i = 1;
        end
        
        nconstraints = numel(integrateFunctions_constraint);
        varargout = cell(fun_i,1);
        out_i = cell(fun_i,1);
        for i = nconstraints:-1:1
            int_i = integr(T, integrateFunctions_constraint{i}(fun_i,:));
            [out_i{:}] = constraintFunctions{i}(constr_obj{i}, int_i);
            for j = 1:fun_i
                if ~isempty(out_i{j})
                    if j <= 2
                        varargout{j}(i,1) = out_i{j} - constr_vals(i);
                    else
                        varargout{j}(:,i) = out_i{j};
                    end
                end
            end
        end
        
    end

    function int = get_memoized_int(T, integrateFunction)
        
        switch func2str(integrateFunction)
            case 'integrateAllSys'
                lastintval = lastint.integrateAllSys;
                lastTval = lastT.integrateAllSys;
            case 'integrateAllSens'
                lastintval = lastint.integrateAllSens;
                lastTval = lastT.integrateAllSens;
            case 'GenerateFIMEigenvalueFunction/integrateAllCurvVecProd'
                lastintval = lastint.integrateAllCurvVecProd;
                lastTval = lastT.integrateAllCurvVecProd;
            otherwise
                if ~displayedMemoWarning
                    warning('Unrecognized integration function. Memoization of results will be disabled until GenerateObjective/integr is updated to support this integration function.');
                    displayedMemoWarning = true;
                end
                lastintval = [];
                lastTval = [];
        end
        
        if isequal(T, lastTval)
            int = lastintval;
        else
            int = [];
        end
    end

    function set_memoized_int(int, T, integrateFunction)
        switch func2str(integrateFunction)
            case 'integrateAllSys'
                lastint.integrateAllSys = int;
                lastT.integrateAllSys = T;
            case 'integrateAllSens'
                lastint.integrateAllSens = int;
                lastint.integrateAllSys = int; % Higher-order integration contains all information for lower-order integration
                lastT.integrateAllSens = T;
                lastT.integrateAllSys = T;
            case 'GenerateFIMEigenvalueFunction/integrateAllCurvVecProd'
                lastint.integrateAllCurvVecProd = int;
                lastint.integrateAllSens = int;
                lastint.integrateAllSys = int;
                lastT.integrateAllCurvVecProd = T;
                lastT.integrateAllSens = T;
                lastT.integrateAllSys = T;
        end
    end

end

function [con,obj,opts,TisWorker] = splitExperiments(con, obj, opts, i_cons, T_experiment)
% T_experiment:
%   nT-by-1 vector of experiment indices indicating which experiment fits
%   each parameter. 0 indicates model parameter fit by all experiments.
% i_cons:
%   Vector of experimental indices to be simulated by the slave function.

    %nT = numel(T_experiment);
    
    % If worker has no experiments assigned, assign a 0-valued objective
    % function to it
%     if isempty(i_cons)
%         % Clear variables to avoid wasting memory in closure
%         clear m con obj opts intOpts T_experiment
%         
%         objectivefun = @emptyobjective;
%         
%         return
%     end

    % Get information on which parameters need to be used for the worker's
    % experiments
    TisWorker = T_experiment == 0 | any(bsxfun(@eq, T_experiment, i_cons(:)'), 2);

    % Filter down input arguments to those relevant to the worker
    con = con(i_cons);
    obj = obj(:, i_cons);
    opts.UseSeeds = opts.UseSeeds(:,i_cons);
    opts.UseInputControls = opts.UseInputControls(i_cons);
    opts.UseDoseControls = opts.UseDoseControls(i_cons);
    opts.ObjWeights = opts.ObjWeights(:,i_cons);
    if iscell(opts.AbsTol)
        opts.AbsTol = opts.AbsTol(i_cons);
    end
    opts.LowerBound = opts.LowerBound(TisWorker);
    opts.UpperBound = opts.UpperBound(TisWorker);
%     intOpts.UseSeeds = intOpts.UseSeeds(:,i_cons);
%     intOpts.UseInputControls = intOpts.UseInputControls(i_cons);
%     intOpts.UseDoseControls = intOpts.UseDoseControls(i_cons);
%     intOpts.ObjWeights = intOpts.ObjWeights(:,i_cons);
%     if iscell(intOpts.AbsTol)
%         intOpts.AbsTol = intOpts.AbsTol(i_cons);
%     end
%     intOpts.LowerBound = intOpts.LowerBound(TisWorker);
%     intOpts.UpperBound = intOpts.UpperBound(TisWorker);
    
%     % Return worker-specific integration function
%     objectivefun = @objective;
%
%     function varargout = objective(T)
%     % Slave objective function
%         
%         % Get portion of T that is relevant to the worker
%         T_worker = T(TisWorker);
%         
%         % Reset answers
%         G = 0;
%         D_worker = zeros(numel(T_worker),1);
%         
%         % Unnormalize
%         if opts.Normalized
%             T_worker = exp(T_worker);
%         else
%             % If fmincon chooses values that violate the lower bounds, force them to be equal to the lower bounds
%             T_worker(T_worker < opts.LowerBound) = opts.LowerBound(T_worker < opts.LowerBound);
%         end
%         
%         % Update parameter sets
%         [m, con] = updateAll(m, con, T_worker, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
%         
%         varargout = cell(nargout,1);
%         [varargout{:}] = myobjective(m, con, obj, intOpts);
%         
% 
%         
%         % Assign back to vector sized for all workers
%         switch opts.Solver
%             case 'fmincon'
%                 D = zeros(nT, 1);
%                 D(TisWorker) = D_worker;
%             case 'lsqnonlin'
%                 D = zeros(size(D_worker,1), nT);
%                 D(:,TisWorker) = D_worker;
%         end
%     end
% 
%     function [G, D] = emptyobjective(T)
%         
%     end

end