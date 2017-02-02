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

% Unify the objective and all the constraints' objectives into one struct
nIntegrateFunctions_constraint = sum(cellfun(@numel, integrateFunctions_constraint));
nIntegrateFunctions = numel(integrateFunctions_objective) + nIntegrateFunctions_constraint;
int_all = struct('index', cell(nIntegrateFunctions,1), ...
    'isobjective', cell(nIntegrateFunctions,1), 'integratefun', cell(nIntegrateFunctions,1));
count = 1;
for ii = 1:size(integrateFunctions_objective,1)
    for jj = 1:size(integrateFunctions_objective,2)
        int_all(count).index = 1;
        int_all(count).isobjective = true;
        int_all(count).integratefun = func2str(integrateFunctions_objective{ii,jj});
        count = count + 1;
    end
end
for kk = 1:nconstraints
    intfuns_kk = integrateFunctions_constraint{kk};
    for ii = 1:size(intfuns_kk,1)
        for jj = 1:size(intfuns_kk,2)
            if isempty(intfuns_kk{ii,jj})
                continue % Empty placeholder. Skip to next
            else
                funcstr = func2str(intfuns_kk{ii,jj});
            end
            if ~any(ismember([int_all(1:count-1).index]', kk) & ismember([int_all(1:count-1).isobjective]', false) & ismember({int_all(1:count-1).integratefun}', funcstr))  
            % Only add unique entries
                int_all(count).index = kk;
                int_all(count).isobjective = false;
                int_all(count).integratefun = funcstr;
                count = count + 1;
            end
        end
    end
end
int_all(count:end) = [];

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
        NumWorkers = 1;
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
[icon_worker,cons,objs,optss,TisWorker,constr_objs] = deal(cell(NumWorkers, 1));
endi = 0;
for ii = 1:NumWorkers
    starti = endi+1;
    endi = sum(numPerWorker(1:ii));
    icon_worker{ii} = (starti:endi)';
    [cons{ii}, objs{ii}, optss{ii}, TisWorker{ii}, constr_objs{ii}] = splitExperiments(con, obj, opts, icon_worker{ii}, T_experiment, constr_obj);
end

% Distribute experiments to workers
if opts.ParallelizeExperiments
    spmd
        cons = codistributed(cons);
        objs = codistributed(objs);
        optss = codistributed(optss);
        TisWorker = codistributed(TisWorker);
        icon_workers = codistributed(icon_worker);
        constr_objs = codistributed(constr_objs);
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

    function int_return = integr(T, integrateFunctions, isobjective_return, index_return)
        
        nIntegrates = sum(~cellfun(@isempty,integrateFunctions));
        
        prevint = [];
        prevint_filtered = [];
        prev_obj_all_cell = [];% Initialize to prevent "Analyzing and transferring files to the workers..." message
        prev_obj_all_cell_filtered = [];
        prev_index = [];
        prev_isobjective = [];
        
        for j = 1:nIntegrates
            
            integrateFunStr_j = func2str(integrateFunctions{j});
            is_int_all = strcmp({int_all.integratefun}, integrateFunStr_j);
            isobjective = vertcat(int_all(is_int_all).isobjective);
            index = vertcat(int_all(is_int_all).index);
            
            int_temp = get_memoized_int(T, integrateFunctions{j});
            if ~isempty(int_temp)
                int_cell = int_temp;
                prevint = int_cell;
                prev_obj_all_cell = cell(numel(index),1);
                prev_obj_all_cell(isobjective) = {obj};
                prev_obj_all_cell(~isobjective) = constr_obj(index(~isobjective));
                prev_index = index;
                prev_isobjective = isobjective;
                continue
            end
            
            % Filter the previous integration struct to those ints
            % needed for this round of integration
            if ~isempty(prevint)
                intstokeep = zeros(numel(index),1);
                for q = 1:numel(index)
                    intstokeep(q) = find(index(q) == prev_index & isobjective(q) == prev_isobjective);
                end
                assert(numel(intstokeep) == numel(index), 'KroneckerBio:GenerateObjective:PrevIntNotFound', ...
                    'Some of the ints needed for the next round of integration were not found. This is a bug.')
                prevint_filtered = prevint(intstokeep);
                prev_obj_all_cell_filtered = prev_obj_all_cell(intstokeep);
            end
            
            if opts.ParallelizeExperiments
                spmd
                    % Get local experiments
                    con_i = getLocalPart(cons);
                    con_i = con_i{1};
                    obj_i = getLocalPart(objs);
                    obj_i = obj_i{1};
                    constr_obj_i = getLocalPart(constr_objs);
                    constr_obj_i = constr_obj_i{1};
                    opts_i = getLocalPart(optss);
                    opts_i = opts_i{1};
                    TisWorker_i = getLocalPart(TisWorker);
                    TisWorker_i = TisWorker_i{1};
                    %icon_worker_i = getLocalPart(icon_workers);
                    
                    % Filter down to the objectives matching the requested
                    % integrate function
                    obj_all_i_cell = cell(numel(index),1);
                    obj_all_i_cell(isobjective) = {obj_i};
                    obj_all_i_cell(~isobjective) = constr_obj_i(index(~isobjective));
                    obj_all_i = vertcat(obj_all_i_cell{:});
                    
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
                                ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_all_i(:,i), opts_ii);
                            else
                                ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_all_i(:,i), opts_ii, prev_obj_all_cell_filtered, prevint_filtered);
                            end
                        end
                    end
                end
            else
                % Get local experiments
                con_i = cons{1};
                obj_i = objs{1};
                constr_obj_i = constr_objs{1};
                opts_i = optss{1};
                TisWorker_i = TisWorker{1};
                
                % Filter down to the objs matching the requested
                % integrate function
                obj_all_i_cell = cell(numel(index),1);
                obj_all_i_cell(isobjective) = {obj_i};
                obj_all_i_cell(~isobjective) = constr_obj_i(index(~isobjective));
                obj_all_i = vertcat(obj_all_i_cell{:});
                
                if isempty(con_i)
                    ints = [];
                else
                    % Update model and experiments with provided parameters
                    [m_i,con_i] = updateAll(m, con_i, T(TisWorker_i), opts_i.UseParams, opts_i.UseSeeds, opts_i.UseInputControls, opts_i.UseDoseControls);
                    
                    % Filter down previous int struct to those applicable to this integration
                    if ~isempty(prevint)
                    end
                    
                    % Integrate the system
                    ncon_i = numel(con_i);
                    for i = ncon_i:-1:1
                        opts_ii = opts_i;
                        opts_ii.AbsTol = opts_ii.AbsTol{i};
                        opts_ii.UseSeeds = opts_ii.UseSeeds(:,i);
                        opts_ii.UseInputControls = opts_ii.UseInputControls{i};
                        opts_ii.UseDoseControls = opts_ii.UseDoseControls{i};
                        if j == 1
                            ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_all_i(:,i), opts_ii);
                        else
                            ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_all_i(:,i), opts_ii, prev_obj_all_cell_filtered, prevint_filtered);
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
            
            % Separate int array into cells representing the objective and
            % different constraints
            obj_all_cell = cell(numel(index),1);
            obj_all_cell(isobjective) = {obj};
            obj_all_cell(~isobjective) = constr_obj(index(~isobjective));
            int_cell = mat2cell(int, cellfun(@(x)size(x,1), obj_all_cell), n_con);
            
            % Add other fields needed
            ObjWeights = num2cell(opts.ObjWeights);
            if any(isobjective)
                [int_cell{isobjective}.ObjWeight] = deal(ObjWeights{:});
            end
            for i = 1:n_con
                for k = 1:size(int_cell,1)
                    [int_cell{k}(:,i).TisExperiment] = deal(TisExperiment(:,i));
                    [int_cell{k}(:,i).T] = deal(T(TisExperiment(:,i)));
                    [int_cell{k}(:,i).nT] = deal(sum(TisExperiment(:,i)));
                    [int_cell{k}.ObjWeight] = deal(1); % Use weights of 1 for constraints. No support for varying ObjWeights for constraints at this point.
                end
            end
            
            prevint = int_cell;
            prev_obj_all_cell = obj_all_cell;
            prev_index = index;
            prev_isobjective = isobjective;
            
            set_memoized_int(int_cell, T, integrateFunctions{j});
            
        end
        
        % Return only the requested objective or constraint's int struct
        if isobjective_return
            int_return = int_cell{isobjective};
        else
            int_return = int_cell{~isobjective & index == index_return};
        end
        
    end

    function varargout = objective(T)
        
        % Unnormalize
        if opts.Normalized
            T = exp(T);
        else
            % If fmincon chooses values that violate the lower bounds, force them to be equal to the lower bounds
            T(T < opts.LowerBound) = opts.LowerBound(T < opts.LowerBound);
        end
        
        % Get the int struct for all the obj structs (objective and all
        % constraints)
        fun_i = nargout;
        if fun_i == 0
            fun_i = 1;
        end
        isobjective = true;
        int = integr(T, integrateFunctions_objective(fun_i,:), isobjective, 1);
        
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
        isobjective = false;
        for i = nconstraints:-1:1
            index = i;
            int_i = integr(T, integrateFunctions_constraint{i}(fun_i,:), isobjective, index);
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

function [con,obj,opts,TisWorker,constr_obj] = splitExperiments(con, obj, opts, i_cons, T_experiment, constr_obj)
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
    for i = 1:numel(constr_obj)
        constr_obj{i} = constr_obj{i}(:,i_cons);
    end
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