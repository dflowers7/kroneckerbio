function [objfun, constrfun, hessianapproxfun, outfun, updateoptsfun] = GenerateObjective(m, con, obj, opts, integrateFunctions_objective, objectiveFunction, integrateFunctions_constraint, constraintFunctions, funopts, other_outfun)
% Input arguments:
%   m
%   con
%   obj
%   opts
%   integrateFunctions_objective [ cell array ]
%       The cell array is of size nargout-by-nintegratefunctions. When the
%       objective function is called with i output arguments, the
%       integration functions in row i are performed one at a time, from
%       left to right, with the data from the previous integration of all
%       experiments accessible to the next integration.
%   objectiveFunction [ function handle ]
%       The function(obj,int) takes the obj struct and the int struct returned by
%       integration and returns a scalar value representing the objective
%       value and a vector representing the objective gradient. The default
%       is the sum of the individual objective values/gradients.
%   integrateFunctions_constraint [ cell vector of cell arrays ]
%       Outer cell element i contains a cell array of integration functions of
%       size nargout-by-nintegratefunctions for constraint i. The inner
%       cell arrays are structured the same way as the
%       integrateFunctions_objective option, described above. The default
%       is system integration for one or two output arguments and
%       sensitivity integration for three or four output arguments.
%   constraintFunctions [ cell vector of function handles ]
%       Element i contains a function handle that reduces the values
%       returned by the different objective functions into one value. The
%       default is to sum the objective values.
%   funopts
%       .ScaleConstraints
%           Scales constraints by dividing by the constraint value such that all
%           constraints will be of the form 
%               f(T)/abs(constraint_value) <= sign(constraint_value) 
%           instead of 
%               f(T) <= constraint_value.
%       .HessianMaximumConditionNumber [ 1000 ]
%           Sets the maximum condition number allowed for the approximated
%           Hessian. Smaller values put a stricter limit on the step size
%           allowed at each optimization iteration.
%       .ApproximateSecondOrderHessianTerm [ true ]
%           If set to true, uses structured BFGS updates to approximate the
%           second-order portions of the Hessian. If set to false, only the
%           first-order least squares terms are used.

if nargin < 10
    other_outfun = [];
    if nargin < 9
        funopts = [];
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
    end
end

if isempty(other_outfun)
    other_outfun = @(x, optimValues, state) false;
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
        case 'sqp'
            objectiveFunction = @sumobjectives_fmincon;
    end
end

nconstraintfuns = numel(opts.ConstraintObj);

if isempty(integrateFunctions_constraint)
    integrateFunctions_constraint = repmat({{@integrateAllSys;@integrateAllSys;@integrateAllSens;@integrateAllSens}}, nconstraintfuns, 1);
end

if isempty(constraintFunctions)
    switch opts.Solver
        case 'fmincon'
            constraintFunctions = repmat({@sumobjectives_fmincon_constraint}, nconstraintfuns, 1);
        case 'lsqnonlin'
            constraintFunctions = cell(0,1);
        case 'sqp'
            constraintFunctions = repmat({@sumobjectives_fmincon_constraint}, nconstraintfuns, 1);
    end
end

% Set default options for objective and constraint functions
if isempty(funopts)
    funopts = struct;
end
defaultFunOpts.ScaleConstraints = false;
defaultFunOpts.HessianMaximumConditionNumber = 1000;
funopts = mergestruct(defaultFunOpts, funopts);

% Extract and check funopts fields
ScaleConstraints = funopts.ScaleConstraints;
assert(islogical(ScaleConstraints) && isscalar(ScaleConstraints), ...
    'KroneckerBio:GenerateObjective:ScaleConstraints', ...
    'ScaleConstraints must be a logical scalar.')

assert(~strcmp(opts.Solver,'lsqnonlin') || nconstraintfuns == 0, ...
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
intfuns_kk = integrateFunctions_objective;
kk = 1;
for ii = 1:size(intfuns_kk,1)
    for jj = 1:size(intfuns_kk,2)
        if isempty(intfuns_kk{ii,jj})
            continue
        else
            funcstr = func2str(intfuns_kk{ii,jj});
        end
        if ~any(ismember([int_all(1:count-1).index]', kk) & ismember([int_all(1:count-1).isobjective]', true) & ismember({int_all(1:count-1).integratefun}', funcstr))
            % Only add unique entries
            int_all(count).index = kk;
            int_all(count).isobjective = true;
            int_all(count).integratefun = funcstr;
            count = count + 1;
        end
    end
end
for kk = 1:nconstraintfuns
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

initializeParallelPool()

    function initializeParallelPool()
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
    icon_worker_ii = (starti:endi)';
    icon_worker(ii) = {icon_worker_ii};
    ans = cell(5,1);
    [ans{:}] = splitExperiments(con, obj, opts, icon_worker_ii, T_experiment, constr_obj);
    cons(ii) = ans(1);
    objs(ii) = ans(2);
    optss(ii) = ans(3);
    TisWorker(ii) = ans(4);
    constr_objs(ii) = ans(5);
    %[cons{ii}, objs{ii}, optss{ii}, TisWorker{ii}, constr_objs{ii}] = splitExperiments(con, obj, opts, icon_worker{ii}, T_experiment, constr_obj);
end

% Distribute experiments to workers
% if opts.ParallelizeExperiments
%     cons = distributed(cons);
% end
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

% Initialize Hessian approximation variables
gs_last_hess_possible = [];
errs_last_hess_possible = [];
derrdTs_last_hess_possible = [];
As_last_hess_possible = [];
T_last_hess = [];
gs_last_hess = [];
errs_last_hess = [];
derrdTs_last_hess = [];
As_last_hess = [];
max_hess_condno = funopts.HessianMaximumConditionNumber;
isExactHessian = [];
hasLeastSquares = [];

updateHessianVarsInHessFun = strcmp(opts.Solver, 'sqp');

objfun = @objective;
if hasconstraint
    constrfun = @constraint;
else
    % Set to empty so that fmincon knows there is no constraint
    constrfun = [];
end
hessianapproxfun = @hessian_approximation;
outfun = @outfun_;
updateoptsfun = @updateOpts;

    function int_return = integr(T, integrateFunctions, isobjective_return, index_return)
        
        timeIntegration = false; % Set to true to time what fraction of the time is spent integrating versus overhead
        
        if timeIntegration
            tottimer = tic;
        end
        
        nIntegrates = sum(~cellfun(@isempty,integrateFunctions));
        
        prevint = [];
        prevint_filtered = [];
        prev_obj_all_cell = [];% Initialize to prevent "Analyzing and transferring files to the workers..." message
        prev_obj_all_cell_filtered = [];
        prev_index = [];
        prev_isobjective = [];
        
        for j = 1:nIntegrates
            
            if timeIntegration
                timer_initint = tic;
            end
            
            integrateFunStr_j = func2str(integrateFunctions{j});
            is_int_all = strcmp({int_all.integratefun}, integrateFunStr_j);
            isobjective = vertcat(int_all(is_int_all).isobjective);
            index = vertcat(int_all(is_int_all).index);
            
            int_cell = get_memoized_int(T, integrateFunctions{j});
%             if ~isempty(int_temp)
%                 
%                 int_cell = int_temp;
%                 prevint = int_cell;
%                 prev_obj_all_cell = cell(numel(index),1);
%                 prev_obj_all_cell(isobjective) = {obj};
%                 prev_obj_all_cell(~isobjective) = constr_obj(index(~isobjective));
%                 prev_index = index;
%                 prev_isobjective = isobjective;
%                 if timeIntegration
%                     inttime_workers = repmat({0},NumWorkers,1);
%                 end
%                 
            if isempty(int_cell)
                
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
                
                if timeIntegration
                    t_initint = toc(timer_initint);
                    fprintf('Initialization of integration %g took %g seconds\n', j, t_initint)
                end
                
                if opts.ParallelizeExperiments
                    if timeIntegration
                        timer_spmd = tic;
                    end
                    spmd
                        if timeIntegration
                            t_spmd = toc(timer_spmd);
                            fprintf('Starting spmd block took %g seconds\n', t_spmd)
                        end
                        if timeIntegration
                            inttime_workers = 0;
                        end
                        
                        % Get local experiments
                        if timeIntegration
                            timer_getlocalparts = tic;
                        end
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
                        if timeIntegration
                            t_getlocalparts = toc(timer_getlocalparts);
                            fprintf('Took %g seconds for worker to get local parts\n', t_getlocalparts)
                        end
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
                                if timeIntegration
                                    fprintf('Integrating experiment %s\n', con_i(i).Name)
                                    inttimer = tic;
                                end
                                if j == 1
                                    ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_all_i(:,i), opts_ii);
                                else
                                    ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_all_i(:,i), opts_ii, prev_obj_all_cell_filtered, prevint_filtered);
                                end
                                if timeIntegration
                                    tcon = toc(inttimer);
                                    fprintf('Experiment %s took %g seconds to integrate\n', con_i(i).Name, tcon)
                                    inttime_workers = inttime_workers + tcon;
                                end
                            end
                        end
                    end
                else
                    if timeIntegration
                        inttime_workers = 0;
                    end
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
                        
                        % Integrate the system
                        ncon_i = numel(con_i);
                        for i = ncon_i:-1:1
                            opts_ii = opts_i;
                            opts_ii.AbsTol = opts_ii.AbsTol{i};
                            opts_ii.UseSeeds = opts_ii.UseSeeds(:,i);
                            opts_ii.UseInputControls = opts_ii.UseInputControls{i};
                            opts_ii.UseDoseControls = opts_ii.UseDoseControls{i};
                            if timeIntegration
                                fprintf('Integrating experiment %s\n', con_i(i).Name)
                                inttimer = tic;
                            end
                            if j == 1
                                ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_all_i(:,i), opts_ii);
                            else
                                ints(:,i) = integrateFunctions{j}(m_i, con_i(i), obj_all_i(:,i), opts_ii, prev_obj_all_cell_filtered, prevint_filtered);
                            end
                            if timeIntegration
                                tcon = toc(inttimer);
                                fprintf('Experiment %s took %g seconds to integrate\n', con_i(i).Name, tcon)
                                inttime_workers = inttime_workers + toc(inttimer);
                            end
                        end
                        
                    end
                    
                    if timeIntegration
                        inttime_workers = {inttime_workers};
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
                
            end
            
            % Add other fields needed
            ObjWeights = num2cell(opts.ObjWeights);
            isobjective_i = find(isobjective);
            isconstraint_i = find(~isobjective);
            % Special fields for objectives
            for oi = isobjective_i(:).'
                [int_cell{oi}.ObjWeight] = deal(ObjWeights{:});
            end
            % Special fields for constraints
            for ci = isconstraint_i(:).'
                [int_cell{ci}.ObjWeight] = deal(1); % Use weights of 1 for constraints. No support for varying ObjWeights for constraints at this point.
            end
            % General fields needed for both objectives and constraints
            for i = 1:n_con
                for k = 1:size(int_cell,1)
                    [int_cell{k}(:,i).TisExperiment] = deal(TisExperiment(:,i));
                    [int_cell{k}(:,i).T] = deal(T(TisExperiment(:,i)));
                    [int_cell{k}(:,i).nT] = deal(sum(TisExperiment(:,i)));
                end
            end
            
            prevint = int_cell;
            prev_obj_all_cell = cell(numel(index),1);
            prev_obj_all_cell(isobjective) = {obj};
            prev_obj_all_cell(~isobjective) = constr_obj(index(~isobjective));
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
        
        if timeIntegration
            tottime = toc(tottimer);
            inttime = 0;
            for i = 1:NumWorkers
                inttime = inttime + inttime_workers{i};
            end
            fprintf('Integration functions took %g percent of the total time\n', inttime./tottime*100)
        end
        
    end

    function [G,D,err,derrdT,Happrox] = objective(T)
        
        % Unnormalize
        if opts.Normalized
            T = exp(T);
        else
            % If fmincon chooses values that violate the lower bounds, force them to be equal to the lower bounds
            T(T < opts.LowerBound) = opts.LowerBound(T < opts.LowerBound);
        end
        
        % Get the int struct for all the obj structs (objective and all
        % constraints)
        nargout_ = nargout;
        if nargout_ == 0
            nargout_ = 1;
        end
        isobjective = true;
        integrfun_i = min(nargout_, 2);
        int = integr(T, integrateFunctions_objective(integrfun_i,:), isobjective, 1);
        
        varargout = cell(nargout_,1);
        [varargout{:}] = objectiveFunction(obj, int);
        G = varargout{1};
        if nargout_ > 1
            D = varargout{2};
            if nargout_ > 2
                err = varargout(3);
                if nargout_ > 3
                    derrdT = varargout(4);
                    if nargout_ > 4
                        Happrox = varargout(5);
                    end
                end
            end
        end
%         if nargout_ > 4
%             varargout{5} = varargout(5);
%         end
        
    end

    function [c,ceq,GC,GCeq,err,derrdT,Happrox] = constraint(T)
        % Outputs: c,ceq,GC,GCeq
        
        % Unnormalize
        if opts.Normalized
            T = exp(T);
        else
            % If fmincon chooses values that violate the lower bounds, force them to be equal to the lower bounds
            T(T < opts.LowerBound) = opts.LowerBound(T < opts.LowerBound);
        end
        
        nargout_ = nargout;
        if nargout_ == 0
            nargout_ = 1;
        end
        
        nconstraintfuns = numel(integrateFunctions_constraint);
        
        % Initialize output cell arrays
        [c,ceq,GC,GCeq,err,derrdT,Happrox] = deal(cell(nconstraintfuns,1));
        
%         varargout = cell(nargout_,1);
%         for j = 1:nargout_
%             varargout{j} = cell(nconstraintfuns,1);
%         end
        
        % Calculate constraint function value for each constraint function
%         out_i = cell(nargout_,1);
        isobjective = false;
        integrfun_i = min(nargout_, 4);
        for i = 1:nconstraintfuns
            index = i;
            int_i = integr(T, integrateFunctions_constraint{i}(integrfun_i,:), isobjective, index);
            if nargout_ < 3
                [c{i},ceq{i}] = constraintFunctions{i}(constr_obj{i}, int_i);
            else
                [c{i},ceq{i},GC{i},GCeq{i},err{i},derrdT{i},Happrox{i}] = constraintFunctions{i}(constr_obj{i}, int_i);
            end
%             for j = 1:nargout_
%                 if ~isempty(out_i{j})
%                     varargout{j}{i} = out_i{j};
%                 end
%             end
        end
        
        % Concatenate constraint function values across constraint
        % functions
        c = vertcat(c{:});
        ceq = vertcat(ceq{:});
        if nargout_ > 2
            GC = [GC{:}];
            GCeq = [GCeq{:}];
        end
%         for j = 1:nargout_
%             if j <= 2
%                 varargout{j} = vertcat(varargout{j}{:});
%             elseif j <= 4
%                 varargout{j} = [varargout{j}{:}];
%             else
%                 varargout{j} = vertcat(varargout{j}{:});
%             end
%         end
        
        % Subtract constraint values from evaluated values, such that
        % all constraints are considered violated when positive
        if ~isempty(c)
            c = c - constr_vals;
        end
        if ~isempty(ceq)
            ceq = ceq - constr_vals;
        end
%         for j = 1:min(nargout_,2)
%             if ~isempty(varargout{j})
%                 varargout{j} = varargout{j} - constr_vals;
%             end
%         end
        
        % Scale by constraint values, if requested
        if ScaleConstraints
            if ~isempty(c)
                c = c./abs(constr_vals);
            end
            if ~isempty(ceq)
                ceq = ceq./abs(constr_vals);
            end
            if nargout_ > 2
                if ~isempty(GC)
                    GC = bsxfun(@rdivide, GC, abs(constr_vals(:).'));
                end
                if ~isempty(GCeq)
                    GCeq = bsxfun(@rdivide, GCeq, abs(constr_vals(:).'));
                end
                if nargout_ > 4
                    for k = 1:nconstraintfuns
                        if ~isempty(err{k})
                            err{k} = err{k}./sqrt(abs(constr_vals(k)));
                        end
                        if ~isempty(derrdT{k})
                            derrdT{k} = derrdT{k}./sqrt(abs(constr_vals(k)));
                        end
                    end
                end
            end
            
%             for j = 1:nargout_
%                 if ~isempty(varargout{j})
%                     if j <= 2
%                         varargout{j} = varargout{j}./abs(constr_vals);
%                     elseif j <= 4
%                         varargout{j} = bsxfun(@rdivide, varargout{j}, abs(constr_vals(:).'));
%                     else
%                         for k = 1:numel(varargout{j})
%                             varargout{j}{k} = varargout{j}{k}./sqrt(abs(constr_vals(k)));
%                         end
%                     end
%                 end
%             end
        end
        
    end

    function H = hessian_approximation(T, lambda)
        
        % If is a struct, extract the constraint
        if isstruct(lambda)
            lc_constr = lambda.ineqnonlin(:);
        else
            lc_constr = lambda(:);
        end
        
        % Calculate new gradients and get Lagrange multipliers
        [~, g_obj, err_obj, derrdT_obj, Happrox_obj] = objective(T);
        g_obj = {g_obj};
        lc = 1;
        if hasconstraint
            [~, ~, g_constr, ~, err_constr, derrdT_constr, Happrox_constr] = constraint(T);
            g_constr = mat2cell(g_constr, nT, ones(size(g_constr,2),1));
            g_constr = g_constr(:);
            lc = [lc; lc_constr];
        else
            g_constr = {};
            Happrox_constr = {};
            err_constr = {};
            derrdT_constr = {};
        end
        gs = [g_obj; g_constr];
        errs = [err_obj; err_constr];
        derrdTs = [derrdT_obj; derrdT_constr];
        Happroxs = [Happrox_obj; Happrox_constr];

        % Determine which objective functions have least squares terms
        % involved
        if isempty(hasLeastSquares)
            hasLeastSquares = cellfun(@(err)~isempty(err), errs);
        end
        
        % Determine which objective functions have exact Hessian
        % approximations (such as for priors, where the calculation of the
        % exact Hessian is fast)
        if isempty(isExactHessian)
            isExactHessian = false(numel(gs), 1);
            for i = 1:numel(gs)
                if i == 1
                    isObjZero = strcmp({obj.Type}, 'Objective.Data.Zero');
                    isExactHessian(i) = all([obj.approximateIsExactHessian] | isObjZero);
                else
                    isObjZero = strcmp({constr_obj.Type}, 'Objective.Data.Zero');
                    isExactHessian(i) = all([constr_obj{i-1}.approximateIsExactHessian] | isObjZero);
                end
            end
        end
        
        % Initialize last iteration terms
        if isempty(T_last_hess)
            gs_last_hess_possible = gs;
            errs_last_hess_possible = errs;
            derrdTs_last_hess_possible = derrdTs;
            As_last_hess_possible = repmat({1e-3*eye(nT)}, numel(gs), 1);
            As_last_hess_possible(isExactHessian) = {zeros(nT)};
            
            H = zeros(nT);
            for i = 1:numel(Happroxs)
                H = H + lc(i)*(Happroxs{i} + As_last_hess_possible{i});
            end
            
            if updateHessianVarsInHessFun
                T_last_hess = T;
                gs_last_hess = gs_last_hess_possible;
                errs_last_hess = errs_last_hess_possible;
                derrdTs_last_hess = derrdTs_last_hess_possible;
                As_last_hess = As_last_hess_possible;
            end
            
            return
        end
        
        % Calculate approximate Hessian for each objective term
        s = T - T_last_hess;
        [As,Bs] = deal(cell(numel(Happroxs),1));
        for i = 1:numel(Happroxs)
            
            if isExactHessian(i)
                y_hash = zeros(nT,1);
                y_s = Happroxs{i}*s;
            elseif hasLeastSquares(i)
                % y_hash is the Hessian-s product approximation for the
                % unknown term of the Hessian
                y_hash = 2*(derrdTs{i} - derrdTs_last_hess{i}).'*errs{i};
                % Modify y_hash to guarantee that y_hash.'*s is greater
                % than 0
                tk = 1e-3 + max(-y_hash.'*s/norm(s).^2,0);
                y_hash = y_hash + tk*eye(nT)*s;
                % y_s is the Hessian-s product approximation for the whole
                % Hessian
                y_s = y_hash + Happroxs{i}*s;
            else
                % If not a least squares objective function, the entirety of the Hessian is
                % unknown, so y_hash contains everything, and y_s = y_hash.
                % In this case we're using the standard BFGS update.
                y_hash = gs{i} - gs_last_hess{i};
                % Modify y_hash to guarantee that y_hash.'*s is greater
                % than 0
                tk = 1e-3 + max(-y_hash.'*s/norm(s).^2,0);
                y_hash = y_hash + tk*eye(nT)*s;
                y_s = y_hash;
            end
            
            % Threshold negative eigenvalues to zero before calculating the update
            B_s = As_last_hess{i} + Happroxs{i};
            [eigvec,eigval] = eig(B_s, 'vector');
            eigval = abs(eigval);
            eigvec_scaled = eigvec*diag(sqrt(eigval));
            B_s = eigvec_scaled*eigvec_scaled.';
            
            if norm(s) < 1e-30
                del = zeros(nT);
            else
                temp = (s.'*eigvec)*diag(sqrt(eigval));
                v = y_s + ((y_s.'*s)/(temp*temp.')).^0.5*B_s*s;

                temp = y_hash - As_last_hess{i}*s;
                temp2 = temp*v.';
                del = (temp2 + temp2.')/(v.'*s) + (temp.'*s)*(v*v.')/(v.'*s).^2;
            end
            
            As{i} = As_last_hess{i} + del;
            
            Bs{i} = As{i} + Happroxs{i};
        end
        
        % Combine individual objective function Hessians
        H = zeros(nT);
        for i = 1:numel(Bs)
            H = H + lc(i)*Bs{i};
        end
        
        % Correct eigenvalues to be above minimum threshold defined by
        % minimum condition number
        [eigvecs,eigvals] = eig(H, 'vector');
        eigvals = abs(real(eigvals));
        min_hess_eigval = max(eigvals)/max_hess_condno;
        eigvals(eigvals < min_hess_eigval) = min_hess_eigval;
        eigvecs_scaled = eigvecs*diag(sqrt(eigvals));
        H = eigvecs_scaled*eigvecs_scaled.';
        %H = (H+H.')/2; % Ensure symmetric enough to avoid warnings
        
        if any(imag(H(:)) > 0)
            warning('Complex values of H detected. Something is wrong.')
        end
        
        % Save variables
        gs_last_hess_possible = gs;
        errs_last_hess_possible = errs;
        derrdTs_last_hess_possible = derrdTs;
        As_last_hess_possible = As;
        
        if updateHessianVarsInHessFun
            T_last_hess = T;
            gs_last_hess = gs_last_hess_possible;
            errs_last_hess = errs_last_hess_possible;
            derrdTs_last_hess = derrdTs_last_hess_possible;
            As_last_hess = As_last_hess_possible;
        end
        
    end

    function stop = outfun_(x, optimValues, state)
        % Output function updates "last" values used for Hessian
        % approximation via BFGS. This is necessary because we don't want
        % to update the values except at the conclusion of an iteration
        switch state
            case 'iter'
                T_last_hess = T;
                gs_last_hess = gs_last_hess_possible;
                errs_last_hess = errs_last_hess_possible;
                derrdTs_last_hess = derrdTs_last_hess_possible;
                As_last_hess = As_last_hess_possible;
        end
        stop = other_outfun(x, optimValues, state);
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

    function updateOpts(newOpts)
        opts = newOpts;
        initializeParallelPool()
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