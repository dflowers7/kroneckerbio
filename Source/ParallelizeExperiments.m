function varargout = ParallelizeExperiments(fun, m, con, obs, opts)

persistent m_i con_i obs_i opts_i k s q h initialized icon_per_worker

%% Initialize variables
% This block only runs if it is the first time it has been run since the
% function was cleared.

nargo = nargout;

if isempty(initialized)
    initialized = false;
end

NumWorkers = initializeParallelPool();
if ~initialized
    fprintf('Initializing parallel variables...\n')
    if NumWorkers > 1
        [m_i, con_i, obs_i, opts_i, k, s, q, h] = deal(Composite());
    else
        varargout = cell(nargo,1);
        [varargout{:}] = fun(m, con, obs, opts);
        return
    end
    
    % Assign experiments to workers
    n_con = numel(con);
    ncon_per_worker = floor(n_con/NumWorkers);
    ncon_per_worker = repmat(ncon_per_worker, NumWorkers, 1);
    nleft = n_con - sum(ncon_per_worker);
    ncon_per_worker(1:nleft) = ncon_per_worker(1:nleft)+1;
    icon_per_worker = mat2cell((1:n_con)', ncon_per_worker, 1);
    
    % Fix options (may need to do this better later)
    switch func2str(fun)
        case 'SimulateSystem'
            derorder = 0;
        case 'SimulateSensitivity'
            derorder = 1;
        case 'SimulateCurvature'
            derorder = 2;
        otherwise
            derorder = 0;
    end
    opts = FixSimulationOpts(m, con, obs, opts, derorder);
    
    % Transfer each experiment's variables into the workers
    for i = 1:NumWorkers
        m_i{i} = m;
        icon_w = icon_per_worker{i};
        con_i{i} = con(icon_w);
        obs_i{i} = obs(:,icon_w);
        
        opts_ii = opts;
        opts_table = {
            'AbsTol'                        @elementassign
            'UseSeeds'                      @columnassign
            'UseInputControls'              @elementassign
            'UseDoseControls'               @elementassign
            'AdjointOutputSensitivities'    @columnassign
            };
        for j = 1:size(opts_table,1)
            fieldj = opts_table{j,1};
            funj = opts_table{j,2};
            opts_ii = funj(opts_ii, fieldj, icon_w);
        end
        opts_i{i} = opts_ii;
    end
    
    initialized = true;
end

%% Run function

% qq = tic;
for i = 1:NumWorkers
    k{i} = m.k;
    con_iw = con(icon_per_worker{i});
    s{i} = [con_iw.s];
    q{i} = {con_iw.q};
    h{i} = {con_iw.h};
end
% fprintf('Updating parameters took %0.3g seconds.\n', toc(qq))

spmd
    % Update parameters
    m_i = m_i.Update(k);
    for i = 1:numel(con_i)
        con_i(i) = con_i(i).Update(s(:,i), q{i}, h{i});
    end
    
    % Run function
    varargout = cell(nargo,1);
    [varargout{:}] = fun(m_i, con_i, obs_i, opts_i);
end

%% Concatenate outputs

% qq = tic;
NumWorkers = numel(varargout);
temp = cell(NumWorkers, nargo);
for i = 1:NumWorkers
    for j = 1:nargo
        varargout_i = varargout{i};
        temp{i,j} = varargout_i{j};
    end
end
varargout = cell(nargo,1);
for i = 1:nargo
    varargout{i} = [temp{:,i}];
end
% fprintf('Assembling outputs took %0.3g seconds.\n', toc(qq))

    function NumWorkers = initializeParallelPool()
        % Check for parallel toolbox if parallel optimization is specified
        if isempty(ver('distcomp'))
            warning('KroneckerBio:FitObjective:ParallelToolboxMissing', ...
                ['Parallel experiments requires the Parallel '...
                'Computing toolbox. Reverting to serial evaluation.'])
            NumWorkers = 1;
            return
        end
        
        % Determine number of workers in the parallel pool
        p = gcp();
        if isempty(p)
            warning('KroneckerBio:FitObjective:NoParallelPool', ...
                ['No parallel pool could be initialized. '...
                'Reverting to serial optimization.']);
            opts.ParallelizeExperiments = false;
            NumWorkers = 1;
        else
            NumWorkers = p.NumWorkers;
        end
    end


end

%% Functions for assigning options to subsets of experiments

function opts = elementassign(opts, field, icon_w)

if isfield(opts, field)
    opts.(field) = opts.(field)(icon_w);
end

end

function opts = columnassign(opts, field, icon_w)

if isfield(opts, field)
    opts.(field) = opts.(field)(:,icon_w);
end

end