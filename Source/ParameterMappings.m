function inds = ParameterMappings(opts)
% Use* fields must be provided in the logical format for this function to
% work

if isfield(opts.UseParams)
    inds.k = find(opts.UseParams);
    nTk = sum(opts.UseParams);
    inds.Tk = 1:nTk;
else
    nTk = 0;
    inds.k = [];
    inds.Tk = [];
end

if isfield(opts.UseSeeds)
    count = nTk+1;
    nTs = 0;
    ncon = size(opts.UseSeeds,2);
    inds.s = cell(ncon,1);
    inds.Ts = cell(ncon,1);
    for icon = 1:ncon
        inds.s{icon} = find(opts.UseSeeds(:,icon));
        inds.Ts{icon} = count : count+numel(inds.s{icon})-1;
        count = count + numel(inds.s{icon});
        nTs = nTs + numel(inds.s{icon});
    end
else
    nTs = 0;
    inds.s = {};
    inds.Ts = {};
end


if isfield(opts.UseInputControls)
    count = nTk + nTs + 1;
    nTq = 0;
    ncon = length(opts.UseInputControls);
    inds.q = cell(ncon,1);
    inds.Tq = cell(ncon,1);
    for icon = 1:ncon
        inds.q{icon} = find(opts.UseInputControls{icon});
        inds.Tq{icon} = count : count+numel(inds.q{icon})-1;
        count = count + numel(inds.q{icon});
        nTq = nTq + numel(inds.q{icon});
    end
else
    nTq = 0;
    inds.q = {};
    inds.Tq = {};
end

if isfield(opts.UseDoseControls)
    count = nTk + nTs + nTq + 1;
    nTh = 0;
    ncon = length(opts.UseDoseControls);
    inds.h = cell(ncon,1);
    inds.Th = cell(ncon,1);
    for icon = 1:ncon
        inds.h{icon} = find(opts.UseDoseControls{icon});
        inds.Th{icon} = count : count+numel(inds.h{icon})-1;
        count = count + numel(inds.h{icon});
        nTh = nTh + numel(inds.h{icon});
    end
else
    nTh = 0;
    inds.h = {};
    inds.Th = {};
end