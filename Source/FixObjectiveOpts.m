function opts = FixObjectiveOpts(m, con, obj, opts, derorder)

% Fix simulation opts
opts = FixSimulationOpts(m, con, obj, opts, derorder);

% Fix additional objective specific options
defaultOpts.ObjWeights       = ones(size(obj));
defaultOpts.Normalized       = true;
defaultOpts.UseAdjoint       = true;

opts = mergestruct(defaultOpts, opts);

end