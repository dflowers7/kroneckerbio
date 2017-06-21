function AbsTolY = fixAbsTolY(AbsTolY, ny, UseParams, UseSeeds, UseInputControls, UseDoseControls, derorder)

if nargin > 2
    nTk = sum(UseParams);
    nTs = sum(sum(UseSeeds));
    nTq = sum(cat(1, UseInputControls{:}));
    nTh = sum(cat(1, UseDoseControls{:}));
    nT = nTk + nTs + nTq + nTh;
    sensitivities = true;
else
    sensitivities = false;
end

defaultVal = Inf;

if isempty(AbsTolY)
    AbsTolY = defaultVal;
end

AbsTolY = AbsTolY(:);

if isscalar(AbsTolY)
    AbsTolY = repmat(AbsTolY, ny, 1);
end

if numel(AbsTolY) == ny && sensitivities
    AbsTolY = [AbsTolY; repmat(defaultVal, ny*nT.^derorder, 1)];
end

if sensitivities
    assert(numel(AbsTolY) == ny+ny*nT, 'KroneckerBio:fixAbsTolY:AbsTolYSize', ...
        'AbsTolY must be empty, a scalar, length ny, or length ny+ny*nT.^(derivative order).')
else
    assert(numel(AbsTolY) == ny, 'KroneckerBio:fixAbsTolY:AbsTolYSize', ...
        'AbsTolY must be empty, a scalar, or of length ny.')
end

end