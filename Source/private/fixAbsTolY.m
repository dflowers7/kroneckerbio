function AbsTolY = fixAbsTolY(AbsTolY, ny, UseParams, UseSeeds, UseInputControls, UseDoseControls, derorder)

if nargin > 2
    nTk = sum(UseParams);
    nTs = sum(sum(UseSeeds));
    nTq = sum(cat(1, UseInputControls{:}));
    nTh = sum(cat(1, UseDoseControls{:}));
    nT = nTk + nTs + nTq + nTh;
end

defaultVal = Inf;

if isempty(AbsTolY)
    AbsTolY = defaultVal;
end

AbsTolY = AbsTolY(:);

if isscalar(AbsTolY)
    AbsTolY = repmat(AbsTolY, ny, 1);
end

if numel(AbsTolY) == ny && derorder >= 1
    AbsTolY = [AbsTolY; repmat(defaultVal, ny*nT, 1)];
end
if numel(AbsTolY) == ny + ny*nT && derorder >= 2
    AbsTolY = [AbsTolY; repmat(defaultVal, ny*nT*nT, 1)];
end

switch derorder
case 0
    nAbsTols = ny;
    nAbsTolStr = '';
case 1
    nAbsTols = ny + ny*nT;
    nAbsTolStr = ', ny+ny*nT';
case 2
    nAbsTols = ny + ny*nT + ny*nT*nT;
    nAbsTolStr = ', ny+ny*nT, ny+ny*nT+ny*nT*nT';
end

assert(numel(AbsTolY) == nAbsTols, 'KroneckerBio:fixAbsTolY:AbsTolYSize', ...
    'AbsTolY must be one of the following lengths: 0, 1, ny%s.', nAbsTolStr)

end