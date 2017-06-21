function relTol = fixRelTolY(relTol)

assert(isnumeric(relTol), 'KroneckerBio:RelTolY:Numeric', 'RelTol must be numeric')
assert(numel(relTol) <= 1, 'KroneckerBio:RelTolY:Length', 'RelTol cannot be provided as a vector')

if isempty(relTol)
    relTol = Inf;
end
