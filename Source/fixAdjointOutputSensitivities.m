function AdjointOutputSensitivities = fixAdjointOutputSensitivities(AdjointOutputSensitivities, ny, nCon)
%fixAdjointOutputSensitivities A helper function for defining which outputs
%   will have sensitivities calculated via the adjoint method. It takes a
%   variety of ways to specify which outputs and standardizes them into
%   matrix of logicals indicating which output sensitivities are to be
%   calculated per experiment.
%
%   AdjointOutputSensitivities = AdjointOutputSensitivities(AdjointOutputSensitivities, ny, nCon)
%
%   Inputs
%   AdjointOutputSensitivities: [ logical matrix ny by nCon | positive integer vector ]
%       1) linear index vector into ny, assumed same for all conditions
%       2) logical index vector ny to indicate that all conditions have the
%           same active parameters
%       3) matrix of logical indexes size ny by nCon
%   ny: [ nonegative integer scalar ]
%       The number of outputs in the model
%   nCon: [ nonnegative integer ]
%       Number of experimental conditions
%
%   Outputs
%   UseSeeds: [ logical matrix ns by nCon ]
%       If UseModelSeeds = true, then UseSeeds will be a logical column
%       vector. Otherwise, it will be a logical matrix.
%   nTs: [ nonnegative integer ]
%       Number of active seed parameters

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if isempty(AdjointOutputSensitivities)
    AdjointOutputSensitivities = false(ny, nCon);
elseif isnumeric(AdjointOutputSensitivities) && all(floor(vec(AdjointOutputSensitivities)) == vec(AdjointOutputSensitivities)) && all(AdjointOutputSensitivities >= 1)
    % Linear index
    AdjointOutputSensitivities = vec(AdjointOutputSensitivities);
    assert(all(AdjointOutputSensitivities <= ny), 'KroneckerBio:AdjointOutputSensitivities:LinearIndexOutOfRange', 'AdjointOutputSensitvities, when a linear index, can have no value larger than m.ny. Use a logical matrix to refer to different seeds on each condition.')
    temp = false(ny,nCon);
    temp(AdjointOutputSensitivities,:) = true;
    AdjointOutputSensitivities = temp;
elseif islogical(AdjointOutputSensitivities)
    % Logical index
    assert(numel(AdjointOutputSensitivities) == ny*nCon, 'KroneckerBio:AdjointOutputSensitivities:InvalidLogicalSize', 'AdjointOutputSensitvities, when a logical index, must have a number of elements equal to ny*nCon')
    AdjointOutputSensitivities = reshape(AdjointOutputSensitivities, ny,nCon);
else
    error('KroneckerBio:UseSeeds:InvalidType', 'AdjointOutputSensitvities must be provided as logical or linear index into m.y')
end
