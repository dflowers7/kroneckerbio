function m = LoadModelSimBioAnalytic(simbio, opts)
%LoadModelSimBioAnalytic Load analytic model from a Matlab SimBiology model
%   Modify the model and add outputs after calling this.
%
%   m = LoadModelSbmlMassAction(SimbioModel, opts)
%
%   Inputs
%   simbio: [ simbio model object ]
%       Matlab SimBiology Model object
%   opts: [ options struct scalar {} ]
%       .Verbose [ logical scalar {false} ]
%       	Print progress to command window

%% Clean up inputs
if nargin < 2
    opts = [];
end

% Default options
opts_.Verbose = 0;
opts_.Validate = false;
opts_.UseNames = false;

opts = mergestruct(opts_, opts);

%% Convert model

symbolic = simbio2symbolic(simbio, opts);

assert(isValidSymbolicModel(symbolic), 'LoadModelSbmlAnalytic:InvalidSymbolicModel', 'Symbolic model intermediate failed validation check')

m = symbolic2analytic(symbolic, opts);

end