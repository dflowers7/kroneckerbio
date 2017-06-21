function initializeOutputOdesSundials(der_y, y0, RelTol, AbsTol_quad)
% initializeOutputOdesSundials(der_y, y0, RelTol, AbsTol_quad)
% Initialize the zeroth- or first-order quadrature ODEs.

% Set quadrature options
useErrControl = ~(all(isinf(AbsTol_quad)) & isinf(RelTol));
options_quad = CVodeQuadSetOptions('ErrControl', useErrControl, 'RelTol', RelTol, 'AbsTol', AbsTol_quad);

% Initialize quadrature solver for outputs. Use this primarily for error
% control on output values.
CVodeQuadInit(der_y, y0, options_quad);

end