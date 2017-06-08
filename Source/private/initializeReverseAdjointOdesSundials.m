function [options_l, options_L, index_b] = initializeReverseAdjointOdesSundials(der, jac, L, tf, xf, lf, Lf, RelTol, AbsTol_l, delta_l, delta_L)

if isempty(delta_l)
    delta_l = @(t, x) 0;
end
if isempty(delta_L)
    delta_L = @(t, x) 0;
end

options_l = CVodeSetOptions('RelTol', RelTol,...
                           'AbsTol', AbsTol_l,...
                           'LinearSolver', 'Dense',...
                           'JacobianFn', jac);

lf = lf - delta_l(tf, xf, lf);
                       
index_b = CVodeInitB(der, 'BDF', 'Newton', tf, lf, options_l);

% Choosing not to error control the quadrature variables because we don't
% normally do so
options_L = CVodeQuadSetOptions('ErrControl',true,...
                               'RelTol',RelTol,'AbsTol',1e-9);

Lf = Lf - delta_L(tf, xf, lf);

%CVodeQuadInitB(index_b, L, Lf);
CVodeQuadInitB(index_b, L, Lf, options_L);

end