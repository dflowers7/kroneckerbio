function int = integrateAdjSundials(m, con, tF, eve, fin, t_get, opts, isObjectiveFunction, F_struct)

% Constants
nx = m.nx;
nu = m.nu;
ny = m.ny;
nk = m.nk;
nTk = sum(opts.UseParams);
nTs = sum(opts.UseSeeds);
nTq = sum(opts.UseInputControls);
nTh = sum(opts.UseDoseControls);
nT  = nTk + nTs + nTq + nTh;

u = con.u;

discontinuities = con.Discontinuities;
RelTol = opts.RelTol;
AbsTol_x = opts.AbsTol(1:nx);

if ~con.SteadyState
    order = 1;
    ic = extractICs(m,con,opts,order);
else
    ic = steadystateSys(m, con, opts);
end
t0 = 0;
x0_ = ic(1:nx);

g0 = F_struct.g0;
nF = F_struct.nF;
g = F_struct.g;
dgdT = F_struct.dgdT;

% Construct forward system RHS
[der_x, jac_x, del_x, quad_x] = constructForwardSystem();

% Initialize forward integration of states and quadrature part of function
options_fwd = initializeForwardAdjointOdesSundials(der_x, jac_x, quad_x, t0, x0_, g0, RelTol, AbsTol_x, del_x, eve);

% Integrate forward. State variables are x, quadratures are quadrature
% portion of provided function.
freeMemoryOnFinish = false;
if isempty(quad_x)
    nForwardQuads = 0;
else
    nForwardQuads = nF + nF*nT;
end
sol_fwd = accumulateOdeFwdSundials(options_fwd, [], tF, t_get, discontinuities, x0_, [], g0, nx, [], nForwardQuads, del_x, [], eve, fin, freeMemoryOnFinish);

% Extract state solutions
int.t = sol_fwd.t;
int.x = sol_fwd.x;
int.ie = sol_fwd.ie;
int.te = sol_fwd.te;
int.xe = sol_fwd.xe;

% Extract quadrature portion of function and its sensitivities from the
% quadratures of the forward integration. If values are empty, set them to
% zero-valued arrays.
nt = length(int.t);
nte = length(int.te);
if nForwardQuads == 0
    int.g = zeros(nF, nt);
    int.dgdT = zeros(nF, nT, nt);
    int.ge = zeros(nF, nte);
    int.dgedT = zeros(nF, nT, nte);
else
    int.g = sol_fwd.g(1:nF,:);
    int.dgdT = sol_fwd.g(nF+1:end,:);
    int.dgdT = reshape(int.dgdT, nF, nT, nt);
    int.ge = sol_fwd.ge(1:nF,:);
    int.dgedT = sol_fwd.ge(nF+1:end,:);
    int.dgedT = reshape(int.dgedT, nF, nT, nte);
end

% Evaluate pointwise part of function at each t_get
G = F_struct.G;
dGdT = F_struct.dGdT;
N = length(t_get);
int.G = zeros(nF, N);
int.dGdT = zeros(nF, nT, N);
for ii = 1:N
    int.G(:,ii) = G(sol_fwd.t(ii), sol_fwd.x(:,ii), u(sol_fwd.t(ii)));
    int.dGdT(:,:,ii) = dGdT(sol_fwd.t(ii), sol_fwd.x(:,ii), u(sol_fwd.t(ii)));
end
Ne = length(sol_fwd.ie);
int.Ge = zeros(nF, Ne);
int.dGedT = zeros(nF, nT, Ne);
for ii = 1:Ne
    int.Ge(:,ii) = G(sol_fwd.te(ii), sol_fwd.xe(:,ii), u(sol_fwd.te(ii)));
    int.dGedT(:,:,ii) = dGdT(sol_fwd.te(ii), sol_fwd.xe(:,ii), u(sol_fwd.te(ii)));
end

% Calculate function value at each t_get
int.F = int.g + int.G;
int.Fe = int.ge + int.Ge;

% Construct reverse system RHS
dgdx = F_struct.dgdx;
dGdx = F_struct.dGdx;
[der_l, jac_l, del_l, quad_l, del_quad_l] = constructReverseSystem();

% Initialize reverse integration
xF = int.x(:,end);
lF = vec(dGdx(tF, xF, u(tF)).');
LF = zeros(nF*nT, 1);
AbsTol_l = opts.AbsTol(nx+1:nx+nx*nF); % TODO: develop a better AbsTol interface. I'm not sure this even works as-is.
[options_rev, options_L, index_b] = initializeReverseAdjointOdesSundials(der_l, jac_l, quad_l, tF, xF, lF, LF, RelTol, AbsTol_l, del_l, del_quad_l);

% For objective functions, different time points need to be added together.
% For output sensitivities, different time points are recorded separately.
if isObjectiveFunction
    error('Not implemented yet!')
else
    
    int.l0 = zeros(nx, nF, N);
    int.L = zeros(nF, nT, N);
    for ii = N:-1:1
        
        ti = t_get(ii);
        xi = sol_fwd.x(:, ii);
        ui = u(ti);
        lf_i = dGdx(ti, xi, ui).';
        lf_i = lf_i - del_l(ti, xi, lf_i);
        Lf_i = -del_quad_l(ti, xi, lf_i);
        
        % Reinitialize backwards problem for this time point
        if ii < N
            CVodeReInitB(index_b, ti, lf_i, options_rev);
            CVodeQuadReInitB(index_b, Lf_i, options_L);
        end
        
        % Integrate backwards. State variables are lambdas (referred to as
        % l), quadratures are integral(l.'*df(t)dT, 0, t)
        freeMemoryOnFinish = false;
        sol_rev = accumulateOdeRevSundials(options_rev, index_b, t0, ti, t0, discontinuities, lf_i, Lf_i, nx*nF, nT*nF, del_l, del_quad_l, freeMemoryOnFinish);
        
        % Extract lambda values at t0 (times are sorted backward in
        % accumulateOdeRevSundials)
        int.l0(:, :, ii) = reshape(sol_rev.x(:, end)-del_l(t0, sol_fwd.x(:,1), sol_rev.x(:,end)), nx, nF);
        
        % Extract integral term values. Negative because we integrated
        % backward in time.
        int.L(:, :, ii) = -(reshape(sol_rev.g(:, end)-del_quad_l(t0, sol_fwd.x(:,1), sol_rev.x(:,end)), nF, nT));
        
    end
    
    int.l0e = zeros(nx, nF, Ne);
    int.Le = zeros(nF, nT, Ne);
    for ii = 1:Ne
        
        ti = int.te(ii);
        xi = int.xe(:, ii);
        ui = u(ti);
        lf_i = dGdx(ti, xi, ui);
        lf_i = lf_i - del_l(ti, xi, lf_i);
        Lf_i = -del_quad_l(ti, xi, lf_i);
        
        % Reinitialize backwards problem for this time point
        %if ii > 1
            CVodeReInitB(index_b, ti, lf_i, options_rev);
            CVodeQuadReInitB(index_b, Lf_i);
        %end
        
        % Integrate backwards. State variables are lambdas (referred to as
        % l), quadratures are integral(l.'*df(t)dT, 0, t)
        freeMemoryOnFinish = false;
        sol_rev = accumulateOdeRevSundials(options_rev, index_b, t0, t0, discontinuities, lf_i, Lf_i, nx*nF, nT*nF, del_l, del_quad_l, freeMemoryOnFinish);
        
        % Extract lambda values at t0 (times are sorted backward in
        % accumulateOdeRevSundials)
        int.l0e(:, :, ii) = reshape(sol_rev.x(:, end), nx, nF);
        
        % Extract integral term values. Negative because we integrated
        % backward in time.
        int.Le(:, :, ii) = -reshape(sol_rev.g(:, end), nF, nT);
        
    end
    
end

CVodeFree;

% Add together terms to calculate sensitivity
dx0dT = reshape(ic(nx+1:nx+nx*nT), nx, nT);
l0T_times_dx0dT = permute(reshape(reshape(int.l0, nx, nF*N).'*dx0dT, nF, N, nT), [1 3 2]);
int.dFdT = int.dGdT + int.dgdT + l0T_times_dx0dT + int.L;
l0eT_times_dx0dT = permute(reshape(reshape(int.l0e, nx, nF*Ne).'*dx0dT, nF, Ne, nT), [1 3 2]);
int.dFedT = int.dGedT + int.dgedT + l0eT_times_dx0dT + int.Le;


% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del, quad] = constructForwardSystem()
        f     = m.f;
        dfdx  = m.dfdx;
        d     = con.d;
        x0    = m.x0;
        nd    = m.ns;
        
        y = m.y;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        if isempty(dgdT) && isempty(g)
            quad = [];
        else
            quad = @quadrature;
        end
        
        % Derivative of x with respect to time
        function [val, flag] = derivative(t, x)
            u_t = u(t);
            val = f(t, x, u_t);
            flag = 0;
        end
        
        % Jacobian of x derivative
        function [val, flag] = jacobian(t, x, dxdt)
            u_t = u(t);
            val = full(dfdx(t, x, u_t));
            flag = 0;
        end
        
        % Dosing
        function val = delta(t, x)
            val = x0(d(t)) - x0(zeros(nd,1));
        end
        
        % Quadrature
        function [val, flag] = quadrature(t, x)
            % This function, evaluated as a pure quadrature in the forward
            % pass, calculates the part of the derivative of g with respect
            % to T that doesn't depend on dxdT.
            u_t = u(t);
            val = vec([g(t, x, u_t) dgdT(t, x, u_t)]);
            flag = 0;
        end
        
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating lambda %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [der, jac, del, quad, del_quad] = constructReverseSystem()
        
        dfdx  = m.dfdx;
        dfdu = m.dfdu;
        dfdk = m.dfdk;
        dudq  = con.dudq;
        
        der = @derivative;
        jac = @jacobian;
        % If an objective function, then we use delta to add in the
        % discrete part of the objective function at each time point of the
        % objective function. If output sensitivities, the only time this
        % term appears is at tF, so we just set the final condition to the
        % corrected value, avoiding the need to create a new delta for each
        % requested time point.
        if isObjectiveFunction
            del = @delta;
        else
            del = @(t,x,l) 0;
        end
        quad = @quadrature;
        del_quad = @quadrature_delta;
        
        function [val, flag] = derivative(t, x, l)
            u_i = u(t);
            if ~isempty(dgdx)
                val = -dfdx(t, x, u_i).'*l - dgdx(t, x, u_i).';
            else
                val = -dfdx(t, x, u_i).'*l;
            end
            val = vec(val);
            flag = 0;
        end
        
        function [val, flag] = jacobian(t, x, l, dldt)
            u_i = u(t);
            val = full(kron(eye(nF), -dfdx(t, x, u_i).'));
            flag = 0;
        end
        
        function val = delta(t, x, l)
            if ismember(t, t_get)
                u_t = u(t);
                val = vec(-dGdx(t, x, u_t).'); % Not sure about the negative sign here
            else
                val = 0;
            end
        end
        
        function [val, flag] = quadrature(t, x, l)
            u_t = u(t);
%             delta_x = collectDoseImpact(m, con, t, 0, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
%             val = vec(l.'*dfdT(t, x-delta_x, u_t))
            val = vec(l.'*dfdT(t, x, u_t));
            flag = 0;
        end
        
        function val = quadrature_delta(t, x, l)
            deltax_joint = collectDoseImpact(m, con, t, 1, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
            dxdT = reshape(deltax_joint(nx+1:nx+nx*nT), nx, nT);
            val = vec(l.'*dxdT);
        end
        
        function val = dfdT(t, x, u_t)
            dfdk_t = dfdk(t, x, u_t);
            dfdq_t = dfdu(t, x, u_t)*dudq(t);
            val = [dfdk_t(:, opts.UseParams) zeros(nx, nTs) dfdq_t(:, opts.UseInputControls) zeros(nx, nTh)];
        end
        
    end

end