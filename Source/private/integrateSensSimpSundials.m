function int = integrateSensSimpSundials(m, con, tF, eve, fin, t_get, opts)

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
ng = 0;

nt = numel(t_get);

normalized = opts.Normalized;
T = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});

y = m.y;
u = con.u;
dudq = con.dudq;
dydx = m.dydx;
dydu = m.dydu;
dydk = m.dydk;

dkdT = sparse(find(opts.UseParams),1:nTk,1,nk,nT);

discontinuities = con.Discontinuities;
RelTol = opts.RelTol;
AbsTol_x = opts.AbsTol(1:nx);
AbsTol_dxdT = opts.AbsTol(nx+1:nx+nx*nT);

if ~con.SteadyState
    order = 1;
    ic = extractICs(m,con,opts,order);
else
    ic = steadystateSys(m, con, opts);
end

t0 = 0;

%%%%%% Construct state system %%%%%%%%%

[der_x, jac_x, del_x] = constructStateSystem();

x0_ = ic(1:nx);

% Initialize integrator

options_x = initializeStateOdesSundials(der_x, jac_x, t0, x0_, RelTol, AbsTol_x, del_x, eve);

%%%%% Construct sensitivity system %%%%%

[der_dxdT, del_dxdT] = constructStateSensitivitySystem();

dx0dT_ = full(reshape(ic(nx+1:nx+nx*nT), nx, nT));

options_dxdT = initializeStateSensitivityOdesSundials(der_dxdT, t0, x0_, dx0dT_, RelTol, AbsTol_dxdT, del_dxdT);

% Integrate f over time
freeMemoryOnFinish = true;
sol = accumulateOdeFwdSundials(options_x, options_dxdT, tF, t_get, discontinuities, x0_, dx0dT_, [], nx, nT, [], del_x, del_dxdT, eve, fin, freeMemoryOnFinish);

% Work down
int.Type = 'Integration.System.Simple';
int.Name = [m.Name ' in ' con.Name];

int.x_names = vec({m.States.Name});
int.u_names = vec({m.Inputs.Name});
int.y_names = vec({m.Outputs.Name});
int.k_names = vec({m.Parameters.Name});
int.s_names = vec({m.Seeds.Name});

int.nx = nx;
int.ny = m.ny;
int.nu = m.nu;
int.nk = m.nk;
int.ns = m.ns;
int.nq = con.nq;
int.nh = con.nh;
int.k = m.k;
int.s = con.s;
int.q = con.q;
int.h = con.h;

int.dydx = m.dydx;
int.dydu = m.dydu;

int.nT = nT;
int.Normalized = opts.Normalized;
int.UseParams = opts.UseParams;
int.UseSeeds = opts.UseSeeds;
int.UseInputControls = opts.UseInputControls;
int.UseDoseControls = opts.UseDoseControls;

int.t = sol.t;
int.x = sol.x;
int.u = con.u(int.t);
int.y = y(int.t, int.x, int.u);

int.dxdT = reshape(sol.dxdT, nx*nT, nt);
int.dudT = zeros(nu*nT,nt);
int.dydT = zeros(ny*nT,nt);
for it = 1:nt
    dudq_i = dudq(int.t(it)); % u_q
    dudq_i = dudq_i(:,opts.UseInputControls); % u_Q
    dudT_i = [sparse(nu,nTk+nTs), reshape(dudq_i, nu,nTq), sparse(nu,nTh)]; % u_Q -> u_T
    int.dudT(:,it) = vec(dudT_i); % u_T -> uT_
    
    dydx_i = dydx(int.t(it), int.x(:,it), int.u(:,it)); % y_x
    dydu_i = dydu(int.t(it), int.x(:,it), int.u(:,it)); % y_u
    dydk_i = dydk(int.t(it), int.x(:,it), int.u(:,it)); % y_k
    dxdT_i = reshape(int.dxdT(:,it), nx,nT); % xT_ -> x_T
    int.dydT(:,it) = vec(dydx_i * dxdT_i + dydu_i * dudT_i + dydk_i * dkdT); % y_x * x_T + y_u * u_T + y_k * k_T -> y_T -> yT_
end

if normalized
    % Normalize sensitivities
    int.dxdT = normalizeDerivatives(T, int.dxdT);
    int.dudT = normalizeDerivatives(T, int.dudT);
    int.dydT = normalizeDerivatives(T, int.dydT);
end

nte = numel(sol.ie);
int.ie = sol.ie;
int.te = sol.te;
int.xe = sol.xe(1:nx,:);
int.ue = u(int.te);
int.ye = y(int.te, int.xe, int.ue);

int.dxedT = reshape(sol.dxedT, nx*nT, nte);
int.duedT = zeros(nu*nT,nte);
int.dyedT = zeros(ny*nT,nte);
for it = 1:nte
    duedq_i = dudq(int.te(it)); % u_q
    duedq_i = duedq_i(:,opts.UseInputControls); % u_Q
    duedT_i = [sparse(nu,nTk+nTs), reshape(duedq_i, nu,nTq), sparse(nu,nTh)]; % u_Q -> u_T
    int.duedT(:,it) = vec(duedT_i); % u_T -> uT_

    dyedx_i = dydx(int.te(it), int.xe(:,it), int.ue(:,it)); % y_x
    dyedu_i = dydu(int.te(it), int.xe(:,it), int.ue(:,it)); % y_u
    dyedk_i = dydk(int.te(it), int.xe(:,it), int.ue(:,it)); % y_k
    dxedT_i = reshape(int.dxedT(:,it), nx,nT); % xT_ -> x_T
    int.dyedT(:,it) = vec(dyedx_i * dxedT_i + dyedu_i * duedT_i + dyedk_i * dkdT); % y_x * x_T + y_u * u_T + y_k * k_T -> y_T -> yT_
end

int.sol = sol;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructStateSystem()
        f     = m.f;
        dfdx  = m.dfdx;
        u     = con.u;
        d     = con.d;
        x0    = m.x0;
        nd    = m.ns;
        
        y = m.y;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of x with respect to time
        function [val, flag, new_data] = derivative(t, x, data)
            u_t = u(t);
            val = f(t, x, u_t);
            flag = 0;
            new_data = [];
        end
        
        % Jacobian of x derivative
        function [val, flag, new_data] = jacobian(t, x, data)
            u_t = u(t);
            val = full(dfdx(t, x, u_t));
            flag = 0;
        end
        
        % Dosing
        function val = delta(t, x)
            val = x0(d(t)) - x0(zeros(nd,1));
        end
    end

    function [der, del] = constructStateSensitivitySystem()
        
        dfdx = m.dfdx;
        dfdk = m.dfdk;
        dfdu = m.dfdu;
        
        der = @derivative;
        del = @delta;
        
        function [val, flag, new_data] = derivative(t, x, f, dxdT)
            u_i = u(t);
            val = dfdx(t, x, u_i)*dxdT + dfdT(t, x, u_i);
            flag = 0;
        end
        
        function val = delta(t, x)
            val = collectDoseImpact(m, con, t, 1, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
            % Convert from [x; dxdT(:)] to nx-by-nT dxdT
            val = reshape(val(nx+1:nx+nx*nT), nx, nT);
        end
        
        function val = dfdT(t, x, u)
            val = dfdk(t,x,u);
            dfdq = dfdu(t,x,u) * dudq(t);
            val = [val(:,opts.UseParams), sparse(nx,nTs), dfdq(:,opts.UseInputControls), sparse(nx,nTh)];
        end
    end
end