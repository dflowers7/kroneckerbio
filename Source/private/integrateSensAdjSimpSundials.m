function int = integrateSensAdjSimpSundials(m, con, tF, eve, fin, t_get, opts)

% Constants
nx = m.nx;
nu = m.nu;
ny = m.ny;
nk = m.nk;
nq = con.nq;
nTk = sum(opts.UseParams);
nTs = sum(opts.UseSeeds);
nTq = sum(opts.UseInputControls);
nTh = sum(opts.UseDoseControls);
nT  = nTk + nTs + nTq + nTh;

AdjointOutputSensitivities = opts.AdjointOutputSensitivities;
nY = sum(AdjointOutputSensitivities);

% Set up dkdT
dkdT = zeros(nk, nT);
useParamIndices = sub2ind([nk, nT], vec(find(opts.UseParams)), vec(1:nTk));
dkdT(useParamIndices) = 1;

% Set up dqdT
dqdT = zeros(nq, nT);
useInputControlIndices = sub2ind([nq, nT], vec(find(opts.UseInputControls)), vec(nTk+nTs+1:nTk+nTs+nTq));
dqdT(useInputControlIndices) = 1;

% Copy model and experiment function handles needed
dydu = m.dydu;
dydk = m.dydk;
dudq = con.dudq;

% Set up functions for outputs
g = [];
dgdx = [];
dgdT = [];
g0 = [];
nF = nY;
G = @G_fun;
dGdx = @dGdx_fun;
dGdT = @dGdT_fun;
F_struct = struct('g', g, 'dgdx', dgdx, 'dgdT', dgdT, 'g0', g0, 'nF', nF, 'G', G, 'dGdx', dGdx, 'dGdT', dGdT);

    function val = G_fun(t,x,u)
        val = m.y(t,x,u);
        val = full(val(AdjointOutputSensitivities,:));
    end

    function val = dGdx_fun(t,x,u)
        val = m.dydx(t,x,u);
        val = full(val(AdjointOutputSensitivities,:));
    end

    function val = dGdT_fun(t, x, u)
        % Returns non-state-dependent terms in full derivative of y
        dGdu = dydu(t, x, u);
        dGdu = dGdu(AdjointOutputSensitivities,:);
        dGdk = dydk(t, x, u);
        dGdk = dGdk(AdjointOutputSensitivities,:);
        dudT = dudq(t)*dqdT;
        val = full(dGdu*dudT + dGdk*dkdT);
    end

% Call adjoint integrator
isObjectiveFunction = false;
int_temp = integrateAdjSundials(m, con, tF, eve, fin, t_get, opts, isObjectiveFunction, F_struct);

% Extract solutions
int.t = int_temp.t;
int.x = int_temp.x;
int.y = int_temp.F;
int.ie = int_temp.ie;
int.te = int_temp.te;
int.xe = int_temp.xe;
int.ye = int_temp.Fe;
nt = length(int.t);
int.dydT = reshape(int_temp.dFdT, nF*nT, nt);
nte = length(int.te);
int.dyedT = reshape(int_temp.dFedT, nF*nT, nte);

int.Type = 'Integration.AdjointSensitivity.Simple';
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

int.u = con.u(int.t);
int.dudT = zeros(nu*nT,nt);
for it = 1:nt
    dudq_i = dudq(int.t(it)); % u_q
    dudq_i = dudq_i(:,opts.UseInputControls); % u_Q
    dudT_i = [sparse(nu,nTk+nTs), reshape(dudq_i, nu,nTq), sparse(nu,nTh)]; % u_Q -> u_T
    int.dudT(:,it) = vec(dudT_i); % u_T -> uT_
end

normalized = opts.Normalized;
T = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});
if normalized
    % Normalize sensitivities
    int.dudT = normalizeDerivatives(T, int.dudT);
    int.dydT = normalizeDerivatives(T, int.dydT);
end

int.ue = con.u(int.te);
int.duedT = zeros(nu*nT,nte);
for it = 1:nte
    duedq_i = dudq(int.te(it)); % u_q
    duedq_i = duedq_i(:,opts.UseInputControls); % u_Q
    duedT_i = [sparse(nu,nTk+nTs), reshape(duedq_i, nu,nTq), sparse(nu,nTh)]; % u_Q -> u_T
    int.duedT(:,it) = vec(duedT_i); % u_T -> uT_
end

if normalized
    % Normalize sensitivities
    int.duedT = normalizeDerivatives(T, int.duedT);
    int.dyedT = normalizeDerivatives(T, int.dyedT);
end

int.dxdT = zeros(0, nt);
int.dxedT = zeros(0, nte);

int.sol = int_temp;

end