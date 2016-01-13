function int = integrateSysSimpSundials(m, con, tF, eve, fin, t_get, opts)

% Constants
nx = m.nx;

y = m.y;
u = con.u;

discontinuities = con.Discontinuities;
RelTol = opts.RelTol;
AbsTol = opts.AbsTol;

% Construct system
[der, jac, del] = constructSystem();

if ~con.SteadyState
    order = 0;
    ic = extractICs(m,con,opts,order);
else
    ic = steadystateSys(m, con, opts);
end

% Integrate f over time
sol = accumulateOdeFwdSundials(der, jac, 0, tF, ic, discontinuities, t_get, RelTol, AbsTol(1:nx), del, eve, fin, nx);

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

int.t = sol.t;
int.x = sol.x;
int.u = con.u(int.t);
int.y = y(int.t, int.x, int.u);

int.ie = sol.ie;
int.te = sol.te;
int.xe = sol.xe(1:nx,:);
int.ue = u(int.te);
int.ye = y(int.te, int.xe, int.ue);

int.sol = sol;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
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
end