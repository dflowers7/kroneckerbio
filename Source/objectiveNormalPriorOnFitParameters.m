function obj = objectiveNormalPriorOnFitParameters(Tbar, VTbar, name, uselogbasis)
%obj = objectiveNormalPriorOnKineticParameters(kbar, Vkbar, name)
% kbar and Vkbar are non-normalized as inputs

% Clean up inputs
if nargin < 4
    uselogbasis = [];
    if nargin < 3
        name = [];
    end
end

% Check inputs
n = numel(Tbar);
assert(ismatrix(VTbar) && all(size(VTbar) == [n,n]), 'KroneckerBio:objectiveNormalPriorOnFitParameters:Vsize', 'Input "Vkbar" must be a square matrix of size numel(Tbar)-by-numel(Tbar)')
FTbar = infoinv(VTbar);
FSTbar = chol(FTbar);

% Fix log basis
if isempty(uselogbasis)
    uselogbasis = true;
end
if isscalar(uselogbasis)
    uselogbasis = repmat(uselogbasis, size(Tbar));
end
Tbar(uselogbasis) = log(Tbar(uselogbasis));

if isempty(name)
	name = 'NormalPriorOnFitParameters';
end

% Objective structure
obj.Type = 'Objective.Information.NormalPriorOnFitParameters';
obj.Name = name;

obj.Continuous    = false;
obj.Complex       = false;
obj.Linked        = 0;
obj.DiscreteTimes = 0;

obj.G      = @G;
obj.dGdT   = @dGdT;
obj.d2GdT2 = @d2GdT2;
obj.err    = @err;
obj.derrdT = @derrdT;

obj.p      = @p;
obj.logp   = @logp;
obj.F      = @F;
obj.Fn     = @Fn;

obj = pastestruct(objectiveZero, obj);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [val, stopTimes] = G(int)
        T = int.T;
        
        % Normalize
        T(uselogbasis) = log(T(uselogbasis));
        
        stopTimes = 0;
        diff = T - Tbar;
        val = diff.' * FTbar * diff;
    end

    function val = dGdT(int)

        T = int.T;

        % Normalize
        T(uselogbasis) = log(T(uselogbasis));

        val = 2*FTbar*(T-Tbar);

        if int.Normalized
            dTdlogT = T;
            val(~uselogbasis) = val(~uselogbasis).*dTdlogT(~uselogbasis);
        else
            dlogTdT = 1./T;
            val(uselogbasis) = val(uselogbasis).*dlogTdT(uselogbasis);
        end
    end

    function val = d2GdT2(int)
         T = int.T;

        % Normalize
        T(uselogbasis) = log(T(uselogbasis));
        
        val = zeros(size(FTbar));
        if int.Normalized
            val(uselogbasis,uselogbasis) = 2*FTbar(uselogbasis, uselogbasis);
            dTdlogT = T(~uselogbasis);
            d2TdlogT2 = -T(~uselogbasis).^2;
            dT = T - Tbar;
            dT = dT(~uselogbasis);
            val(~uselogbasis,~uselogbasis) = 2*diag(dTdlogT)*FTbar(~uselogbasis,~uselogbasis)*diag(dTdlogT) ...
                + diag((2*FTbar(~uselogbasis,~uselogbasis)*dT).*d2TdlogT2);
        else
            val(~uselogbasis,~uselogbasis) = 2*FTbar(~uselogbasis,~uselogbasis);
            dlogTdT = 1./T(uselogbasis);
            d2logTdT2 = -1./T(uselogbasis).^2;
            dT = T - Tbar;
            dT = dT(uselogbasis);
            val(uselogbasis,uselogbasis) = 2*diag(dlogTdT)*FTbar(uselogbasis,uselogbasis)*diag(dlogTdT) ...
                + diag((2*FTbar(uselogbasis,uselogbasis)*dT).*d2logTdT2);
        end

    end

    function val = err(int)
        T = int.T;
        T(uselogbasis) = log(T(uselogbasis));
        val = FSTbar*(T - Tbar);
    end

    function val = derrdT(int)
        T = int.T;
        T(uselogbasis) = log(T(uselogbasis));
        val = zeros(size(FSTbar));
        if int.Normalized
            val(uselogbasis,uselogbasis) = FSTbar(uselogbasis,uselogbasis);
            dTdlogT = T(~uselogbasis);
            val(~uselogbasis,~uselogbasis) = FSTbar(~uselogbasis,~uselogbasis)*diag(dTdlogT);
        else
            val(~uselogbasis,~uselogbasis) = FSTbar(~uselogbasis,~uselogbasis);
            dlogTdT = 1./T(uselogbasis);
            val(uselogbasis,uselogbasis) = FSTbar(uselogbasis,uselogbasis)*diag(dlogTdT);
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Information theory %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Likelihood function
    function val = p(sol)
        error('Not implemented yet!')

        % Old code from objectiveLogNormalPriorOnKineticParameters.m:

        nk = numel(sol.k);
        nTk = nnz(sol.UseParams);
        Tk = sol.k(sol.UseParams);
        Tbark = kbar(sol.UseParams);
        VTbark = VTbar(sol.UseParams,sol.UseParams);
        
        % Normalize
        logTk = log(Tk);
        logTbark = log(Tbark);
        VlogTbark = spdiags(Tbark.^(-1),0,nTk,nTk) * VTbark * spdiags(Tbark.^(-1),0,nTk,nTk);
        
        FlogTbark = infoinv(VlogTbark);

        diff = logTk - logTbark;
        val = (2*pi).^(-nTk/2) * det(VlogTbark).^(-1/2) * exp(-1/2 * diff.' * FlogTbark * diff);
    end

%% Log likelihood
    function val = logp(sol)
        error('Not implemented yet!')

        % Old code from objectiveLogNormalPriorOnKineticParameters.m:
        nk = numel(sol.k);
        nTk = nnz(sol.UseParams);
        Tk = sol.k(sol.UseParams);
        Tbark = kbar(sol.UseParams);
        VTbark = VTbar(sol.UseParams,sol.UseParams);
        
        % Normalize
        logTk = log(Tk);
        logTbark = log(Tbark);
        VlogTbark = spdiags(Tbark.^(-1),0,nTk,nTk) * VTbark * spdiags(Tbark.^(-1),0,nTk,nTk);
        
        FlogTbark = infoinv(VlogTbark);

        diff = logTk - logTbark;
        val = (-nTk/2) * log(2*pi) + -1/2 * sum(log(infoeig(VlogTbark))) + -1/2 * diff.' * FlogTbark * diff;
    end

%% Fisher information
    function val = F(int)
        error('Not implemented yet!')

        % Old code from objectiveLogNormalPriorOnKineticParameters.m:
        nTk = nnz(int.UseParams);
        VTbark = VTbar(int.UseParams,int.UseParams);

        if int.Normalized
            Tbark = kbar(int.UseParams);
            VTbark = spdiags(Tbark.^(-1),0,nTk,nTk) * VTbark * spdiags(Tbark.^(-1),0,nTk,nTk);
        end
        
        % This objective only provides information on the first nTk parameters
        T = [int.k(int.UseParams); int.s(int.UseSeeds); int.q(int.UseInputControls); int.h(int.UseDoseControls)];
        nT = numel(T);
        val = zeros(nT,nT);
        
        val(1:nTk,1:nTk) = infoinv(VTbark);
    end
end
