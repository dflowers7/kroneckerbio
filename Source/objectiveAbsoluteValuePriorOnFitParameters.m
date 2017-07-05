function obj = objectiveAbsoluteValuePriorOnFitParameters(Tbar, Tslope, name, uselogbasis)
%obj = objectiveAbsoluteValuePriorOnFitParameters(Tbar, Tslope, name, uselogbasis)
% kbar and Vkbar are non-normalized as inputs

% Clean up inputs
if nargin < 4
    uselogbasis = [];
    if nargin < 3
        name = [];
        if nargin < 2
            Tslope = [];
        end
    end
end

% Check inputs
n = numel(Tbar);
if isempty(Tslope)
    Tslope = 1;
end
if isscalar(Tslope)
    Tslope = repmat(Tslope, n, 1);
end
assert(numel(Tslope) == n, 'KroneckerBio:objectiveAbsoluteValuePriorOnFitParameters:TslopeSize', 'Input "Tslope" must be empty, a scalar, or a vector of the same length as Tbar.')

% Fix log basis
if isempty(uselogbasis)
    uselogbasis = true;
end
if isscalar(uselogbasis)
    uselogbasis = repmat(uselogbasis, size(Tbar));
end
Tbar(uselogbasis) = log(Tbar(uselogbasis));

if isempty(name)
	name = 'AbsoluteValuePriorOnFitParameters';
end

% Set value below which a quadratic function will be used, to avoid
% nonsmoothness at 0
% Not implemented yet!
% fun = quadcoeff*dT^2 + c
% No linear term because slope needs to be zero at dT = 0
threshvals = repmat(0, n, 1); % Disabled when set to zeros
% 2*quadcoeff*threshval = Tslope (slope at threshval needs to be Tslope)
quadcoeffs = 0.5./threshvals.*Tslope;
% quadcoeff*threshval^2 + quadconst = Tslope*threshval (value at threshval needs
% to be Tslope*threshval)
quadconsts = Tslope.*threshvals - quadcoeffs.*threshvals.^2;

% Objective structure
obj.Type = 'Objective.Information.AbsoluteValuePriorOnFitParameters';
obj.Name = name;

obj.Continuous    = false;
obj.Complex       = false;
obj.Linked        = 0;
obj.DiscreteTimes = 0;

obj.G      = @G;
obj.dGdT   = @dGdT;
obj.d2GdT2 = @d2GdT2;
% obj.err    = @err;
% obj.derrdT = @derrdT;

% obj.p      = @p;
% obj.logp   = @logp;
% obj.F      = @F;
% obj.Fn     = @Fn;

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
        val = Tslope.*abs(diff);
        isquad = abs(diff) < threshvals;
        val(isquad) = quadcoeffs(isquad).*diff(isquad).^2 + quadconsts(isquad);
        val = sum(val);
    end

    function val = dGdT(int)

        T = int.T;

        % Normalize
        T(uselogbasis) = log(T(uselogbasis));

        diff = T - Tbar;
        
        val = Tslope.*sign(diff);
        
        isquad = abs(diff) < threshvals;
        val(isquad) = 2*quadcoeffs(isquad).*diff(isquad);

        if int.Normalized
            dTdlogT = T;
            val(~uselogbasis) = val(~uselogbasis).*dTdlogT(~uselogbasis);
        else
            dlogTdT = 1./T;
            val(uselogbasis) = val(uselogbasis).*dlogTdT(uselogbasis);
        end
        
        % Fix quadratic portion (not implemented!)
    end

    function val = d2GdT2(int)
         T = int.T;

        % Normalize
        T(uselogbasis) = log(T(uselogbasis));
       
        val = zeros(n);
        
        diff = T - Tbar;
        if int.Normalized
            d2TdlogT2 = -T.^2;
            val(~uselogbasis,~uselogbasis) = diag(Tslope(~uselogbasis).*sign(diff(~uselogbasis)).*d2TdlogT2(~uselogbasis));
        else
            d2logTdT2 = -1./T.^2;
            val(uselogbasis,uselogbasis) = diag(Tslope(uselogbasis).*sign(diff(uselogbasis)).*d2logTdT2(uselogbasis));
        end

        % Fix quadratic portion (not implemented!)
    end

% This function is not a valid least squares type function, so it does not
% have err and derrdT functions. You can't use the square root of the
% absolute value because the derivative is infinite at zero.

%     function val = err(int)
%         T = int.T;
%         
%         % Normalize
%         T(uselogbasis) = log(T(uselogbasis));
%         
%         diff = T - Tbar;
%         val = sqrt(Tslope.*abs(diff));
%     end
% 
%     function val = derrdT(int)
%         T = int.T;
%         T(uselogbasis) = log(T(uselogbasis));
%         diff = T - Tbar;
%         
%         val = 0.5*diag((Tslope.*abs(diff)).^-0.5);
%         if int.Normalized
%             val(uselogbasis,uselogbasis) = STbar(uselogbasis,uselogbasis);
%             dTdlogT = T(~uselogbasis);
%             val(~uselogbasis,~uselogbasis) = STbar(~uselogbasis,~uselogbasis)*diag(dTdlogT);
%         else
%             val(~uselogbasis,~uselogbasis) = STbar(~uselogbasis,~uselogbasis);
%             dlogTdT = 1./T(uselogbasis);
%             val(uselogbasis,uselogbasis) = STbar(uselogbasis,uselogbasis)*diag(dlogTdT);
%         end
%     end

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
        VTbark = Tslope(sol.UseParams,sol.UseParams);
        
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
        VTbark = Tslope(sol.UseParams,sol.UseParams);
        
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
        VTbark = Tslope(int.UseParams,int.UseParams);

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
