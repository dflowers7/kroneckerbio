function [intfun,objfun] = GenerateFIMEigenvalueFunction(eigindex, isconstraint)

if isconstraint
    intfun = {{
        @integrateAllSens,[]
        @integrateAllSens,[]
        @integrateAllSens,@integrateAllCurvVecProd
        @integrateAllSens,@integrateAllCurvVecProd
        }};
    objfun = {@fimeigenvalue_fmincon_constraint};
else
    intfun = {
        @integrateAllSens,[]
        @integrateAllSens,@integrateAllCurvVecProd
        };
    objfun = @fimeigenvalue_fmincon;
end


    function int = integrateAllCurvVecProd(m, con, obj, opts, objs, int)
        % Does not yet support complex integration
        
        % Get eigendirection from provided sensitivity integration
        [~,v] = getFeig(objs, int);
        
        % Extract T
        T = [m.k(opts.UseParams); con.s(opts.UseSeeds); con.q(opts.UseInputControls); con.h(opts.UseDoseControls)];
        %T = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
        
        % Add imaginary step
        v_mag = norm(v);
        v_dir = v./v_mag;
        step_mag = 1e-8;
        if opts.Normalized
            step = T.*v_dir*step_mag*1i;
        else
            step = v_dir*step_mag*1i;
        end
        T_dv = T + step;
        
        % Update parameters
        [m,con] = updateAll(m, con, T_dv, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
        
        % Integrate directional curvature
        int_dv = integrateAllSens(m, con, obj, opts);
        
        % Reformat fields and multiply by magnitude of v to get
        % curvature-vector product
        oldfields = {'dydT'};
        newfields = {'d2ydT2_v'};
        for i = 1:numel(oldfields)
            for j = 1:numel(int_dv)
                int_dv(j).(newfields{i}) = imag(int_dv(j).(oldfields{i}))/step_mag*v_mag;
            end
            int_dv = rmfield(int_dv, oldfields{i});
        end
        
        % Merge original int and new int
        int = mergestruct(int_dv, int);
        
        % Add v field
        [int.v] = deal(v);
        
    end

    function [G,D] = fimeigenvalue_fmincon(obj, int)
        
        [l,v] = getFeig(obj, int);
        G = -log(l);
        
        if nargout > 1
            % Collect derrdT and d2errdT2*v
            for i = numel(int):-1:1
                derrdT{i} = obj(i).derrdT(int(i));
                
                % Calculate d2errdT2*v by swapping d2ydT2*v in for dydT
                % and passing it into the same derrdT function
                int_temp = int(i);
                int_temp.dydT = int_temp.d2ydT2_v;
                d2errdT2_v{i} = obj(i).derrdT(int_temp);
            end
            derrdT = vertcat(derrdT{:});
            d2errdT2_v = vertcat(d2errdT2_v{:});
            
            D = 2*d2errdT2_v.'*derrdT*v;
%             if int(1).Normalized
%                 D = D + 2*(derrdT.'*derrdT)*v;
%             end
            D = -1./l*D;
        end
        
    end

    function [c,ceq,GC,GCeq] = fimeigenvalue_fmincon_constraint(obj, int)
        
        if nargout <= 2
            c = fimeigenvalue_fmincon(obj, int);
        else
            [c,GC] = fimeigenvalue_fmincon(obj, int);
            GCeq = zeros(size(GC,1),0);
        end
        ceq = zeros(0,1);

    end

    function [l,v] = getFeig(obj, int)
        for i = numel(int):-1:1
            derrdT{i} = obj(i).derrdT(int(i));
        end
        derrdT = vertcat(derrdT{:});
        F = derrdT.'*derrdT;
        [Feigvec,Feigval] = eig(F);
        [Feigval,sorti] = sort(diag(Feigval), 1, 'descend');
        Feigvec = Feigvec(:,sorti);
        l = Feigval(eigindex);
        v = Feigvec(:,eigindex);
    end

end