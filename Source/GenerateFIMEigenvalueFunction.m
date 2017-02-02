function [intfun,objfun] = GenerateFIMEigenvalueFunction(eigindex, isconstraint)
% WARNING: Currently there is a flaw in the design that does not allow
% multiple eigenvalue functions to be used in the same optimization with different
% eigenvalue indices. The cause is that the memoization scheme depends on
% the function name to distinguish between separate integration calls.
% Since different eigindex values use the same function name, each
% call overwrites the memoized results of the others. I'm not fixing this
% right now because I have no intent to use multiple eigindex values in the
% same optimization.

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


    function int_dv_all = integrateAllCurvVecProd(m, con, obj, opts, objs_cell, int_all_cell)
        % Does not yet support complex integration
        % objs_cell: cell array of obj structs
        % int_all_cell: cell array of int structs
        % The cell arrays should be the same size. They are evaluated one
        % at a time, matching each objs with each int_all.
        % Resulting int_dv's will be vertically concatenated to return a
        % final result similar to that returned by other integrate*()
        % functions.
        
        assert(iscell(objs_cell), 'KroneckerBio:integrateAllCurvVecProd:ObjsNotCell', ...
            'objs_cell should be a cell vector. This is a bug in GenerateObjective.')
        nEvals = numel(objs_cell);
        int_dv_cell = cell(nEvals,1);
        for k = 1:nEvals
            
            int_all = int_all_cell{k};
            objs = objs_cell{k};
            
            % Get eigendirection from provided sensitivity integration
            [~,v] = getFeig(objs, int_all);
            
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
            
            % Add curvature-vector product to int struct and only keep real
            % components of other fields
            fields = fieldnames(int_dv);
            for i = 1:numel(fields)
                for j = 1:numel(int_dv)
                    switch fields{i}
                        case 'dydT'
                            int_dv(j).d2ydT2_v = imag(int_dv(j).dydT)/step_mag*v_mag;
                            int_dv(j).dydT = real(int_dv(j).dydT);
                        otherwise
                            if isnumeric(int_dv(j).(fields{i}))
                                int_dv(j).(fields{i}) = real(int_dv(j).(fields{i}));
                            end
                    end
                end
                %int_dv = rmfield(int_dv, oldfields{i});
            end
            
            % Merge original int and new int
            %int = mergestruct(int_dv, int);
            
            % Add v field
            [int_dv.v] = deal(v);
            
            int_dv_cell{k} = int_dv;
        end
        
        int_dv_all = vertcat(int_dv_cell{:});
        
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