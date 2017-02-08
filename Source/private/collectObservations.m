function [tF, eve, fin, t_get] = collectObservations(m, con, obs)
nx = m.nx;
u = con.u;
y = m.y;

n_obs = numel(obs);
nes = vec([obs.ne]);
n_eve = sum(nes);

iszeroobs = strcmp({obs.Type}, 'Objective.Data.Zero');
obslist_i = find(~iszeroobs);
obsnonzero = obs(obslist_i);
n_nonzeroobs = numel(obslist_i);
nes_nonzeroobs = vec([obsnonzero.ne]);

tF = max([obs.tF]);
eve = @events_combined;
fin = @is_finished_combined;
t_get = unique([obs.DiscreteTimes]);

    function [val, is_terminal, direction] = events_combined(t, joint)
        x = joint(1:nx);
        
        val = zeros(n_eve,1);
        is_terminal = zeros(n_eve,1);
        direction = zeros(n_eve,1);
        
        u_t = u(t);
        y_t = y(t,x,u_t);
        
        i_eve = 0;
        
        % To avoid significant overhead of evaluating obs.Events, only
        % calculate on objectives that are non-zero
        for iobsnonzero = 1:n_nonzeroobs
            i_start = i_eve + 1;
            i_end = i_eve + nes_nonzeroobs(iobsnonzero);
            i_int = i_start:i_end;
            [val(i_int), is_terminal(i_int), direction(i_int)] = obsnonzero(iobsnonzero).Events(t,y_t);
        end
    end

    % Is finished
    function finished = is_finished_combined(sol)
        finished = true;
        for iobs = 1:n_obs
            if obs(iobs).IsFinished(sol)
                finished = false;
                break
            end
        end
    end
end
