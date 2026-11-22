function [z, v, dz] = advectParticles(z, velFun, dt, maxHop)
%PF.ADVECTPARTICLES Move particles one time step with the flow.
%   [Z, V, DZ] = PF.ADVECTPARTICLES(Z, VELFUN, DT, MAXHOP) advances the
%   complex positions Z by DT with the midpoint rule. VELFUN(z) returns
%   u - i v. No particle moves further than MAXHOP in one step, so fast
%   particles near a pole do not jump over it. V is the velocity u + i v
%   used for the step and DZ the displacement.

    hop = @(d) d.*min(1, maxHop./max(abs(d), eps));
    zm = z + 0.5*hop(dt*conj(velFun(z)));
    v  = conj(velFun(zm));
    dz = hop(dt*v);
    z  = z + dz;
end
