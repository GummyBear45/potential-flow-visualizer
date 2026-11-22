function [u, v] = cappedArrows(w, cap, len)
%PF.CAPPEDARROWS Arrow components whose length follows speed, up to a cap.
%   [U, V] = PF.CAPPEDARROWS(W, CAP, LEN) turns complex velocities
%   W = u - i v into QUIVER components. An arrow is LEN long at speed CAP
%   and above and proportionally shorter below, so the arrows next to a
%   singularity do not swamp the plot. Non-finite velocities give no arrow.

    w(~isfinite(w)) = 0;
    sp = abs(w);
    sc = len*min(sp, cap)/cap./max(sp, eps);
    u = real(w).*sc;
    v = -imag(w).*sc;
end
