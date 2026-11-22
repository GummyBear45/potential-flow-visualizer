function [x, y] = traceStreamlines(velFun, z0, sg, h, nstep, keepFun)
%PF.TRACESTREAMLINES Streamlines through a velocity field, by arc length.
%   [X, Y] = PF.TRACESTREAMLINES(VELFUN, Z0, SG, H, NSTEP, KEEPFUN) integrates
%   one streamline from every seed in the row Z0 (complex positions) with
%   the midpoint rule and step length H along the curve.
%     VELFUN(z)   returns u - i v at the complex points z
%     SG          row of +1 (follow the flow) or -1 (go upstream)
%     KEEPFUN(z)  returns a logical row, false where a line must end
%   X and Y hold all lines one after another, separated by NaN, ready for
%   a single LINE object.

    n = numel(z0);
    if n == 0, x = NaN; y = NaN; return; end
    Z = nan(nstep+2, n);                % last row stays NaN: line separator
    Z(1,:) = z0;  z = z0;  alive = true(1, n);
    for i = 1:nstep
        zm = z + 0.5*h*sg.*unitDir(velFun(z));
        z  = z + h*sg.*unitDir(velFun(zm));
        z(~alive) = NaN;
        Z(i+1,:) = z;
        alive = alive & keepFun(z);
        if ~any(alive), break; end
    end
    x = real(Z(:));  y = imag(Z(:));
end

function v = unitDir(w)
    % unit vector along the flow, as a complex number
    v = conj(w)./abs(w);
end
