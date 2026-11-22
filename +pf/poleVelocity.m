function w = poleVelocity(z, Uc, zk, A, B)
%PF.POLEVELOCITY Complex velocity of a freestream plus simple and double poles.
%   W = PF.POLEVELOCITY(Z, UC, ZK, A, B) returns u - i v at the points Z for
%       w = UC + sum A(k)/(z - ZK(k)) + sum B(k)/(z - ZK(k))^2
%   UC is the complex freestream U*exp(-1i*alpha). For a source of strength m
%   A = m/(2*pi), for a vortex of circulation Gamma (counter-clockwise
%   positive) A = -1i*Gamma/(2*pi), for a doublet of strength kappa and axis
%   angle alpha B = -kappa*exp(1i*alpha)/(2*pi).

    w = Uc + zeros(size(z));
    for k = 1:numel(zk)
        d = z - zk(k);
        w = w + A(k)./d + B(k)./d.^2;
    end
end
