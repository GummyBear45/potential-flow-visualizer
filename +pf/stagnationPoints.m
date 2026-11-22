function zs = stagnationPoints(Uc, zk, A, B, isDouble)
%PF.STAGNATIONPOINTS All zeros of a freestream-plus-poles velocity field.
%   ZS = PF.STAGNATIONPOINTS(UC, ZK, A, B, ISDOUBLE) returns every point
%   where w = UC + sum A/(z - ZK) + sum B/(z - ZK)^2 vanishes (see
%   PF.POLEVELOCITY). Multiplying w by all the denominators gives a
%   polynomial, so ROOTS finds every stagnation point at once. Element k
%   is a double pole when ISDOUBLE(k) is true; B(k) is used only then.
%   The caller filters the result (domain, distance from the poles).

    zs = zeros(0, 1);
    K = numel(zk);
    if K == 0, return; end
    owner = repelem(1:K, 1 + isDouble);     % double poles appear twice
    r = zk(owner);
    N = numel(r);
    P = Uc*poly(r);
    for k = 1:K
        i1 = find(owner == k, 1);
        t = A(k)*poly(r([1:i1-1, i1+1:N]));
        if isDouble(k)
            t2 = B(k)*poly(r(owner ~= k));
            t = t + [zeros(1, numel(t) - numel(t2)), t2];
        end
        P = P + [zeros(1, N+1-numel(t)), t];
    end
    if ~any(P), return; end
    zs = roots(P);
end
