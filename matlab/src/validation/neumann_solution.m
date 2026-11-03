function [T, s] = neumann_solution(x, t, Tw, Ti, Tm, alpha, lam)
%NEUMANN_SOLUTION  Analytical two-phase melting of a semi-infinite slab.
%   [T, s] = neumann_solution(x, t, Tw, Ti, Tm, alpha, lam)
%   Wall at x = 0 held at Tw, initial temperature Ti, equal properties.
%   Returns the temperature profile T(x) and front position s at time t.
eta = x ./ (2*sqrt(alpha*t));
s = 2*lam*sqrt(alpha*t);
T = zeros(size(x));
liq = x <= s;
T(liq)  = Tw - (Tw - Tm) * erf(eta(liq)) / erf(lam);
T(~liq) = Ti + (Tm - Ti) * erfc(eta(~liq)) / erfc(lam);
end
