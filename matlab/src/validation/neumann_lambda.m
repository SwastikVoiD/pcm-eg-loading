function lam = neumann_lambda(St_l, St_s)
%NEUMANN_LAMBDA  Root of the two-phase Neumann transcendental equation.
%   lam = neumann_lambda(St_l, St_s) for equal properties in both phases:
%   St_l*exp(-lam^2)/erf(lam) - St_s*exp(-lam^2)/erfc(lam) = lam*sqrt(pi)
%   St_l = cp (Tw - Tm)/L,  St_s = cp (Tm - Ti)/L.
g = @(x) St_l*exp(-x.^2)./erf(x) - St_s*exp(-x.^2)./erfc(x) - x*sqrt(pi);
lam = fzero(g, [1e-6 5]);
end
