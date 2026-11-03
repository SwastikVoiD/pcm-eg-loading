function [w_opt, E_opt, interior] = optimum_loading(w, E, w_min)
%OPTIMUM_LOADING  Loading that maximises stored energy.
%   [w_opt, E_opt, interior] = optimum_loading(w, E)          unconstrained
%   [w_opt, E_opt, interior] = optimum_loading(w, E, w_min)   w >= w_min
%   Uses the maximum of a shape-preserving (pchip) interpolant on a 1e-4 grid,
%   the same interpolant used for the design penalty, so both are consistent.
%   interior = false when the maximum sits on a bound (a corner solution).

if nargin < 3, w_min = -inf; end
w = w(:); E = E(:);
wf = (max(w(1), w_min):1e-4:w(end))';
Ef = interp1(w, E, wf, 'pchip');
[E_opt, i] = max(Ef);
w_opt = wf(i);
interior = i > 1 && i < numel(wf);
end
