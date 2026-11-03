function [eps, nliq] = convection_multiplier(T, f, grid, Tm, conv)
%CONVECTION_MULTIPLIER  Effective-conductivity factor for the molten layer.
%   [eps, nliq] = convection_multiplier(T, f, grid, Tm, conv)
%   Raithby & Hollands (1975) correlation for a horizontal concentric annulus,
%   applied quasi-steadily to the melt layer between the tube wall and the melt
%   front (cells with f >= 0.5 contiguous from the wall):
%     k_eff/k = 0.386 (Pr/(0.861+Pr))^(1/4) Ra_c^(1/4),  k_eff/k >= 1
%     Ra_c = [ln(Do/Di)]^4 / (L^3 (Di^-3/5 + Do^-3/5)^5) * Ra_L
%   EG suppresses convection: the enhancement falls linearly to zero at w = w_s.
%   Applied to a full annular melt layer this overestimates real convection,
%   so it is used as the upper bound and conduction-only as the lower bound.
%   eps : B x 1 multiplier;  nliq : B x 1 number of liquid cells from the wall

liq = double(f >= 0.5);
nliq = sum(cumprod(liq, 2), 2);
ri = grid.rf(1);
s = grid.rf(nliq + 1); s = s(:);
delta = s - ri;
dT = max(T(:, 1) - Tm, 0);
Di = 2*ri; Do = 2*s;
RaL = conv.Ra_coef .* dT .* delta.^3;
Rac = log(Do ./ Di).^4 ./ (delta.^3 .* (Di^-0.6 + Do.^-0.6).^5) .* RaL;
ratio = 0.386 * (conv.Pr / (0.861 + conv.Pr))^0.25 * max(Rac, 0).^0.25;
ratio(nliq == 0 | ~isfinite(ratio)) = 1;
ratio = max(ratio, 1);
phi = min(max(1 - conv.w(:) / conv.w_s, 0), 1);
eps = 1 + phi .* (ratio - 1);
end
