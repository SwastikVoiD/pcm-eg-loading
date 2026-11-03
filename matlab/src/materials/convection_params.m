function conv = convection_params(name, w, w_s, ra_mult)
%CONVECTION_PARAMS  Inputs for the natural-convection enhancement model.
%   conv = convection_params(pcm_name, w)                 w_s = 0.10, ra_mult = 1
%   conv = convection_params(pcm_name, w, w_s, ra_mult)
%   w       : EG mass fraction of every batch member
%   w_s     : loading at which EG fully suppresses melt convection (form-stable limit)
%   ra_mult : multiplier on the Rayleigh number (property-uncertainty sensitivity)
if nargin < 3, w_s = 0.10; end
if nargin < 4, ra_mult = 1; end
L = liquid_properties(name);
nu = L.mu / L.rho; alpha = L.k / (L.rho * L.cp);
conv = struct('Pr', nu/alpha, 'Ra_coef', 9.81*L.beta/(nu*alpha)*ra_mult, 'w', w(:), 'w_s', w_s);
end
