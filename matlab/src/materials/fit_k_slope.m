function [beta, data] = fit_k_slope()
%FIT_K_SLOPE  Reproduce the three conductivity slopes from published data.
%   [beta, data] = fit_k_slope() returns beta.low/mid/high and the raw points.

% Sari & Karaipekli (2007): paraffin (n-docosane)/EG, measured k
data.sari.w = [0 0.02 0.04 0.07 0.10];
data.sari.k = [0.22 0.40 0.52 0.68 0.82];
% Li Y. et al. (2024): paraffin wax/EG, vacuum adsorption
data.li2024.w = [0 0.20];
data.li2024.k = [0.348 6.003];
% Yuan et al. (2019): erythritol/EG, compressed and sintered (17.38x at 10 wt%)
data.yuan2019.w = [0 0.10];
data.yuan2019.k = [12.51/17.38 12.51];

% least squares with the intercept fixed at the pure-PCM value
w = data.sari.w(:); dk = data.sari.k(:) - data.sari.k(1);
beta.low  = (w' * dk) / (w' * w);
beta.mid  = diff(data.li2024.k)   / diff(data.li2024.w);
beta.high = diff(data.yuan2019.k) / diff(data.yuan2019.w);
end
