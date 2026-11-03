function k = k_maxwell(km, kf, phi)
%K_MAXWELL  Maxwell-Eucken conductivity of dispersed spheres (lower bound).
%   k = k_maxwell(km, kf, phi): matrix km, filler kf, filler volume fraction phi.
k = km .* (kf + 2*km + 2*phi.*(kf - km)) ./ (kf + 2*km - phi.*(kf - km));
end
