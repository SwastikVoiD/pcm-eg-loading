function [E, fl, j] = energy_at_window(res, t_c)
%ENERGY_AT_WINDOW  Stored energy and melt fraction at the end of a charging window.
%   [E, fl, j] = energy_at_window(res, t_c)  -> B x 1 vectors at res.t(j) ~ t_c
[~, j] = min(abs(res.t - t_c));
E  = res.E(:, j);
fl = res.fl(:, j);
end
