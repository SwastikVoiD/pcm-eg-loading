function xf = melt_front_position(x, f)
%MELT_FRONT_POSITION  Position where the liquid fraction crosses 0.5.
%   xf = melt_front_position(x, f) for a single profile (vectors).
i = find(f < 0.5, 1, 'first');
if isempty(i), xf = x(end); return; end
if i == 1,    xf = x(1);   return; end
xf = x(i-1) + (0.5 - f(i-1)) * (x(i) - x(i-1)) / (f(i) - f(i-1));
end
