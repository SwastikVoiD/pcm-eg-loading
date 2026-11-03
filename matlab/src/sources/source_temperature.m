function T = source_temperature(t, src)
%SOURCE_TEMPERATURE  HTF temperature at time t for every case in the batch.
%   T = source_temperature(t, src)
%   src.Tmean     scalar mean temperature [degC]
%   src.A, src.P, src.phase, src.shape_id : B x 1
%   shape_id 0 = constant, 1 = sine, 2 = square wave (50 % duty)
%   All profiles share the same time-averaged temperature, so fluctuating and
%   constant sources offer the same mean driving force.

s = sin(2*pi*t ./ src.P(:) + src.phase(:));
sq = ones(size(s));
sq(s < 0) = -1;
y = zeros(size(s));
id = src.shape_id(:);
y(id == 1) = s(id == 1);
y(id == 2) = sq(id == 2);
T = src.Tmean + src.A(:) .* y;
end
