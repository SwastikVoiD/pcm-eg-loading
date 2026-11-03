function eg = eg_properties()
%EG_PROPERTIES  Expanded graphite skeleton properties.
%   rho = 2200 kg/m3 (graphite true density), cp = 710 J/kgK (graphite at room
%   temperature; slightly low at 120 C, a minor effect at w <= 0.25),
%   k_particle = 150 W/mK (only used by the Maxwell lower bound).
eg = struct('rho', 2200.0, 'cp', 710.0, 'k_particle', 150.0);
end
