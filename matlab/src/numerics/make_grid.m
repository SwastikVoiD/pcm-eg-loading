function grid = make_grid(geom, r_in, r_out, N)
%MAKE_GRID  Uniform finite-volume grid for the PCM layer.
%   grid = make_grid('radial' | 'planar', r_in, r_out, N)
%   Radial: annulus per unit tube length (V and face areas include 2*pi).
%   Planar: slab per unit face area.

rf = linspace(r_in, r_out, N + 1);
grid.rf = rf;
grid.rc = 0.5 * (rf(1:end-1) + rf(2:end));
grid.dr = rf(2) - rf(1);
grid.N = N;
grid.geom = geom;
switch lower(geom)
    case 'radial'
        grid.V  = pi * (rf(2:end).^2 - rf(1:end-1).^2);   % 1 x N  [m3/m]
        grid.Af = 2 * pi * rf;                              % 1 x N+1 [m2/m]
    case 'planar'
        grid.V  = grid.dr * ones(1, N);
        grid.Af = ones(1, N + 1);
    otherwise
        error('make_grid:geom', 'Unknown geometry "%s".', geom);
end
end
