function E = stored_energy(T, f, props, grid, Tref)
%STORED_ENERGY  Sensible + latent energy relative to Tref, per case [J/m].
%   E = stored_energy(T, f, props, grid, Tref)   T, f are B x N; E is B x 1
E = sum((props.rho(:) .* grid.V) .* (props.cp(:) .* (T - Tref) + f .* props.L(:)), 2);
end
