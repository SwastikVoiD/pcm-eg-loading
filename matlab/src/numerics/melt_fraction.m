function fl = melt_fraction(f, grid)
%MELT_FRACTION  Volume-averaged liquid fraction, per case (B x 1).
fl = (f * grid.V(:)) / sum(grid.V);
end
