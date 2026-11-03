function [props, src, meta] = build_batch(pcm, beta, cfg)
%BUILD_BATCH  Every (loading, source case) pair as one solver batch.
%   [props, src, meta] = build_batch(pcm, beta, cfg)
%   meta.w, meta.case index each batch member; meta.cases is the case table.

cases = build_cases(cfg);
nw = numel(cfg.w); nc = size(cases, 1);
[ww, cc] = ndgrid(cfg.w(:), 1:nc);       % loading varies fastest
ww = ww(:); cc = cc(:);
props = composite_properties(pcm, ww, beta);
src.Tmean    = pcm.Tm + cfg.dT_super;
src.shape_id = cases(cc, 1);
src.A        = cases(cc, 2);
src.P        = cases(cc, 3);
src.phase    = cases(cc, 4);
meta.w = ww; meta.case = cc; meta.cases = cases; meta.nw = nw; meta.nc = nc;
end
