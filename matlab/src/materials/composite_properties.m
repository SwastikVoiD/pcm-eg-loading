function props = composite_properties(pcm, w, beta)
%COMPOSITE_PROPERTIES  Effective properties of a PCM/expanded-graphite composite.
%   props = composite_properties(pcm, w, beta)
%   pcm  : struct from pcm_properties
%   w    : EG mass fraction(s) [-], any shape (returned as column vectors)
%   beta : conductivity slope [W/mK per unit mass fraction] (k_scenario_beta)
%
%   rho : inverse mass-weighted mixing rule (non-porous composite)
%   cp  : mass-weighted mixing rule
%   L   : (1 - w) * L_pcm  (latent heat carried by the PCM only; confirmed by
%         DSC in Sari & Karaipekli 2007 and Lee et al. 2014)
%   k   : k_pcm + beta * w  (linear percolating-network fit; Maxwell's
%         dispersed-sphere model is kept only as a lower bound, k_maxwell.m)

eg = eg_properties();
w = w(:);
props.w   = w;
props.rho = 1 ./ ((1 - w) / pcm.rho + w / eg.rho);
props.cp  = (1 - w) * pcm.cp + w * eg.cp;
props.L   = (1 - w) * pcm.L;
props.k   = pcm.k + beta * w;
props.Tm  = pcm.Tm;
props.dTm = pcm.dTm;
end
