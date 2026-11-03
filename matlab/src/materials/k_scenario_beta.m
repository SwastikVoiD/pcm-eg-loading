function beta = k_scenario_beta(scenario)
%K_SCENARIO_BETA  Conductivity-enhancement slope for each preparation route.
%   beta = k_scenario_beta('low' | 'mid' | 'high')   [W/(m K) per unit EG mass fraction]
%
%   The composite conductivity is modelled as k = k_pcm + beta * w
%   (a percolating-network model; see composite_properties.m). The slope is set
%   by how the composite is made, not by the PCM, so three literature-anchored
%   scenarios bracket the uncertainty:
%
%   low   6.38  loose mixing/absorption: least-squares fit (intercept fixed at
%               k_pcm) to paraffin/EG data at 0-10 wt%  -> Sari & Karaipekli (2007)
%   mid  28.28  vacuum impregnation: paraffin/EG, 0.348 -> 6.003 W/mK at 20 wt%
%               -> Li Y. et al., ChemistrySelect 9 (2024)
%   high 117.9  impregnation + compression + sintering: erythritol/EG,
%               0.72 -> 12.51 W/mK at 10 wt% -> Yuan et al. (2019)
%   (fit_k_slope.m reproduces these numbers from the raw data.)

switch lower(scenario)
    case 'low',  beta = 6.38;
    case 'mid',  beta = 28.28;
    case 'high', beta = 117.9;
    otherwise
        error('k_scenario_beta:unknown', 'Unknown scenario "%s".', scenario);
end
end
