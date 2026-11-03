function pcm = pcm_properties(name)
%PCM_PROPERTIES  Thermophysical properties of the pure phase change materials.
%   pcm = pcm_properties('erythritol' | 'paraffin' | 'li2020_salt')
%
%   Fields: Tm [degC] melting point, L [J/kg] latent heat, cp [J/kgK],
%   rho [kg/m3], k [W/mK], dTm [K] width of the numerical mushy band.
%   Constant (phase-averaged) cp, rho and k are used, as in Li et al. (2020).
%
%   Sources (verify against the originals before publication):
%   erythritol  Tm = 117 C, L = 344 kJ/kg ........ Kaizawa et al., Heat Mass Transfer 44 (2008) 763-769
%               k  = 0.72 W/mK (pure) ............. Yuan et al., Renew. Energy 136 (2019) 211-222
%               cp = 2070 J/kgK (mean of solid 1380 and liquid 2760) and
%               rho = 1390 kg/m3 (mean of 1480 solid, 1300 liquid): commonly
%               reported values, e.g. Agyenim et al., Renew. Energy 36 (2011) 108-117  [VERIFY]
%   paraffin    n-docosane, the paraffin used by Sari & Karaipekli (2007)
%               Tm = 44 C (NIST triple point 315-317 K), L = 256 kJ/kg (Li, SEMSC 128, 2014)
%               k  = 0.22 W/mK ........................ Sari & Karaipekli (2007)
%               cp = 1950 J/kgK (mean of NIST 1.51 solid / 2.38 liquid kJ/kgK)
%               rho = 800 kg/m3 (approximate phase mean)  [VERIFY]
%   li2020_salt LiNO3-NaNO3-KCl, Table 1 of Li et al., Energy (2020) 117400

switch lower(name)
    case 'erythritol'
        pcm = struct('name', 'erythritol', 'Tm', 117.0, 'L', 344e3, ...
                     'cp', 2070.0, 'rho', 1390.0, 'k', 0.72, 'dTm', 1.0);
    case 'paraffin'
        pcm = struct('name', 'paraffin', 'Tm', 44.0, 'L', 256e3, ...
                     'cp', 1950.0, 'rho', 800.0, 'k', 0.22, 'dTm', 1.0);
    case 'li2020_salt'
        pcm = struct('name', 'li2020_salt', 'Tm', 160.0, 'L', 266e3, ...
                     'cp', 1330.0, 'rho', 2297.0, 'k', 0.88, 'dTm', 1.0);
    otherwise
        error('pcm_properties:unknown', 'Unknown PCM "%s".', name);
end
end
