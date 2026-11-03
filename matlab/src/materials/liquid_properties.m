function L = liquid_properties(name)
%LIQUID_PROPERTIES  Molten-PCM properties for the natural-convection model.
%   k [W/mK], rho [kg/m3], cp [J/kgK], mu [Pa s], beta [1/K] (thermal expansion).
%   k, rho, cp: liquid values from the sources listed in pcm_properties.m.
%   mu and beta are typical literature values [VERIFY]; their uncertainty is
%   bracketed by the Rayleigh multiplier in convection_params.m.
switch lower(name)
    case 'erythritol'
        L = struct('k', 0.326, 'rho', 1300, 'cp', 2760, 'mu', 0.016,  'beta', 5.0e-4);
    case 'paraffin'
        L = struct('k', 0.15,  'rho', 778,  'cp', 2380, 'mu', 0.0036, 'beta', 9.0e-4);
    otherwise
        error('liquid_properties:unknown', 'No liquid data for "%s".', name);
end
end
