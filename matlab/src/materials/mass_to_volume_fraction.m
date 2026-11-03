function phi = mass_to_volume_fraction(w, rho_pcm, rho_eg)
%MASS_TO_VOLUME_FRACTION  Convert filler mass fraction to volume fraction.
phi = (w ./ rho_eg) ./ (w ./ rho_eg + (1 - w) ./ rho_pcm);
end
