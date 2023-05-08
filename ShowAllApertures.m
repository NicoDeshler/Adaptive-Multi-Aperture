DS = DSAdaptive();

t = tiledlayout(3,2);

for a= 1:numel(DS.apertures)
    nexttile;
    hold on
    VisualizeGaussianAperture(DS.apertures{a})
    hold off
    title(DS.aperture_names{a})
end
