
addpath('Adaptive Algorithm')
DS = DSAdaptive();

t = tiledlayout(2,5);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for a = 1:numel(DS.apertures)
    nexttile;
    hold on
    VisualizeGaussianAperture(DS.apertures{a},DS.R_max)
    hold off
    title(DS.aperture_names{a})
    set(gcf,'color','k');
end
