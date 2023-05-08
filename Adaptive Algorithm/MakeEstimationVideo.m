

% set interpreter to latex
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% load in the reconstruction
load 'out/est.mat'

delta_min = min(pdist(scene(:,2:3)));

v = VideoWriter('SLD_Est_2ap_4src_v2','Uncompressed AVI');
v.FrameRate = 5;
open(v)

fig = figure;
for i = 1:numel(cand)
    
    est = cand{i};
    
    scatter(scene(:,2),scene(:,3),'filled','black')
    
    hold on
    scatter(est(:,3),est(:,4),'filled','red')
    hold off
    
    xlim([-0.5,.5])
    ylim([-0.5,.5])
    
    xlabel('$x / \sigma$')
    ylabel('$y / \sigma$')
    title({'2-Aperture Adaptive Imaging',['$\Delta_{min} = \sigma / ',sprintf('%0.3g',1/delta_min),'$']})
    
    legend({'Scene','Estimate'})
    frame = getframe(gcf);
    writeVideo(v,frame)
    
end

close(v)