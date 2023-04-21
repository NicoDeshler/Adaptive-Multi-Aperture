load SLD_Est_1ap-2src_example1.mat

v = VideoWriter('ap1_2src_example1.avi','Motion JPEG AVI');
v.Quality = 100;
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
    
    xlabel('x [\sigma]')
    ylabel('y [\sigma]')
    title('1-Aperture Adaptive Imaging')
    
    legend({'Scene','Estimate'})
    frame = getframe(gcf);
    writeVideo(v,frame)
    
end

close(v)