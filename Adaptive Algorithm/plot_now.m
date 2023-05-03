function plot_now(scene, est)

FWHM = 2.354820;
rl = FWHM/2;
zoom_rat = 8;

range_x = 4;
x = linspace(-range_x,range_x,1001)/rl;

est = sortrows(est, 'descend');
est = est(1:nnz(est(:,1)),:);

axis_vec = [-max(x), max(x), -max(x), max(x)]/zoom_rat;

rel_size = 1000;

%clf
figure
   
scatter( scene(:,2)/rl,scene(:,3)/rl,scene(:,1)*rel_size,'k','filled')
hold on;

scatter( est(:,2)/rl, est(:,3)/rl,est(:,1)*rel_size,'b','LineWidth',2)
hold on;


xlabel('x (rl)');
ylabel('y (rl)');
axis(axis_vec);
box on;
set(gca,'FontSize',14);
legend('Ground Truth', 'Estimation')

drawnow;

end