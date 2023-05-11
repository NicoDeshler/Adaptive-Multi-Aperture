A = 7;                  % [m^2] collection area
num_src = 3:6;          % numbers of sources
mx = 10:19;             % visual magnitudes  
pho_flux = VisualMagnitude_to_PhotonFlux(mx); % photon flux for each visual magnitude
T = [1, 10, 30, 60, 300, 1200]; % [s] integration times
T = [1,2,5,10,30, 60, 300, 1200];
%T = 100:100:1200;

%num_src = reshape(num_src, [numel(num_src),1,1]);
%mx = reshape(mx, [1,numel(mx),1]);
%T =  reshape(T, [1,1,numel(T)]);

pho_lwbnd = 1e4;
pho_upbnd = 1e7;

q = 1;

for i = 1:numel(num_src)
   for j = 1:numel(mx)
       for k = 1:numel(T)
           x(q) = num_src(i);
           y(q) = mx(j);
           z(q) = T(k);
           p(q) = A * x(q) * VisualMagnitude_to_PhotonFlux(y(q)) * z(q);
           q = q+1;
       end
   end
end

tiledlayout(1,2)

nexttile
colormap hot
scatter3(x,y,z,50*ones(size(p)),log10(p),'filled')
zlim([0,1e4])
ax = gca;
set(ax,'zscale','log')
xlabel('Number of Sources')
ylabel('Visual Magnitude')
zlabel('Integration Time [s]')
cbar = colorbar;
ylabel(cbar,'$\log_{10}($ Photons $)$','interpreter','latex')
title({'Parameter Space Photon Counts','(No Bounds)'},'Color','w')
axis('square')

set(gca,'Color','none')
set(gca,'YColor','w')
set(gca,'XColor','w')
set(gca,'ZColor','w')
set(cbar,'Color','w')

% thresholded version
nexttile
bnd = pho_upbnd >= p & p>=pho_lwbnd;
scatter3(x(bnd),y(bnd),z(bnd),50*ones(size(p(bnd))),log10(p(bnd)),'filled');
zlim([0,1e4])
ax = gca;
set(ax,'zscale','log')
xlabel('Number of Sources')
ylabel('Visual Magnitude')
zlabel('Integration Time [s]')
cbar = colorbar;
ylabel(cbar,'$\log_{10}($ Photons $)$','interpreter','latex')
title({'Parameter Space Photon Counts',['10^',num2str(log10(pho_lwbnd)),' \leq # Photons \leq', '10^',num2str(log10(pho_upbnd))]},'Color','w')
axis('square')

set(gca,'Color','none')
set(gca,'YColor','w')
set(gca,'XColor','w')
set(gca,'ZColor','w')
set(cbar,'Color','w')

set(gcf,'Color','none')
set(gcf, 'InvertHardcopy', 'off');
set(gcf,'Renderer','painter');
