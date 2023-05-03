function VisualizeGaussianAperture(aperture)
    % demonstrates the gaussian aperture configuration
    % aperture : [kx,ky,sigma] 
    
    n_ap = size(aperture,1);
    kx  = aperture(:,1);
    ky  = aperture(:,2);
    sig = aperture(:,3);
    
    % get max extent of the aperture for plotting
    r = 3*sig;
    R = sqrt(sum(([kx,ky]- mean([kx,ky],1)).^2,2));
    R_eff = max(R+r);
    
    % aperture space over which to plot
    [Kx,Ky] = meshgrid(linspace(-R_eff,R_eff,1001));
    
    % initialize the gaussian aperture
    gauss_aperture = zeros(numel(Kx),1);
    
    % take the max contribution from each aperture
    for s = 1:n_ap
        gauss_s = mvnpdf([Kx(:),Ky(:)],[kx(s),ky(s)],eye(2)*sig(s)^2);
        gauss_aperture = max([gauss_aperture,gauss_s],[],2);
    end
    
    % reshape back to a square
    gauss_aperture = reshape(gauss_aperture,size(Kx));
    
    
    % get the effective Gaussian aperture
    sig_eff = R_eff;
    gauss_eff = reshape(mvnpdf([Kx(:),Ky(:)],mean([kx,ky],1),eye(2)*sig_eff^2),size(Kx));
    %levels = mvnpdf([zeros(5,1),linspace(0,1,5).'*sig_eff],mean([kx,ky],1),eye(2)*sig_eff^2);
    levels = linspace(mvnpdf([0,0],mean([kx,ky],1),eye(2)*sig_eff^2),mvnpdf([0,sqrt(2)*sig_eff],mean([kx,ky],1),eye(2)*sig_eff^2),20);
       
    % plot
    figure
    
    % plot the gaussian multi-aperture with effective gaussian contours and
    % the sigma_eff edge.
    hold on
    imagesc([min(Kx(:)),max(Kx(:))],[min(Ky(:)),max(Ky(:))],gauss_aperture) % plot contour of 
    %contour(Kx,Ky,gauss_eff,levels,'red')                                   % plot contour of effective aperture standard deviation
    circles(mean(kx),mean(ky),sig_eff,'facecolor','none','linewidth',2)
    circles(kx,ky,r,'facecolor','none','edgecolor','green','linewidth',1)
    hold off
    
    % labelling
    axis 'square'
    xlabel('$k_x \, [\tilde{\sigma}]$','interpreter','latex')
    ylabel('$k_y \, [\tilde{\sigma}]$','interpreter','latex')
    title('Gaussian Multi-Aperture')
    xlim([min(Kx(:)),max(Kx(:))])
    ylim([min(Ky(:)),max(Ky(:))])
    legend({'$\tilde{\sigma}_{eff}$','$r = 3\tilde{\sigma}$'},'interpreter','latex')
    
    
end