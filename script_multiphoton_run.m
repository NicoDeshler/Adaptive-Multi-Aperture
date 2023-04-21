clc
addpath('MultiAperture/');
addpath('MultiAperture/circles_v1.1.1');


% load test scene
%{
load('test_scene.mat')
scene = scene(:,:,1);
figure
scatter(scene(:,2,1),scene(:,3,1),'filled','black')
xlabel('x')
ylabel('y')
%}

% ---------------------- aperture --------------------------------- 
sigma = 1;                        % single sub-aperture stdev
r = 3*sigma;                        % max radius of sub-aperture locations

ap1 = [0,0,sigma];
ap2 = [Polygon(2,pi/2,'separation',2*r),sigma*ones(2,1)];
ap3 = [Polygon(3,0,'separation',2*r),sigma*ones(3,1)];
ap4 = [Polygon(4,0,'separation',2*r),sigma*ones(4,1)];
golay9 = [Golay9(6*r),sigma*ones(9,1)];

% set the aperture
aperture = ap3;
ap_num = size(aperture,1);
aper_coords = aperture(:,1:2);
aper_sigs = aperture(:,3);

% get effective aperture standard deviation
if ap_num == 1
    sigma_eff = sigma;
else
    R = sqrt(sum((aper_coords - mean(aper_coords,1)).^2,2));
    sigma_eff = max(R+3*aper_sigs);     % effective gaussian aperture stdev    
end
psf_sig = 1/(2*sigma_eff);          % rl = FWHM of gaussian
rl = (2*sqrt(2*log(2))) * psf_sig;  % gaussian aperture rayleigh length [length units]


% rescale aperture units to be fractional units of the sub-aperture standard deviation.
aperture = aperture / sigma; % this step is critical - the reference unit in Kwan's code is the sigma of an individual sub-aperture



% visualize the aperture
VisualizeGaussianAperture(aperture)

% ---------------------- scene --------------------------------- 
load('test_scene.mat')
scene = scene(:,:,1);
src_coords = scene(:,2:3)...
              * sigma / sigma_eff;
              
num_src = size(src_coords,1);
scene = [ones(num_src,1)/num_src,src_coords];

%{
% get the effective aperture diameter
if ap_num>1
    B = squareform(pdist(aper_coords));             % sub-aperture pairwise centroid distances
    D = aper_rads + aper_rads';                     % sub-aperture pairwise radius sums
    assert(all(triu(B,1) >= triu(D,1),'all'));      % check if any apertures overlap

    % set the effective aperture diameter to the minimum enclosing circle diameter
    cm_coords = aper_coords - mean(aper_coords,1);                          % get centered coordinates (with respect to centroid of centroids -- still need to prove with this works)
    tangent_pts = cm_coords + cm_coords.*aper_rads./vecnorm(cm_coords,2,2); % candidate tangent points where the enclosing circle might touch
    [~,~,R_eff] = SmallestEnclosingCircle(tangent_pts(:,1)',tangent_pts(:,2)'); % effective aperture radius
    D_eff = 2*R_eff;
else
    R_eff = aper_rads(1);
    D_eff = 2*R_eff;                         % set the effective aperture diameter to that of the input aperture
end
%}


% The rayleigh length of the multi-aperture system is defined to be the
% same rayleigh length as a single circular hard aperture with diameter
% D_eff, where D_eff is the diameter of the minimum enclosing circle 
% that contains all of the sub-apertures.
%rl = 2*pi * 1.2197/D_eff; % rayleigh length in units of [rads/length]

% ----------------- Hermite-Gaussian Local-Aperture Modes----------------
% KWAN WORKS IN COORDINATES WITH DIMENSIONS OF RAYLEIGH LENGTH. SO WE MUST
% DEFINE THE GS BASIS MODES AS FUNCTIONS OF POSITION SPACE COORDINATES WITH
% UNITS OF RAYLEIGH. 
max_order = 4; 
n_HG_modes = max_order+1;                                       % number of 1D HG modes
N_modes = ap_num*(max_order)*(max_order+1)/2;           % number of local aperture 2D HG modes
[pj,qj,uj] = Indices_HG_GMMAperture(max_order,ap_num);  % vector of linear index map for each mode
U = dftmtx(ap_num)/sqrt(ap_num);                        % unitary matrix

save('aperture.mat','aperture','U','n_HG_modes','N_modes','pj','qj','uj','psf_sig')


measurement_multipho(scene(:,:,1), ...
                     'n_max',num_src, ...
                     'n_pho_group', 10000, ...
                     'n_pho_SLD', 500000, ... % number of photons used throughout bayesian update
                     'n_imag_mu', 1000,...
                     'n_pri',num_src,...
                     'bri_known',1,...
                     'proj_method','Personick',...
                     'per_eps', 0,...% 0<=per_eps<=1; HG 0 mode background
                     'n_HG_modes', max_order ...
                    ); 
%{
 measurement_multipho(scene(:,:,1), ...
                     'n_max', 4, ...
                     'n_imag_mu', 500000, ...
                     'use_SLD', 0);

%}
