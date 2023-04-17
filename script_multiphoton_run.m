clc
addpath('MultiAperture/');
addpath('Parabolic_Cylinder_Functions/')

% load test scene
load('test_scene.mat')
scene = scene(:,:,1);
figure
scatter(scene(:,2,1),scene(:,3,1),'filled','black')
xlabel('x')
ylabel('y')

% ---------------------- aperture --------------------------------- 

D = 30;
R = D/2;
d = 3;
r = d/2;
psf_tilde_sig = r/3; % standard deviation of a single aperture
ap1 = [0,0,R];
ap2 = [Polygon(2,0,'radius',R-r),r/3*ones(2,1)];

%aper_coords = Golay9(R); %aper_coords = [0,0];
%ap_num = size(aper_coords,1);
%aper_rads = psf_tilde_sig*ones(ap_num,1);  %aper_rads = R; 
%aperture = [aper_coords,aper_rads];
aperture = ap1;
ap_num = size(aperture,1);
aper_coords = aperture(:,1:2);
aper_rads = aperture(:,3);


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


% The rayleigh length of the multi-aperture system is defined to be the
% same rayleigh length as a single circular hard aperture with diameter
% D_eff, where D_eff is the diameter of the minimum enclosing circle 
% that contains all of the sub-apertures.
%rl = 2*pi * 1.2197/D_eff; % rayleigh length in units of [rads/length]
rl = 1;


% ----------------- Hermite-Gaussian Local-Aperture Modes----------------
% KWAN WORKS IN COORDINATES WITH DIMENSIONS OF RAYLEIGH LENGTH. SO WE MUST
% DEFINE THE GS BASIS MODES AS FUNCTIONS OF POSITION SPACE COORDINATES WITH
% UNITS OF RAYLEIGH. 
max_order = 3; 
n_HG_modes = max_order + 1;                            % number of 1D HG modes
N_modes = ap_num*(max_order+1)*(max_order+2)/2;     % number of local aperture 2D HG modes
[pj,qj,uj] = Indices_HG_GMMAperture(max_order,ap_num); % vector of linear index map for each mode
U = dftmtx(ap_num);
psf_sig = rl;


save('aperture.mat','aperture','U','n_HG_modes','N_modes','pj','qj','uj','psf_sig')
scene = [.5,0,0;
         .5,0,0];


measurement_multipho(scene(:,:,1), ...
                     'n_max',2, ...
                     'n_pho_group', 10000, ...
                     'n_pho_SLD', 200000, ... % number of photons used throughout bayesian update
                     'n_imag_mu', 10000,...
                     'n_pri',2,...
                     'bri_known',1,...
                     'proj_method','Personick',...
                     'per_eps', 0,...% 0<=per_eps<=1; HG 0 mode background
                     'n_HG_modes', N_modes...
                    ); 

 measurement_multipho(scene(:,:,1), ...
                     'n_max', 4, ...
                     'n_imag_mu', 500000, ...
                     'use_SLD', 0);


