clc
addpath('circles_v1.1.1')

% ---------------------- aperture --------------------------------- 

% constraint parameters
A = 7;             % area budget alotted for each aperture configuration [length^2]
D_max = 20;        % multi-aperture max diameter [length]
R_max = D_max/2;   % multi-aperture max radius   [length]

% small function handles
r = @(n) sqrt(A/n/pi); % sub-apeture radius          [length]         
tilde_sigma = @(n) r(n)/3;  % sub-aperture standard deviation for gaussian sub-apertures [length]
 
% the apertures
mono = [0,0,r(1)];
ring3 = [Polygon(3,0,'radius',R_max-r(3)),r(3)*ones(3,1)];
ring4 = [Polygon(4,0,'radius',R_max-r(4)),r(4)*ones(4,1)];
ring6 = [Polygon(6,0,'radius',R_max-r(6)),r(6)*ones(6,1)];
ring7 = [Polygon(7,0,'radius',R_max-r(7)),r(7)*ones(7,1)];
ring9 = [Polygon(9,0,'radius',R_max-r(9)),r(9)*ones(9,1)];
golay4 = [Golay4(R_max-r(4)),r(4)*ones(4,1)];
golay6 = [Golay6(R_max-r(6)),r(6)*ones(6,1)];
golay7 = [Golay7(R_max-r(7)),r(7)*ones(7,1)];
golay9 = [Golay9(R_max-r(9)),r(9)*ones(9,1)];

% set the aperture
aperture = ring3;
ap_num = size(aperture,1);
aper_coords = aperture(:,1:2);
aper_rads = aperture(:,3);

%------------ Rescale dimensions ---------------- %
    
% Gaussian widths
tilde_sigma = min(aperture(:,3))/3;  % stdev of gaussian sub-aperture
tilde_sigma_eff = R_max/3;        % stdev of the effective gaussian aperture
sigma = 1/2/tilde_sigma;             % stdev of sub-aperture gaussian PSF
sigma_eff = 1/2/tilde_sigma_eff;     % effective stdev of effective gaussian PSF

% make units of the aperture coordinates equal to [(2 sigma)^-1] to
% correspond to the fact that Kwan defines the image-space coordinates
% to be [x/(2 sigma)].
aperture = aperture ./ tilde_sigma;   

% visualize the aperture
%VisualizeGaussianAperture(aperture)

% ---------------------- scene --------------------------------- 
% make the min_sep_frac of the effective rayleigh length scaled in
% terms of the minimum separation fraction of the sub-aperture
% rayleigh length
min_sep_frac = 1/4;
min_sep_frac = min_sep_frac * sigma_eff / sigma;
centroid_aligned = 1;

src_brites = [1e5;1];
src_brites = src_brites / sum(src_brites);
src_coords = genMinDistConstellation(src_brites, min_sep_frac,centroid_aligned);
src_coords = src_coords - src_coords(1,:);

num_src = size(src_coords,1);
scene = [src_brites, src_coords];


% The rayleigh length of the multi-aperture system is defined to be the
% same rayleigh length as a single circular hard aperture with diameter
% D_eff, where D_eff is the diameter of the minimum enclosing circle 
% that contains all of the sub-apertures.

% ----------------- Hermite-Gaussian Local-Aperture Modes----------------
% KWAN WORKS IN COORDINATES WITH DIMENSIONS OF RAYLEIGH LENGTH. SO WE MUST
% DEFINE THE GS BASIS MODES AS FUNCTIONS OF POSITION SPACE COORDINATES WITH
% UNITS OF RAYLEIGH. 
max_order = 4; 
n_HG_modes = max_order+1;                                       % number of 1D HG modes
N_modes = ap_num*(max_order)*(max_order+1)/2;           % number of local aperture 2D HG modes
[pj,qj,uj] = Indices_HG_GMMAperture(max_order,ap_num);  % vector of linear index map for each mode
U = dftmtx(ap_num)/sqrt(ap_num);                        % unitary matrix

% ----------------- Photons ----------------
n_pho_imag_mu = 1000;
n_pho_group = 2e4;
n_pho_SLD = 5e5;


% ----------------- Run Estimation ----------------
[est, est_trace] = measurement_multipho(scene, ...
                     'n_max',num_src, ...
                     'n_pho_group', n_pho_group, ...
                     'n_pho_SLD', n_pho_SLD, ...
                     'n_imag_mu', n_pho_imag_mu,...
                     'n_pri',num_src,...
                     'bri_known',1,...
                     'proj_method','Personick',...
                     'per_eps', 0,...% 0<=per_eps<=1; HG 0 mode background
                     'n_HG_modes', max_order, ...
                     'aperture',aperture,...
                     'U',U,...
                     'visualize',1, ...
                     'rel_sigma',sigma/sigma_eff ...
                    ); 
