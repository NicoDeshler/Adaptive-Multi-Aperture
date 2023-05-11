function DS = DSAdaptive()


%%%%%%%%%%  APERTURES %%%%%%%%%%%%%%

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

% collect all the apertures
apertures = {ring3,ring4,ring6,ring7,ring9,mono,golay4,golay6,golay7,golay9};
aperture_names = {'Ring-3','Ring-4','Ring-6','Ring-7','Ring-9','Monolith','Golay-4','Golay-6','Golay-7','Golay-9'};
%apertures = {ring3,mono,ring4,golay4,ring6,golay6,ring7,golay7,ring9,golay9};
%aperture_names = {'Ring-3','Monolith','Ring-4','Golay-4','Ring-6','Golay-6','Ring-7','Golay-7','Ring-9','Golay-9'};


% Datastore variables
DS.A = A;
DS.R_max = R_max;
DS.max_order = 4;
DS.trials = 50;        % trials per configuration
DS.min_sep_frac = 2.^linspace(-6,-3,4);
DS.num_src = 3:6;
DS.num_pho = VisualMagnitude_to_PhotonFlux(fliplr(10:3:19))*A; 
DS.apertures = apertures;
DS.aperture_names = aperture_names;
DS.save_dir = 'out';
DS.cfg_size = [numel(DS.apertures),numel(DS.num_src),numel(DS.min_sep_frac),numel(DS.num_pho)]; % the dimensionality of the parameter space range
DS.data = cell(DS.cfg_size);

end

%{
% rayleigh lengths and gaussian standard devs
    sigma = zeros(1,numel(DS.apertures));
    sigma_eff = zeros(1,numel(DS.apertures));


    for a = 1:numel(DS.apertures)
        ap = DS.apertures{a};
        n_ap = size(ap,1);
        aper_coords = ap(:,1:2);
        aper_rads = ap(:,3);

        % get the effective radius of the multi-aperture system by the smallest
        % enclosing circle
        if n_ap>1
            B = squareform(pdist(aper_coords));             % sub-aperture pairwise centroid distances
            D = aper_rads + aper_rads';                     % sub-aperture pairwise radius sums
            assert(all(triu(B,1) >= triu(D,1),'all'));      % check if any apertures overlap

            % set the effective aperture diameter to the minimum enclosing circle diameter
            cm_coords = aper_coords - mean(aper_coords,1);                          % get centered coordinates (with respect to centroid of centroids -- still need to prove with this works)
            tangent_pts = cm_coords + cm_coords.*aper_rads./vecnorm(cm_coords,2,2); % candidate tangent points where the enclosing circle might touch
            [~,~,R_eff] = SmallestEnclosingCircle(tangent_pts(:,1)',tangent_pts(:,2)'); % effective aperture radius
        else
            R_eff = aper_rads(1);
        end

        % Gaussian widths
        tilde_sigma = r(n_ap)/3;                % stdev of gaussian sub-aperture
        %tilde_sigma_eff = R_eff/3;             % effective stdev of effective gaussian aperture
        tilde_sigma_eff = R_max/3;              % %%%UPDATE%% - use the max radius in order for the scene generation to generate the same scale scenes
        sigma(a) = 1/2/tilde_sigma;             % stdev of sub-aperture gaussian PSF
        sigma_eff(a) = 1/2/tilde_sigma_eff;     % effective stdev of effective gaussian PSF
    end



    % rescale dimensions 
    %(Kwan's algorithm assumes that all units are in terms of the stdev of the sub-aperture)
    for a = 1:numel(apertures)
        apertures{a} = ap .* sigma(a);
    end
%}