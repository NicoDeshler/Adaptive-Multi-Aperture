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
ring5 = [Polygon(5,0,'radius',R_max-r(5)),r(5)*ones(5,1)];
ring7 = [Polygon(7,0,'radius',R_max-r(7)),r(7)*ones(7,1)];
golay5 = [Golay5(R_max-r(5)),r(5)*ones(5,1)];
golay7 = [Golay7(R_max-r(7)),r(7)*ones(7,1)];

% collect all the apertures
apertures = {mono,ring3,ring5,ring7,golay5,golay7};
aperture_names = {'Monolith','Ring-3','Ring-5','Ring-7','Golay-5','Golay-7'};


% Datastore variables
DS.A = A;
DS.R_max = R_max;
DS.max_order = 5;
DS.trials = 50;        % trials per configuration
DS.min_sep_frac = 2.^linspace(-6,-3,4);
DS.num_src = 3:6;
DS.num_pho = VisualMagnitude_to_PhotonFlux(10:2:20)*A; 
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