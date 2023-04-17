function sence_gen = scene_gen(varargin)

    FWHM = 2.354820;

%%%%%%%%%%%%%%%%% Paser %%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultseed = 1;
    defaultmin_sp = 0.01; % /4 is due to normalization 
    defaultn_s_mu = 8;
	defaults_random = 1;
    defaultn_max = 10;
    defaultn_rl = 2;
    defaultbri_ratio = [];
    defaultshape = 'square';
    
   
    p = inputParser;
   
    addOptional(p,'seed',defaultseed);
    addOptional(p,'min_sp',defaultmin_sp);
    addOptional(p,'n_s_mu',defaultn_s_mu);
	addOptional(p,'s_random',defaults_random);
    addOptional(p,'n_max',defaultn_max);
    addOptional(p,'n_rl',defaultn_rl);
    addOptional(p,'bri_ratio',defaultbri_ratio);
    addOptional(p,'shape',defaultshape);
   
    parse(p, varargin{:});
   
    seed = p.Results.seed;
    min_sp = p.Results.min_sp;
    n_s_mu = p.Results.n_s_mu;
	s_random = p.Results.s_random;
    n_max = p.Results.n_max;
    n_rl = p.Results.n_rl;
    bri_ratio = p.Results.bri_ratio;
    shape = p.Results.shape;
    
    min_sp = FWHM/2*min_sp;
    
%%%%%%%%%%%%%%%%%%%%%%%% Generate # Sources
rng(seed);

if s_random

    n_s_sigma = n_s_mu/10;

while 1
    
  n_s = fix(normrnd(n_s_mu,n_s_sigma));
  
  if n_s > 0 && n_s <= n_max
      
      break;
      
  end
  
end

else

	if n_s_mu <= n_max
		n_s = n_s_mu;
	else
		n_s = n_max;
		warning('Number of input point sources is greater than the N_max. Number of point soureces is set to be N_max')
	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Random distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n_rl = 2; % Field of view

while 1==1

switch shape
    
    case 'square'
        
        tp0 = FWHM/2*(rand(n_s, 2)-1/2)*n_rl;
        
    case 'circle'
   
        randr = FWHM/2*n_rl*rand(n_s,1); %radius
        
        %randr = FWHM/2*n_rl; %radius
        
        randa = rand(n_s,1)*2*pi; %angle
        
        tp0 = repmat( randr, [1,2] ).*[ cos(randa), sin(randa) ];
        
end

% Minimum distance check

dis_check = FWHM*ones(n_s);

for i = 1:n_s
    for j = 1:i
   
        dis_check(i,j) = norm(tp0(i,:)-tp0(j,:));
        
    end
    dis_check(i,i) = FWHM;
end

if any(~(dis_check>min_sp)) == 0

    break;

end

%　Generation of brightness 

end

if isempty(bri_ratio) || bri_ratio < 1

    m = rand(size(tp0,1));
    
else
    
    m = flip( (1:(bri_ratio-1)/(n_s-1):bri_ratio) )';
    
end


m = m/sum(m);

sence_gen = [m, tp0];


end
