function pho_imag = imag_pho_gen(scene, varargin)

    %%%%%%%%%%%%%%%%% Paser %%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultseed = 1;
    defaultn_exp = 1;
    defaultrat_imag = 0.05;
    defaultn_imag_mu = 500;
    defaultsig = [1,1];
    defaultrange_x = max(defaultsig)*4;
    defaultper_eps = 0;
    
    
    p = inputParser;
   
    addRequired(p,'scene')
    
    addOptional(p,'seed',defaultseed);
    addOptional(p,'n_exp',defaultn_exp);
    addOptional(p,'rat_imag',defaultrat_imag);
    addOptional(p,'n_imag_mu',defaultn_imag_mu);
    addOptional(p,'sig',defaultsig);
    addOptional(p,'range_x',defaultrange_x);
    addOptional(p,'per_eps',defaultper_eps);
   
    parse(p, scene, varargin{:});
   
    scene = p.Results.scene; 
    n_exp = p.Results.n_exp;
    seed = p.Results.seed;
    rat_imag = p.Results.rat_imag;
    n_imag_mu = p.Results.n_imag_mu;
    sig = p.Results.sig;
    range_x = p.Results.range_x;
    per_eps = p.Results.per_eps;

    %%%%%%%%%%%%%%%%%%%%%%% Generate # pho_imagtons
    
    rng(seed);
    
    m = scene(:,1);
    
    tp0 = scene(:,2:3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display range(Field of view)

    x = linspace(-range_x,range_x,1001);
    dx = x(2) - x(1);
    [X,Y] = meshgrid(x,x);
    
    n_imag = binornd(n_imag_mu/rat_imag, rat_imag); % true number of photons
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Generate intensity ~ probability
    I = 0;

    for i = 1:size(tp0,1)

        I = I + m(i)*( exp( -(X-tp0(i,1)).^2 -(Y-tp0(i,2)).^2  )...
                    * sqrt( sig(1)*sig(2)*2/pi ) ).^2*dx^2;  
           
    end
    
    if per_eps ~=0
       
        I = (1-per_eps)*I + per_eps*( exp( -X.^2 - Y.^2  )...
                                   * sqrt( sig(1)*sig(2)*2/pi ) ).^2*dx^2; 
        
    end

    pho_imag = zeros(2,n_imag);

        
    for i = 1:n_imag
    
        [pho_imag(1,i),pho_imag(2,i)] = pinky(x,x,I); % generate arrival photons

    end


end