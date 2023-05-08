 function prob_gen = prob_gen_2nd_mom(scene, L, varargin)

format short
 
% Parser

   defaultn_modes = 2;
   defaultHG_flag = 0;
   defaultper_eps = 0;
   defaultaperture = [0,0,1];
   defaultU = 1;
   
   p = inputParser;

   addRequired(p,'scene');
   addRequired(p,'L');
   
   addOptional(p,'n_modes',defaultn_modes);
   addOptional(p,'HG_flag',defaultHG_flag);
   addOptional(p,'per_eps',defaultper_eps);
   addOptional(p,'aperture',defaultaperture);
   addOptional(p,'U',defaultU);
   
   parse(p, scene, L, varargin{:});
   
   scene = p.Results.scene;
   L = p.Results.L;
   
   n_modes = p.Results.n_modes;
   HG_flag = p.Results.HG_flag;
   per_eps = p.Results.per_eps;
   aperture = p.Results.aperture;
   U = p.Results.U;
   %%%%%%%%%%%%%%%%%%%%%%%%%%
   
%   scene_HG = rho_HG(n_modes, scene);
   scene_HG = rho_HG_GMM(n_modes, scene, aperture, U);
   
   if per_eps~=0
       
       temp = zeros(size(scene_HG));
       
       temp(1,1) = 1;
      
       scene_HG = (1-per_eps)*scene_HG + per_eps*temp(1,1);
       
   end
   
   if ~HG_flag
   
       [V_L, D_L] = eig(L);
   
       D_L = real(D_L);
   
       [V_L, D_L] = sortem(V_L, D_L);
          
       for i = 1:size(L,2)
      
           prob_gen(i) = (V_L(:,i))'*scene_HG*V_L(:,i);
       
       end
   
       prob_gen = real(prob_gen);
   
   else
       
       disp('Use HG modes')
       prob_gen = diag(scene_HG);
       
   end
   
   prob_gen = prob_gen.*(prob_gen >= 0);
   
   prob_gen = [max(0, 1 - sum(prob_gen)), prob_gen];


end