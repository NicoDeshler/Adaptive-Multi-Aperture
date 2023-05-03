function rho_HG = rho_HG_Bay(n_modes, scene, V, Mom, varargin)


   % Parser
   defaultbri_flag = 0;
   
   p = inputParser;

   addRequired(p,'n_modes');
   addRequired(p,'scene');
   addRequired(p,'V');
   addRequired(p,'Mom');
   
   addOptional(p,'bri_flag',defaultbri_flag );
   
   parse(p, n_modes, scene, V, Mom, varargin{:});
   
   n_modes = p.Results.n_modes;
   scene = p.Results.scene;
   V = p.Results.V;
   Mom = p.Results.Mom;
   
   bri_flag = p.Results.bri_flag;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    scene = sortrows(scene, 'descend');

    scene = scene(1:nnz(scene(:,1)), :);
    
    %N_modes = 3*n_modes*(3*n_modes+1)/2;
    N_modes = n_modes;
    
    if Mom == 1
        
        if bri_flag
            h = [V(1:size(scene,1)), V(size(scene,1)+1:2*size(scene,1))];
            b = V( 2*size(scene,1)+1:3*size(scene,1) );
            
        else
            h = [V(1:size(scene,1)), V(size(scene,1)+1:2*size(scene,1))];
        end
            
    end
    
    alp0 = sum(scene(:,2),1);
    alp = scene(:,2)/alp0;
    Alp = diag(alp)/(alp0+1) + alp*alp'*( 1 - 1/(alp0+1) );

    rho_HG = zeros(N_modes*(N_modes+1)/2);
    
    count = 1;

    for i = 1:N_modes
        for j = 1:i
    
            HG_proj(count).ind_x = i-j;
            HG_proj(count).ind_y = i-HG_proj(count).ind_x-1;  
            
            count = count + 1;
            
        end
    end
    
    for i = 1:N_modes*(N_modes+1)/2
        for j = 1:N_modes*(N_modes+1)/2

            for k = 1:size(scene,1)
                
                if Mom == 0
                    
                    if bri_flag
                        
                    % Gamma0 if bri is unknown
                        
                    rho_HG(i,j) = rho_HG(i,j) + alp(k)...
                                *mom_Bay(scene(k,3), scene(k,5), HG_proj(i).ind_x + HG_proj(j).ind_x)...
                                *mom_Bay(scene(k,4), scene(k,6), HG_proj(i).ind_y + HG_proj(j).ind_y)...
                                /sqrt( factorial(HG_proj(i).ind_x)*factorial(HG_proj(j).ind_x) ...
                                      *factorial(HG_proj(i).ind_y)*factorial(HG_proj(j).ind_y));    
                        
                    else
                        
                    % Gamma0 if bri is known
                    
                    rho_HG(i,j) = rho_HG(i,j) + scene(k,1)...
                                *mom_Bay(scene(k,3), scene(k,5), HG_proj(i).ind_x + HG_proj(j).ind_x)...
                                *mom_Bay(scene(k,4), scene(k,6), HG_proj(i).ind_y + HG_proj(j).ind_y)...
                                /sqrt( factorial(HG_proj(i).ind_x)*factorial(HG_proj(j).ind_x) ...
                                      *factorial(HG_proj(i).ind_y)*factorial(HG_proj(j).ind_y));
                    
                    % GMM rho
                                  
                    end
                               
                elseif Mom == 1

                    for t = 1:size(scene,1)
                    
                        if t == k
                        
                        rho_HG(i,j) = rho_HG(i,j) + alp(k)...
                                    *( h(t,1)*mom_Bay(scene(k,3), scene(k,5), HG_proj(i).ind_x + HG_proj(j).ind_x + 1)...
                                      *mom_Bay(scene(k,4), scene(k,6), HG_proj(i).ind_y + HG_proj(j).ind_y)...
                                      +h(t,2)*mom_Bay(scene(k,3), scene(k,5), HG_proj(i).ind_x + HG_proj(j).ind_x)...
                                      *mom_Bay(scene(k,4), scene(k,6), HG_proj(i).ind_y + HG_proj(j).ind_y + 1) )...
                                    /sqrt( factorial(HG_proj(i).ind_x)*factorial(HG_proj(j).ind_x) ...
                                          *factorial(HG_proj(i).ind_y)*factorial(HG_proj(j).ind_y));
                               
                        else
                            
                            rho_HG(i,j) = rho_HG(i,j) + alp(k)...
                                    *mom_Bay(scene(k,3), scene(k,5), HG_proj(i).ind_x + HG_proj(j).ind_x)...
                                    *mom_Bay(scene(k,4), scene(k,6), HG_proj(i).ind_y + HG_proj(j).ind_y)...
                                    *( h(t,1)*scene(t,3) + h(t,2)*scene(t,4) )...
                                    /sqrt( factorial(HG_proj(i).ind_x)*factorial(HG_proj(j).ind_x) ...
                                          *factorial(HG_proj(i).ind_y)*factorial(HG_proj(j).ind_y));
                            
                        end
                        
                               
                    end
                        %%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                    if bri_flag
  
                    % Gamma1 if bri is unknown
                    
                    rho_HG(i,j) = rho_HG(i,j)...
                                    +mom_Bay(scene(k,3), scene(k,5), HG_proj(i).ind_x + HG_proj(j).ind_x)...
                                    *mom_Bay(scene(k,4), scene(k,6), HG_proj(i).ind_y + HG_proj(j).ind_y)...
                                    *sum( b.*Alp(:,k), 1)...
                                    /sqrt( factorial(HG_proj(i).ind_x)*factorial(HG_proj(j).ind_x) ...
                                          *factorial(HG_proj(i).ind_y)*factorial(HG_proj(j).ind_y));
 
                    end
                    
                    % Gamma1 if bri is known --> Let go
                                   
                end
            
            end
        
        end
    end


end