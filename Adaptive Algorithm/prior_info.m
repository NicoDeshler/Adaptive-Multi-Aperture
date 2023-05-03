function mod_pri = prior_info(models,scene, varargin)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    defaultn_pri = 0;
    defaultpos_pri = 0;
    defaultbri_pri = 0;
    defaultpos_var_pri = 0;
    defaultbri_var_pri = 0;
   
    p = inputParser;
    
    addRequired(p,'models');
    addRequired(p,'scene');
    
    addOptional(p,'n_pri',defaultn_pri);
    addOptional(p,'pos_pri',defaultpos_pri);
    addOptional(p,'bri_pri',defaultbri_pri);
    addOptional(p,'pos_var_pri',defaultpos_var_pri);
    addOptional(p,'bri_var_pri',defaultbri_var_pri);
   
    parse(p, models, scene, varargin{:});
   
    models = p.Results.models;
    scene = p.Results.scene;
    n_pri = p.Results.n_pri;
    pos_pri = p.Results.pos_pri;
    bri_pri = p.Results.bri_pri;
    pos_var_pri = p.Results.pos_var_pri;
    bri_var_pri = p.Results.bri_var_pri;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    mod_pri = models;
    
    n = size(scene,1);
    
    
    if pos_pri || bri_pri
        
        n_pri = 1;
        
    end

    if n_pri
       
        for i = 1:size(mod_pri,3)
            
            if nnz(mod_pri(:,1,i)) == n
            
                mod_pri = models( 1:nnz(models(:,1,i)), :, i);
                
            end
            
        end
        
    end
    
    if bri_pri
       
        mod_pri(:,1) = scene(:,1);
        
    end
    
    if pos_pri
       
        mod_pri(:,3:4) = scene(:,2:3);
        
    end
    
    if bri_var_pri >1
       
        for i = 1:size(mod_pri,3)
            
            mod_pri(:,2,i) = mod_pri(:,2,i)*bri_var_pri;
            
        end
        
    end
    
    if pos_var_pri > 1
       
        for i = 1:size(mod_pri,3)
            
            mod_pri(:,5:6,i) = mod_pri(:,5:6,i)/pos_var_pri;
            
        end
        
    end
        
    

end