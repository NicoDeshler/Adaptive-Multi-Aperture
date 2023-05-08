function L_out = SLD_delta_HG(n_modes, scene, varargin)

   % Parser

   defaultmethod = 'all';
   defaultbri_flag = 1;
   defaultpos_var = [];
   defaultmul = 1;
   defaultaperture = [0,0,1];       % aperture configuration [kx,ky,r]
   defaultU = 1;                    % mixing matrix n_ap x n_ap
   
   p = inputParser;

   addRequired(p,'n_modes');
   addRequired(p,'scene');
   
   addOptional(p,'method',defaultmethod);
   addOptional(p,'bri_flag',defaultbri_flag );
   addOptional(p,'pos_var',defaultpos_var );
   addOptional(p,'mul',defaultmul );
   addOptional(p,'aperture',defaultaperture );
   addOptional(p,'U',defaultU );
  
   
   parse(p, n_modes, scene, varargin{:});
   
   n_modes = p.Results.n_modes;
   scene = p.Results.scene;
   
   method = p.Results.method;
   bri_flag = p.Results.bri_flag;
   pos_var = p.Results.pos_var;
   mul = p.Results.mul;
   aperture = p.Results.aperture;
   U = p.Results.U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_modes = n_modes*(n_modes+1)/2;

if strcmp( method, 'HG' )
    
    L_out = eye(N_modes);
    
    return;
    
end

switch method
    
    case 'all'
        
        n_L = size(scene,1);
        
        del_vec = eye(n_L);
        
    case 'difference'
        
        n_L = size(scene,1)-1;
        
        del_vec = del_mat(scene);
        
end

switch method
    
    case {'all','difference'}
        
        if size(scene,2) == 6
            
            scene_temp(:,1) = scene(:,1);
            
            scene_temp(:,2:3) = scene(:,3:4);
            
            scene = scene_temp;
            
        end
        
        rho = rho_HG(n_modes, scene);
        
        [V_rho, D_rho] = eig(rho);
        
        [V_rho, D_rho] = sortem(V_rho, D_rho);
        
        [D1,D2] = meshgrid(diag(D_rho) ,diag(D_rho));
        
        [db, dx, dy] = d_rho_HG(n_modes, scene);
        
        for i = 1:n_L
            
            L(:,:,i) =  V_rho'*sum( repmat( permute( del_vec(:,i), [3,2,1] ), [N_modes,1,1] ).*dx, 3)*V_rho ...
                ./(D1+D2 + (D1+D2 == 0) ).*(D1+D2 ~= 0);
            
            L(:,:,n_L+i) =  V_rho'*sum( repmat( permute( del_vec(:,i), [3,2,1] ), [N_modes,1,1] ).*dy , 3)*V_rho ...
                ./(D1+D2 + (D1+D2 == 0) ).*(D1+D2 ~= 0);
            
        end
        
        if bri_flag
            
            for i = 1:size(scene,1)
                
                L(:,:,end+1) = V_rho'*db(:,:,i)*V_rho./(D1+D2 + (D1+D2 == 0) ).*(D1+D2 ~= 0) ;
                
            end
            
        end
        
        n_now = size(L,3);
        
        for i = 1:n_now
            for j = 1:n_now
                
                H(i,j) = trace( D_rho*( L(:,:,i)*L(:,:,j) ...
                    + L(:,:,j)*L(:,:,i) )/2 ); % tgammas is the QFI matrix elements
                
            end
        end
        
        [V_H, D_H] = eig(H);
        
        [V_H, ~] = sortem(V_H, D_H);
        
        L_out = sum( repmat( permute( V_H(:,1), [3,2,1] ), [N_modes,N_modes,1] ).*L, 3);
        
        %trace(L_out^2*D_rho)
        
    case 'Personick'
        
        % avoid 0 variance
        
        scene(:,5:6) = max(scene(:,5:6), 1e-10);
        
        %rho = rho_HG_Bay(n_modes, scene, [], 0, 'bri_flag', bri_flag);
        rho = rho_HG_Bay_GMM(n_modes, scene, [], 0,...
                            'bri_flag', bri_flag,...
                            'aperture',aperture,...
                            'U',U);
        
        [V_rho, D_rho] = eig(rho);
        
        [V_rho, D_rho] = sortem(V_rho, D_rho);
        
        [D1,D2] = meshgrid(diag(D_rho) ,diag(D_rho));
        
        if bri_flag
            
            n_L = size(scene,1)*3;
            
            para_temp = [scene(:,3); scene(:,4); scene(:,2)/sum(scene(:,2))];
            
        else
            
            n_L = size(scene,1)*2;
            
        end
        
        for i = 1:n_L
            
            %disp('Calculate L')
            temp = zeros(n_L,1);
            temp(i) = 1;
            
            %gamma1 = rho_HG_Bay(n_modes, scene, temp, 1, 'bri_flag', bri_flag);
            gamma1 = rho_HG_Bay_GMM(n_modes, scene, temp, 1, 'bri_flag', bri_flag,'aperture',aperture,'U',U);
            
            L(:,:,i) =  V_rho'*gamma1*V_rho ...
                ./(D1+D2 + (D1+D2 == 0) ).*(D1+D2 ~= 0);
            
        end
        
        n_now = size(L,3);
        
        for i = 1:n_now
            for j = 1:n_now
                
                H(i,j) = trace( D_rho*( L(:,:,i)*L(:,:,j) ...
                    + L(:,:,j)*L(:,:,i) )/2 ); % This is the G matrix used in the Quantum CRLB matrix
                
            end
        end
        
        if bri_flag
            
            Sig_prior(1:size(scene,1)*2, 1:size(scene,1)*2) = [scene(:,3); scene(:,4)]*[scene(:,3); scene(:,4)]' ...
                +diag([scene(:,5); scene(:,6)]);
            
            alp0 = sum(scene(:,2),1);
            alp = scene(:,2)/alp0;
            
            Sig_prior(size(scene,1)*2+1:size(scene,1)*3, size(scene,1)*2+1:size(scene,1)*3) ...
                = diag(alp)/(alp0+1) + alp*alp'*( 1 - 1/(alp0+1) );
            
            Sig_prior(size(scene,1)*2+1:size(scene,1)*3, 1:size(scene,1)*2) = alp*[scene(:,3); scene(:,4)]';
            
            Sig_prior(1:size(scene,1)*2, size(scene,1)*2+1:size(scene,1)*3) ...
                = Sig_prior(size(scene,1)*2+1:size(scene,1)*3, 1:size(scene,1)*2)';
            
            Sigma_Q = Sig_prior - H;
            
        else
            
            Sigma_Q = [scene(:,3); scene(:,4)]*[scene(:,3); scene(:,4)]' ...
                +diag([scene(:,5); scene(:,6)])- H; % This is the Quantum CRLB matrix
            
        end
        
        % Use Sigma_Q
        [V_H, D_H] = eig( Sigma_Q );
        
        [V_H, D_H] = sortem(V_H, D_H);
        
        %Gamma1 = rho_HG_Bay(n_modes, scene, V_H(:,end), 1, 'bri_flag', bri_flag); % this is the density first-moment operator dotted with the projection vector Gamma_1 = ( h_vec dot Gamma_1 vec) which in this case is 
        Gamma1 = rho_HG_Bay_GMM(n_modes, scene, V_H(:,end), 1, 'bri_flag', bri_flag,'aperture',aperture,'U',U); % this is the density first-moment operator dotted with the projection vector Gamma_1 = ( h_vec dot Gamma_1 vec) which in this case is 
        
        L_out =  V_rho'*Gamma1*V_rho ...
            ./(D1+D2 + (D1+D2 == 0) ).*(D1+D2 ~= 0); % this is the joint parameter estimator B_gamma
        
        L_out =  V_rho*L_out*V_rho';  % %%%%%%%%% transforms this back (without this the SLD equation is not satisfied)%%%%%%%%%%% //TODO
        
        
end


end