function Pri_up = pri_up_3N(pho_get, models, mod_pick, dtp, dm)

global n_nax;

dim = 21;
n = size(dtp,1);

mod_pick = mod_pick(1:nnz(mod_pick(:,1)),:);

dtp = dtp(1:nnz(mod_pick(:,1)),:);

dm = dm(1:nnz(mod_pick(:,1)),:);

x = mod_pick(:,3);
y = mod_pick(:,4);

[Xq, Xr] = meshgrid(x, x);
[Yq, Yr] = meshgrid(y, y);

dx = dtp(:,1);
dy = dtp(:,2);

dX = (repmat(dx,[1,size(dtp,1)]))';
dY = (repmat(dy,[1,size(dtp,1)]))';

L = SLD_3N(mod_pick(:,1), mod_pick(:,3:4), dim, 'dm', dm, 'dtp', dtp, 'out', 'SLD = L;');
C = SLD_3N(mod_pick(:,1), mod_pick(:,3:4), dim, 'dm', dm, 'dtp', dtp, 'out', 'SLD = V_rho;');
Aper = SLD_3N(mod_pick(:,1), mod_pick(:,3:4), dim, 'dm', dm, 'dtp', dtp, 'out', 'SLD = Aper;'); % Chnage for multi parameters
Apara = SLD_3N(mod_pick(:,1), mod_pick(:,3:4), dim, 'dm', dm, 'dtp', dtp, 'out', 'SLD = Apara;');
    
Omega = -Aper/Apara;

[V_L, D_L] = eig(L(:,:,1)); % Chnage for multi parameter
[V_L, D_L] = sortem(V_L, D_L);

pri_up = zeros(size(models));

b_del = 1e-6;

%{
l_rho = V_L(1:n,:);
l_mu = V_L(n+1:end,:);
%}

for i = 1:size(models,3)
    
    model_temp = models(:,:,i);
    
    model_temp = model_temp(1:nnz(model_temp(:,1)),1:6);

    % Model: j parameter
    
    b = model_temp(:,1);
    
    a = model_temp(:,2);
    
    xj = model_temp(:,3);
    
    yj = model_temp(:,4);
    
    sigx = model_temp(:,5);
    
    sigy = model_temp(:,6);
    
    Xq_3D = repmat(Xq, [1,1,size(xj,1)]);
    Xr_3D = repmat(Xr, [1,1,size(xj,1)]);
    Yq_3D = repmat(Yq, [1,1,size(xj,1)]);
    Yr_3D = repmat(Yr, [1,1,size(xj,1)]);
    
    % Moments
    
    %{
    check_prob = prob_gen(b, b, mod_pick(:,3:4), dtp, mod_pick(:,3:4))
    
    sum(check_prob)
    %}
    
    m_x_0 = mom(0, x, xj, sigx);
    
    m_y_0 = mom(0, y, yj, sigy);
    
    m_x_1 = mom(1, x, xj, sigx);
    
    m_y_1 = mom(1, y, yj, sigy);
    
    m_x_2 = mom(2, x, xj, sigx);
    
    m_y_2 = mom(2, y, yj, sigy);
    
    m_x_3 = mom(3, x, xj, sigx);
    
    m_y_3 = mom(3, y, yj, sigy);
    
    m_x_4 = mom(4, x, xj, sigx);
    
    m_y_4 = mom(4, y, yj, sigy);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% K0 matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    K0_xx = ( m_x_2 - (Xr_3D + Xq_3D).*m_x_1 + Xr_3D.*Xq_3D.*m_x_0 ).*m_y_0;
    
    K0_xy = ( m_x_1 - Xr_3D.*m_x_0 ).*( m_y_1 - Yq_3D.*m_y_0 );
    
    K0_yx = ( m_x_1 - Xq_3D.*m_x_0 ).*( m_y_1 - Yr_3D.*m_y_0 );
    
    K0_yy = ( m_y_2 - (Yr_3D + Yq_3D).*m_y_1 + Yr_3D.*Yq_3D.*m_y_0 ).*m_x_0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% K1 matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    K1X_xx = ( m_x_3 - (Xr_3D + Xq_3D).*m_x_2 + Xr_3D.*Xq_3D.*m_x_1 ).*m_y_0;
    
    K1X_xy = ( m_x_2 - Xr_3D.*m_x_1 ).*( m_y_1 - Yq_3D.*m_y_0 );
    
    K1X_yx = ( m_x_2 - Xq_3D.*m_x_1 ).*( m_y_1 - Yr_3D.*m_y_0 );
    
    K1X_yy = ( m_y_2 - (Yr_3D + Yq_3D).*m_y_1 + Yr_3D.*Yq_3D.*m_y_0 ).*m_x_1;    
    

    
    K1Y_xx = ( m_x_2 - (Xr_3D + Xq_3D).*m_x_1 + Xr_3D.*Xq_3D.*m_x_0 ).*m_y_1; 
    
    K1Y_xy = ( m_x_1 - Xr_3D.*m_x_0 ).*( m_y_2 - Yq_3D.*m_y_1 );
    
    K1Y_yx = ( m_x_1 - Xq_3D.*m_x_0 ).*( m_y_2 - Yr_3D.*m_y_1 );
    
    K1Y_yy = ( m_y_3 - (Yr_3D + Yq_3D).*m_y_2 + Yr_3D.*Yq_3D.*m_y_1 ).*m_x_0;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j = 1:size(xj,1)
        
    % p(l_i)
        
        for p = 1:size(x,1)
            for q = 1:size(x,1)
                %{
                K0(p,q,j) =  [dx(p),dy(p)] ...
                            *( eye(2)...
                              -[K0_xx(p,q,j), K0_xy(p,q,j); ...
                                K0_yx(p,q,j), K0_yy(p,q,j)] ) ...
                            *[dx(q);dy(q)];
                %}
                
                K0(p,q,j) =  [dx(p),dy(p)] ...
                            *( [K0_xx(p,q,j), K0_xy(p,q,j); ...
                                K0_yx(p,q,j), K0_yy(p,q,j)] ) ...
                            *[dx(q);dy(q)];
            end
        end
    
        p_l_11(:,:,j) = C'*( m_x_0(:,:,j).*m_y_0(:,:,j) )*C;
        
        p_l_12(:,:,j) =  C'*( m_x_0(:,:,j).*m_y_0(:,:,j) )*C*Omega ...
                        +C'*( dX.*(m_x_1(:,:,j) - Xq.*m_x_0(:,:,j)).*m_y_0(:,:,j) ...
                            + dY.*(m_y_1(:,:,j) - Yq.*m_y_0(:,:,j)).*m_x_0(:,:,j) )/Apara;
                        
        p_l_21(:,:,j) = Omega'*C'*( m_x_0(:,:,j).*m_y_0(:,:,j) )'*C ...
                        +(inv(Apara))'*( dX'.*(m_x_1(:,:,j) - Xr.*m_x_0(:,:,j)).*m_y_0(:,:,j) ...
                                       + dY'.*(m_y_1(:,:,j) - Yr.*m_y_0(:,:,j)).*m_x_0(:,:,j) )*C;
        
        p_l_22(:,:,j) =   Omega'*C'*( m_x_0(:,:,j).*m_y_0(:,:,j) )*C*Omega ...
                        + Omega'*C'*( dX.*(m_x_1(:,:,j) - Xq.*m_x_0(:,:,j)).*m_y_0(:,:,j) ...
                                    + dY.*(m_y_1(:,:,j) - Yq.*m_y_0(:,:,j)).*m_x_0(:,:,j) )/Apara ...
                        +(Omega'*C'*( dX.*(m_x_1(:,:,j) - Xq.*m_x_0(:,:,j)).*m_y_0(:,:,j) ...
                                    + dY.*(m_y_1(:,:,j) - Yq.*m_y_0(:,:,j)).*m_x_0(:,:,j) )/Apara)' ...
                        +(inv(Apara))'*K0(:,:,j)/Apara;
        
        % p_l should be positive semidefinite
        
        p_l(:,:,j) = [  p_l_11(:,:,j),    p_l_12(:,:,j);...
                      ( p_l_12(:,:,j) )', p_l_22(:,:,j)];

        
          %{               
        p_l(:,:,j) = [ p_l_11(:,:,j), p_l_12(:,:,j);...
                       p_l_21(:,:,j), p_l_22(:,:,j) ];
         %}         
        p_li(j,i) = ( V_L(:,pho_get) )'*p_l(:,:,j)*V_L(:,pho_get);
                        
    end
    
    p_li(:,i);
    
    p_norm(i) = ( (a-1)/( sum(a) - numel(a) ) )'*p_li(:,i);
    
    pri_up(:,1:2,i) = zero_padding(intensity_up(b, a, p_li(:,i), b_del));
    

    %%%%%%%%%%%%%% int \theta p(l_i|\theta)*p(\theta) d\theta %%%%%%%%%%%%%
    
    for k = 1:size(xj,1)
        
        clear K1_x K1_y
        
        for j = 1:size(xj,1)
            
            if j == k
                
                for p = 1:size(x,1)
                    for q = 1:size(x,1)
            
                        K1_x(p,q,j) =  [dx(p),dy(p)] ...
                                      *( [K1X_xx(p,q,j), K1X_xy(p,q,j); ...
                                          K1X_yx(p,q,j), K1X_yy(p,q,j)] ) ...
                                      *[dx(q);dy(q)];
                                  
                        K1_y(p,q,j) =  [dx(p),dy(p)] ...
                                      *( [K1Y_xx(p,q,j), K1Y_xy(p,q,j); ...
                                          K1Y_yx(p,q,j), K1Y_yy(p,q,j)] ) ...
                                      *[dx(q);dy(q)];
            
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%% x_bar %%%%%%%%%%%%%%%%%%%%%%%%%%
            
                x_bar_11(:,:,j,k) = C'*( m_x_1(:,:,j).*m_y_0(:,:,j) )*C;
                
                x_bar_12(:,:,j,k) = C'*( m_x_1(:,:,j).*m_y_0(:,:,j) )*C*Omega ...
                                   +C'*( dX.*(m_x_2(:,:,j) - Xq.*m_x_1(:,:,j)).*m_y_0(:,:,j) ...
                                       + dY.*(m_y_1(:,:,j) - Yq.*m_y_0(:,:,j)).*m_x_1(:,:,j) )/Apara;
                
                x_bar_22(:,:,j,k) =  Omega'*C'*( m_x_1(:,:,j).*m_y_0(:,:,j ) )*C*Omega ...
                                   + Omega'*C'*( dX.*(m_x_2(:,:,j) - Xq.*m_x_1(:,:,j)).*m_y_0(:,:,j) ...
                                               + dY.*(m_y_1(:,:,j) - Yq.*m_y_0(:,:,j)).*m_x_1(:,:,j) )/Apara ...
                                   +(Omega'*C'*( dX.*(m_x_2(:,:,j) - Xq.*m_x_1(:,:,j)).*m_y_0(:,:,j) ...
                                               + dY.*(m_y_1(:,:,j) - Yq.*m_y_0(:,:,j)).*m_x_1(:,:,j) )/Apara)' ...
                                   +(inv(Apara))'*K1_x(:,:,j)/Apara;
                               
                x_bar(:,:,j,k) = [  x_bar_11(:,:,j,k),    x_bar_12(:,:,j,k);...
                                  ( x_bar_12(:,:,j,k) )', x_bar_22(:,:,j,k)];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%% y_bar %%%%%%%%%%%%%%%%%%%%%%%%%%
                
                y_bar_11(:,:,j,k) = C'*( m_x_0(:,:,j).*m_y_1(:,:,j) )*C;
                
                y_bar_12(:,:,j,k) = C'*( m_x_0(:,:,j).*m_y_1(:,:,j) )*C*Omega ...
                                   +C'*( dX.*(m_x_1(:,:,j) - Xq.*m_x_0(:,:,j)).*m_y_1(:,:,j) ...
                                       + dY.*(m_y_2(:,:,j) - Yq.*m_y_1(:,:,j)).*m_x_0(:,:,j) )/Apara;
                
                y_bar_22(:,:,j,k) =  Omega'*C'*( m_x_0(:,:,j).*m_y_1(:,:,j ) )*C*Omega ...
                                   + Omega'*C'*( dX.*(m_x_1(:,:,j) - Xq.*m_x_0(:,:,j)).*m_y_1(:,:,j) ...
                                               + dY.*(m_y_2(:,:,j) - Yq.*m_y_1(:,:,j)).*m_x_0(:,:,j) )/Apara ...
                                   +(Omega'*C'*( dX.*(m_x_1(:,:,j) - Xq.*m_x_0(:,:,j)).*m_y_1(:,:,j) ...
                                               + dY.*(m_y_2(:,:,j) - Yq.*m_y_1(:,:,j)).*m_x_0(:,:,j) )/Apara)' ...
                                   +(inv(Apara))'*K1_y(:,:,j)/Apara;
                               
                y_bar(:,:,j,k) = [  y_bar_11(:,:,j,k),    y_bar_12(:,:,j,k);...
                                  ( y_bar_12(:,:,j,k) )', y_bar_22(:,:,j,k)];
                               
            
            else
                
                x_bar(:,:,j,k) = xj(k)*p_l(:,:,j);
                
                y_bar(:,:,j,k) = yj(k)*p_l(:,:,j);
                
            end
            
            pri_up(k,3,i) = pri_up(k,3,i)+ ( (a(j)-1)/( sum(a) - numel(a) ) ) ...
                                          *( V_L(:,pho_get) )'*x_bar(:,:,j,k)*V_L(:,pho_get)/p_norm(i);
            
            pri_up(k,4,i) = pri_up(k,4,i)+ ( (a(j)-1)/( sum(a) - numel(a) ) ) ...
                                          *( V_L(:,pho_get) )'*y_bar(:,:,j,k)*V_L(:,pho_get)/p_norm(i);
        
        end
        
    end
    
    %%%%% int ( \theta-bar{\theta} )^2 p(l_i|\theta)*p(\theta) d\theta %%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% K2 matriices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X_bar = repmat( permute( nonzeros(pri_up(:,3,i)), [3,2,1] ), [size(x,1),size(x,1),1] );
    
    Y_bar = repmat( permute( nonzeros(pri_up(:,4,i)), [3,2,1] ), [size(x,1),size(x,1),1] );
    
    K2X_xx = ( m_x_4 - (Xr_3D + Xq_3D + 2*X_bar).*m_x_3 ...
           + (Xr_3D.*Xq_3D + 2*Xr_3D.*X_bar + 2*Xq_3D.*X_bar + X_bar.^2).*m_x_2 ...
           - X_bar.*(2*Xr_3D.*Xq_3D + Xr_3D.*X_bar + Xq_3D.*X_bar).*m_x_1 ...
           + X_bar.^2.*Xr_3D.*Xq_3D.*m_x_0 ).*m_y_0;
    
    K2X_xy = ( m_x_3 - (Xr_3D + 2*X_bar).*m_x_2 + X_bar.*(2*Xr_3D + X_bar).*m_x_1 ...
             - X_bar.^2.*Xr_3D.*m_x_0 ).*( m_y_1 - Yq_3D.*m_y_0 );
    
    K2X_yx = ( m_x_3 - (Xq_3D + 2*X_bar).*m_x_2 + X_bar.*(2*Xq_3D + X_bar).*m_x_1 ...
             - X_bar.^2.*Xq_3D.*m_x_0 ).*( m_y_1 - Yr_3D.*m_y_0 );
    
    K2X_yy = ( m_x_2 - 2*X_bar.*m_x_1 + X_bar.^2.*m_x_0 ) ...
           .*( m_y_2 - (Yr_3D + Yq_3D).*m_y_1 + Yr_3D.*Yq_3D.*m_y_0 ) ;   
       
       
   
    K2Y_xx = ( m_y_2 - 2*Y_bar.*m_y_1 + Y_bar.^2.*m_y_0 ) ...
           .*( m_x_2 - (Xr_3D + Xq_3D).*m_x_1 + Xr_3D.*Xq_3D.*m_x_0 ) ;
    
    K2Y_xy = ( m_y_3 - (Yq_3D + 2*Y_bar).*m_y_2 + Y_bar.*(2*Yq_3D + Y_bar).*m_y_1 ...
             - Y_bar.^2.*Yq_3D.*m_y_0 ).*( m_x_1 - Xr_3D.*m_x_0 );
    
    K2Y_yx = ( m_y_3 - (Yr_3D + 2*Y_bar).*m_y_2 + Y_bar.*(2*Yr_3D + Y_bar).*m_y_1 ...
             - Y_bar.^2.*Yr_3D.*m_y_0 ).*( m_x_1 - Xq_3D.*m_x_0 );
    
    K2Y_yy = ( m_y_4 - (Yr_3D + Yq_3D + 2*Y_bar).*m_y_3 ...
           + (Yr_3D.*Yq_3D + 2*Yr_3D.*Y_bar + 2*Yq_3D.*Y_bar + Y_bar.^2).*m_y_2 ...
           - Y_bar.*(2*Yr_3D.*Yq_3D + Yr_3D.*Y_bar + Yq_3D.*Y_bar).*m_y_1 ...
           + Y_bar.^2.*Yr_3D.*Yq_3D.*m_y_0 ).*m_x_0;
    
    for k = 1:size(xj,1)
        
        clear K2_x K2_y
        
        for j = 1:size(xj,1)
            
            if j == k
                
                for p = 1:size(x,1)
                    for q = 1:size(x,1)
            
                        K2_x(p,q,j) =  [dx(p),dy(p)] ...
                                      *( [K2X_xx(p,q,j), K2X_xy(p,q,j); ...
                                          K2X_yx(p,q,j), K2X_yy(p,q,j)] ) ...
                                      *[dx(q);dy(q)];
                                  
                        K2_y(p,q,j) =  [dx(p),dy(p)] ...
                                      *( [K2Y_xx(p,q,j), K2Y_xy(p,q,j); ...
                                          K2Y_yx(p,q,j), K2Y_yy(p,q,j)] ) ...
                                      *[dx(q);dy(q)];
            
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%% sigx_bar %%%%%%%%%%%%%%%%%%%%%%%%%%
            
                sigx_bar_11(:,:,j,k) = C'*( (m_x_2(:,:,j) - 2*X_bar(:,:,j).*m_x_1(:,:,j) ...
                                           + X_bar(:,:,j).^2.*m_x_0(:,:,j)).*m_y_0(:,:,j) )*C;
                
                sigx_bar_12(:,:,j,k) = C'*( (m_x_2(:,:,j) - 2*X_bar(:,:,j).*m_x_1(:,:,j) ...
                                           + X_bar(:,:,j).^2.*m_x_0(:,:,j)).*m_y_0(:,:,j) )*C*Omega ...
                                      +C'*( dX.*(m_x_3(:,:,j) - (Xq + 2*X_bar(:,:,j)).*m_x_2(:,:,j) ...
                                               + X_bar(:,:,j).*(2*Xq + X_bar(:,:,j)).*m_x_1(:,:,j) ...
                                               - X_bar(:,:,j).^2.*Xq.*m_x_0(:,:,j) ).*m_y_0(:,:,j) ...
                                          + dY.*(m_y_1(:,:,j) - Yq.*m_y_0(:,:,j))...
                                              .*(m_x_2(:,:,j) - 2*X_bar(:,:,j).*m_x_1(:,:,j) ...
                                               + X_bar(:,:,j).^2.*m_x_0(:,:,j)) )/Apara;
                
                sigx_bar_22(:,:,j,k) =  Omega'*C'*( (m_x_2(:,:,j) - 2*X_bar(:,:,j).*m_x_1(:,:,j) ...
                                                   + X_bar(:,:,j).^2.*m_x_0(:,:,j)).*m_y_0(:,:,j) )*C*Omega ...
                                      + Omega'*C'*( dX.*(m_x_3(:,:,j) - (Xq + 2*X_bar(:,:,j)).*m_x_2(:,:,j) ...
                                                       + X_bar(:,:,j).*(2*Xq + X_bar(:,:,j)).*m_x_1(:,:,j) ...
                                                       - X_bar(:,:,j).^2.*Xq.*m_x_0(:,:,j) ).*m_y_0(:,:,j) ...
                                                  + dY.*(m_y_1(:,:,j) - Yq.*m_y_0(:,:,j))...
                                                      .*(m_x_2(:,:,j) - 2*X_bar(:,:,j).*m_x_1(:,:,j) ...
                                                       + X_bar(:,:,j).^2.*m_x_0(:,:,j)) )/Apara ...
                                      +(Omega'*C'*( dX.*(m_x_3(:,:,j) - (Xq + 2*X_bar(:,:,j)).*m_x_2(:,:,j) ...
                                                       + X_bar(:,:,j).*(2*Xq + X_bar(:,:,j)).*m_x_1(:,:,j) ...
                                                       - X_bar(:,:,j).^2.*Xq.*m_x_0(:,:,j) ).*m_y_0(:,:,j) ...
                                                  + dY.*(m_y_1(:,:,j) - Yq.*m_y_0(:,:,j))...
                                                      .*(m_x_2(:,:,j) - 2*X_bar(:,:,j).*m_x_1(:,:,j) ...
                                                       + X_bar(:,:,j).^2.*m_x_0(:,:,j)) )/Apara)' ...
                                      +(inv(Apara))'*K2_x(:,:,j)/Apara;
                               
                sigx_bar(:,:,j,k) = [  sigx_bar_11(:,:,j,k),    sigx_bar_12(:,:,j,k);...
                                     ( sigx_bar_12(:,:,j,k) )', sigx_bar_22(:,:,j,k)];
                
                %{
                                     for checking
                                     
                check_11(j) = ( l_rho(:,pho_get) )'*sigx_bar_11(:,:,j,k)*l_rho(:,pho_get)/p_norm(i); 
            
                check_12(j) = ( l_rho(:,pho_get) )'*sigx_bar_12(:,:,j,k)*l_mu(:,pho_get)/p_norm(i);
            
                check_22(j) = ( l_mu(:,pho_get) )'*sigx_bar_22(:,:,j,k)*l_mu(:,pho_get)/p_norm(i);
                
                check_sum(j) = check_11(j) + 2*check_12(j) + check_22(j);
                %}
                
                %%%%%%%%%%%%%%%%%%%%%%%%%% sigy_bar %%%%%%%%%%%%%%%%%%%%%%%%%%
                
                sigy_bar_11(:,:,j,k) = C'*( (m_y_2(:,:,j) - 2*Y_bar(:,:,j).*m_y_1(:,:,j) ...
                                           + Y_bar(:,:,j).^2.*m_y_0(:,:,j)).*m_x_0(:,:,j) )*C;
                
                sigy_bar_12(:,:,j,k) = C'*( (m_y_2(:,:,j) - 2*Y_bar(:,:,j).*m_y_1(:,:,j) ...
                                           + Y_bar(:,:,j).^2.*m_y_0(:,:,j)).*m_x_0(:,:,j) )*C*Omega ...
                                      +C'*( dY.*(m_y_3(:,:,j) - (Yq + 2*Y_bar(:,:,j)).*m_y_2(:,:,j) ...
                                                + Y_bar(:,:,j).*(2*Yq + Y_bar(:,:,j)).*m_y_1(:,:,j) ...
                                                - Y_bar(:,:,j).^2.*Yq.*m_y_0(:,:,j) ).*m_x_0(:,:,j) ...
                                          + dX.*(m_x_1(:,:,j) - Xq.*m_x_0(:,:,j))...
                                              .*(m_y_2(:,:,j) - 2*Y_bar(:,:,j).*m_y_1(:,:,j) ...
                                               + Y_bar(:,:,j).^2.*m_y_0(:,:,j)) )/Apara;         
                
                sigy_bar_22(:,:,j,k) =  Omega'*C'*( (m_y_2(:,:,j) - 2*Y_bar(:,:,j).*m_y_1(:,:,j) ...
                                                   + Y_bar(:,:,j).^2.*m_y_0(:,:,j)).*m_x_0(:,:,j) )*C*Omega ...
                                      + Omega'*C'*( dY.*(m_y_3(:,:,j) - (Yq + 2*Y_bar(:,:,j)).*m_y_2(:,:,j) ...
                                                       + Y_bar(:,:,j).*(2*Yq + Y_bar(:,:,j)).*m_y_1(:,:,j) ...
                                                       - Y_bar(:,:,j).^2.*Yq.*m_y_0(:,:,j) ).*m_x_0(:,:,j) ...
                                                  + dX.*(m_x_1(:,:,j) - Xq.*m_x_0(:,:,j))...
                                                      .*(m_y_2(:,:,j) - 2*Y_bar(:,:,j).*m_y_1(:,:,j) ...
                                                       + Y_bar(:,:,j).^2.*m_y_0(:,:,j)) )/Apara ...
                                      +(Omega'*C'*( dY.*(m_y_3(:,:,j) - (Yq + 2*Y_bar(:,:,j)).*m_y_2(:,:,j) ...
                                                       + Y_bar(:,:,j).*(2*Yq + Y_bar(:,:,j)).*m_y_1(:,:,j) ...
                                                       - Y_bar(:,:,j).^2.*Yq.*m_y_0(:,:,j) ).*m_x_0(:,:,j) ...
                                                  + dX.*(m_x_1(:,:,j) - Xq.*m_x_0(:,:,j))...
                                                      .*(m_y_2(:,:,j) - 2*Y_bar(:,:,j).*m_y_1(:,:,j) ...
                                                       + Y_bar(:,:,j).^2.*m_y_0(:,:,j)) )/Apara)' ...
                                      +(inv(Apara))'*K2_y(:,:,j)/Apara;
                               
                sigy_bar(:,:,j,k) = [  sigy_bar_11(:,:,j,k),    sigy_bar_12(:,:,j,k);...
                                     ( sigy_bar_12(:,:,j,k) )', sigy_bar_22(:,:,j,k)];
                               
            
            else
                
                sigx_bar(:,:,j,k) = ( sigx(k) + (pri_up(k,3,i) - xj(k))^2 )*p_l(:,:,j);
                
                sigy_bar(:,:,j,k) = ( sigy(k) + (pri_up(k,4,i) - yj(k))^2 )*p_l(:,:,j);
                
            end
            
            
            
            check(j,k) = ( V_L(:,pho_get) )'*sigx_bar(:,:,j,k)*V_L(:,pho_get)/p_norm(i);
            
            pri_up(k,5,i) = pri_up(k,5,i)+ ( (a(j)-1)/( sum(a) - numel(a) ) ) ...
                                          *( V_L(:,pho_get) )'*sigx_bar(:,:,j,k)*V_L(:,pho_get)/p_norm(i);
            
            pri_up(k,6,i) = pri_up(k,6,i)+ ( (a(j)-1)/( sum(a) - numel(a) ) ) ...
                                          *( V_L(:,pho_get) )'*sigy_bar(:,:,j,k)*V_L(:,pho_get)/p_norm(i);
        
        end
        
    end
    
    clearvars -except dim n x y Xq Xr Yq Yr dx dy dX dY L C ...
               Aper Apara Omega V_L D_L pri_up pri_up models ...
               pho_get b_del dtp p_norm

end

Pri_up{1} = real(pri_up);
Pri_up{2} = real(p_norm');

end


