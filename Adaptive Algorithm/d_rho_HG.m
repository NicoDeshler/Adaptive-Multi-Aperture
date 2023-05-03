function [db, dx, dy] = d_rho_HG(n_modes, scene)

    HG_mode = HG_projection(n_modes, scene);
    
    rho = rho_HG(n_modes, scene);

    xy_ind = [[HG_mode.ind_x]',[HG_mode.ind_y]'];
    
    for i = 1:size(xy_ind,1)
        for j = 1:i
           
           d_HG(i,j).ind = [xy_ind(i,:);xy_ind(j,:)];

           db(i,j,:) = permute( ...
                                  prod( scene(:,2:3).^repmat(sum(d_HG(i,j).ind,1), [size(scene,1),1]) ...
                                      .*exp( -scene(:,2:3).^2 ) , 2) - rho(i,j), ...
                               [3,2,1] )/sqrt(prod(factorial(d_HG(i,j).ind(:)),1));
                           
           dx(i,j,:) = permute( ...
                                    prod( scene(:,2:3).^repmat(sum(d_HG(i,j).ind,1), [size(scene,1),1]) , 2)...
                                  .*( sum(d_HG(i,j).ind(:,1), 1)./( scene(:,2) + ( scene(:,2) == 0 ) ).*( scene(:,2) ~= 0 ) ...
                                    - 2*scene(:,2) ).*scene(:,1).*exp( -sum( scene(:,2:3).^2, 2) ), ...
                               [3,2,1] )/sqrt(prod(factorial(d_HG(i,j).ind(:)),1));
                           
          dy(i,j,:) = permute( ...
                                    prod( scene(:,2:3).^repmat(sum(d_HG(i,j).ind,1), [size(scene,1),1]) , 2)...
                                   .*( sum(d_HG(i,j).ind(:,2), 1)./( scene(:,3) + ( scene(:,3) == 0 ) ).*( scene(:,3) ~= 0 ) ...
                                    - 2*scene(:,3) ).*scene(:,1).*exp( -sum( scene(:,2:3).^2, 2 ) ), ...
                               [3,2,1] )/sqrt(prod(factorial(d_HG(i,j).ind(:)),1));
           
            db(j,i,:) = db(i,j,:);
            
            dx(j,i,:) = dx(i,j,:);
            
            dy(j,i,:) = dy(i,j,:);
                                  

        end
    end

end