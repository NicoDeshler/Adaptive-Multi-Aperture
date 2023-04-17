function [like_temp, a] = reSacle_like(log_like, n_sam, scale)

            disp('Now rescale')
        
            a = fix(log(realmax))*scale - ( max(log_like) - log(n_sam) );
        
            if sum(exp(log_like + a - log(n_sam))) == inf
            
                like_temp = reSacle_like(log_like, n_sam, scale*0.99);
           
            else
            
                log_like = log_like + a;
            
                like_temp = exp(log_like - log(n_sam));
            
            end



end