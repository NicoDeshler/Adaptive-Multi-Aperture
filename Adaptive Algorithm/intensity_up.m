function int_new = intensity_up(b, a, g, b_del)

count = 0;

b_old = b;

a_scale = 1/60;
%a_scale = min(sum(a)/5000,1);


epi = 1e-4;


while 0 == 0
    
    count = count + 1;
    
    bg = g.*b;

    b_new = ( a - 1 + bg/sum(bg) )/( 1+ sum(a) - size(b,1));
    
    %{
    
    bg = g'*b;
    
    b_new = ( a - 1 + ( sum( diag(n_mode_count)*g'*diag(b)./repmat(bg, [1, size(b,1)]), 1) )' )...
            /( sum(n_mode_count) + sum(a) - size(b,1));
    
    %}
    
    if sum( abs(b-b_new) <= b_del )
        
        %if all(b_new > 0)
        
            %int_new = [b_new, a + b_new*a_scale];
    
            int_new = [b_new, b_new*(sum(a) + a_scale - size(a,1))+1];
        
            break;
            
        %else
            
            
            
        %end
        
    else
        
        b = b_new;
        
    end
    
end

end
