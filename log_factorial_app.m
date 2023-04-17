function out = log_factorial_app(n, n_max)
    
    n_app = n > n_max;
    
    if all(all(~n_app))
        
        %disp('All good')
        
        out = log(factorial(n));
        
    elseif all(all(n_app))
        
        %disp('All bad')
        
        out = (n.*log(n + (n == 0)) - n + log(n + (n == 0))/2 + log(2*pi)/2).*(n ~= 0) ;
             
    else

        out = ( log(factorial( n.*(~n_app) ) ) ...
            + (n.*log(n + (n == 0)) - n + log(n + (n == 0))/2 + log(2*pi)/2).*n_app ).*(n ~= 0);
        
    end


end