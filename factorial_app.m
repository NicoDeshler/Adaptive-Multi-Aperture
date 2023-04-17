function out = factorial_app(n, n_max)

    %n_max = 120;
    
    n = [sum(n), n];
    
    n_app = n( n > n_max );
    
    n_exact = n( n <= n_max );
    
    if isempty(n_app)
        
        %disp('All good')
        
        out = factorial(n(1))/prod(factorial(n(2:end)));
        
    elseif all( n > n_max )
        
        %disp('All bad')
        
        out = exp( log(2*pi*n(1)) + n(1)*log(n(1)/exp(1)) ...
                  -log(sum(2*pi*prod(n(2:end)) + n(2:end).*log(n(2:end)/exp(1)),2) ) );
             
    else
        
        %disp('Some bad')
        
        if numel(n_app) == 1
            
            %disp("one")
            
            out = exp( log(2*pi*n_app(1)) + n_app(1)*log(n_app(1)/exp(1)) ) ...
                 /prod(factorial(n_exact));
            
        else
        
            out = exp( log(2*pi*n_app(1)) + n_app(1)*log(n_app(1)/exp(1)) ...
                     -(log(2*pi*prod(n_app(2:end))) + n_app(2:end).*log(n_app(2:end)/exp(1)) ) ) ...
                  /prod(factorial(n_exact));
             
        end
        
    end


end