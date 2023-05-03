function mom_Bay = mom_Bay(t,s,mom)
 % Computes the mom^th non-central moment of the normal distribution.
 % The following integral is computed with the Kummer complex function:
 %
 % integral  Normal(x;t,s) * exp(-x^2 / (2 * delta)^2 ) * (x / (2 * delta)) ^ (mom)  dx
 %
 % Reference: https://arxiv.org/pdf/1209.4340.pdf
 %
 % where mom is related to the index of the HG mode and the moment of the
 % Gamma_ik operator.
 %
 % x is a source position coordinate (either xi or yi)
 % with Gaussian posterior parameters mu_i = t and sigma_i = s.
 % delta is the Gaussian PSF width.
 
 
 
 s_bar = % stdev of the normal distribution of the prior
 s = % stdev of the normal distribution of the constituent gaussian psfs
  
 
 
 

    T = t/(2*s^2+1);
    S = sqrt( s^2/(2*s^2+1) );
        %{
        if mod(mom,2) == 0
            
            mom_Bay = exp(-t^2)/sqrt(2*s^2+1) ...
                    *hypergeom(-mom/2, 1/2, -T^2/2/S^2 ) ...
                    *S^mom*2^(mom/2)*gamma((mom+1)/2)/sqrt(pi);
    
        else
        
            mom_Bay = exp(-t^2)/sqrt(2*s^2+1) ...
                    *hypergeom((1-mom)/2, 3/2, -T^2/2/S^2 ) ...
                    *T*S^(mom-1)*2^((mom+1)/2)*gamma(mom/2+1)/sqrt(pi);

        end
        %}
    
        if mod(mom,2) == 0
            
            mom_Bay = exp(-t^2)/sqrt(2*s^2+1) ...
                    *KummerComplex(-mom/2, 1/2, -T^2/2/S^2 ) ...
                    *S^mom*2^(mom/2)*gamma((mom+1)/2)/sqrt(pi);
    
        else
        
            mom_Bay = exp(-t^2)/sqrt(2*s^2+1) ...
                    *KummerComplex((1-mom)/2, 3/2, -T^2/2/S^2 ) ...
                    *T*S^(mom-1)*2^((mom+1)/2)*gamma(mom/2+1)/sqrt(pi);

        end

end