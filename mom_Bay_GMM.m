function mom_Bay = mom_Bay_GMM(t,s,mom)
 % Computes the mom^th non-central moment of the normal distribution.
 % The following integral is computed with the Parabolic Cylider function:
 %
 % integral  Normal(x;t_bar,s_bar) * exp(-x^2 / (2 * s)^2 ) * (x / (2 * s)) ^ (mom)  dx
 %
 % Reference: https://arxiv.org/pdf/1209.4340.pdf
 %
 % where mom is related to the index of the HG mode and the moment of the
 % Gamma_ik operator.
 %
 % x is a source position coordinate (either xi or yi) 
 % with Gaussian posterior parameters t_bar = mu_i and s_bar = sigma_i. 
 % s is the Gaussian PSF standard deviation.
 
 
 b = % baseline component (x or y)
 s = % stdev of the normal distribution of the constituent gaussian psfs
 s_bar = % stdev of the normal distribution of the prior
 t_bar = % mean of the normal distribution of the prior
 
 T = 2*s^2*(t_bar - 1i*b*s_bar^2) /(s_bar^2 + 2*s^2);
 S = sqrt(2)*s*s_bar / sqrt(s_bar^2 + 2*s^2);
 
 
 mom_Bay = (1i * S)^(mom) * exp(- (T / (2*S))^2 ) * pu(-1i * (T/S));

end