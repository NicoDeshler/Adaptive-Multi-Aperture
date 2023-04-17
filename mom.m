function mom = mom(n, x, xj, sig)

[Xq, Xr, Xj] = meshgrid(x, x, xj); % r var in dim 2; q var in dim 1
Sig = repmat( permute(sig, [3,2,1]), [size(x,1), size(x,1), 1] ); % here sig is the variance not sd, so no need square

mu = exp( - ( (Xj - Xr).^2 + (Xj - Xq).^2 + (Xr - Xq).^2.*Sig )./(2*Sig + 1)/2 );

switch n
    
    case 0
        
        lam = 1./sqrt(2*Sig + 1);
 
        
    case 1
        
        lam = ( (Xr + Xq).*Sig + Xj )./sqrt(2*Sig + 1).^3;

        
    case 2
        
        lam = ( ( (Xr + Xq).*Sig + Xj ).^2 + Sig.*(2*Sig + 1) )...
             ./sqrt(2*Sig + 1).^5;

        
    case 3
        
        lam = ( ( (Xr + Xq).*Sig + Xj ).^3 ...
             +3*( (Xr + Xq).*Sig + Xj ).*Sig.*(2*Sig + 1) )...
             ./sqrt(2*Sig + 1).^7;

        
    case 4
        
        lam = ( ( (Xr + Xq).*Sig + Xj ).^4 ...
             +6*( (Xr + Xq).*Sig + Xj ).^2.*Sig.*(2*Sig + 1) ...
             +3*Sig.^2.*(2*Sig + 1).^2 )...
             ./sqrt(2*Sig + 1).^9;
end

mom = lam.*mu;

end



