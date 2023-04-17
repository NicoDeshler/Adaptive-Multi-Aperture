function [pj,qj,uj] = Indices_HG_GMMAperture(n_max,n_apertures)
    % computes the local aperture HG mode indices for a system of
    % n apertures.
    
    pj = [];
    qj = [];
    uj = [];
    
    % aperture, radial, and azimuthal index list     
    for v = 1:n_apertures
        for p = 0:n_max
            for q = 0:p
                pj = [pj, p];
                qj = [qj, q];
                uj = [uj, v];
            end
        end
    end    
end
