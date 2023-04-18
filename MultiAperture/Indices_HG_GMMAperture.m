function [pj,qj,uj] = Indices_HG_GMMAperture(n_max,n_apertures)
    % computes the local aperture HG mode indices for a system of
    % n apertures.
    
    pj = [];
    qj = [];
    uj = [];
    
    % aperture, radial, and azimuthal index list     
    for v = 1:n_apertures
        for p = 1:n_max
            for q = 1:p
                pj = [pj, p-q];
                qj = [qj, q-1];
                %pj = [pj, p];
                %qj = [qj, q];
                uj = [uj, v];
            end
        end
    end    
end




