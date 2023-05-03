function Gamma_nm = corrFn_GramSchmidt_pos(Xq,Yq,X,Y,GS_basis_pos,A_tot)
    % Computes the correlation function for the Gram schmidt basis
    % 
    % Xq,Yq - correlation function query coordinates
    % X,Y - image plane meshgrids
    % GS_modes  - a matrix stack representing the GS modes over X,Y
    % A_tot - total area of the aperture
    
    Gamma_nm = 2*pi/sqrt(A_tot) * conj(Basis_GramSchmidt_pos(Xq,Yq,X,Y,GS_basis_pos));
   
end
