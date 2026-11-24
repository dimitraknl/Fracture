function id = pointIndex(k,ix,iz,Nx,Nz)
    id = (k-1)*Nx*Nz+(iz-1)*Nx+ix;
end
