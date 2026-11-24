function id = stressIndex(k,ix,iz,comp,Nx,Nz)
    p = pointIndex(k,ix,iz,Nx,Nz);
    id = 3*(p-1)+comp;
end
