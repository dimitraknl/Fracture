function [Psi,dPsidz,psi1,dpsi1dz,psi3,dpsi3dz] = ...
    hermiteExtended(zeta,hk)

    psi1 = 0.25*(1-zeta)^2*(2+zeta);
    psi2 = 0.25*(1-zeta)^2*(1+zeta);
    psi3 = 0.25*(1+zeta)^2*(2-zeta);
    psi4 = 0.25*(1+zeta)^2*(zeta-1);

    dpsi1 = 0.25*(-2*(1-zeta)*(2+zeta)+(1-zeta)^2);
    dpsi2 = 0.25*(-2*(1-zeta)*(1+zeta)+(1-zeta)^2);
    dpsi3 = 0.25*( 2*(1+zeta)*(2-zeta)-(1+zeta)^2);
    dpsi4 = 0.25*( 2*(1+zeta)*(zeta-1)+(1+zeta)^2);

    Psi = [psi1 psi3 hk/2*psi2 hk/2*psi4];
    dPsidz = (2/hk)*[dpsi1 dpsi3 hk/2*dpsi2 hk/2*dpsi4];

    dpsi1dz = (2/hk)*dpsi1;
    dpsi3dz = (2/hk)*dpsi3;
end
