function [B0H,BxH,B0UD,BxUD,B0WD,BxWD] = ...
    localExtendedStrainMatrices(zeta,hk,k,nLayers,nQnodeH)

    [Psi,dPsidz,psi1,dpsi1dz,psi3,dpsi3dz] = ...
        hermiteExtended(zeta,hk);

    B0H = zeros(3,nQnodeH);
    BxH = zeros(3,nQnodeH);

    iu1 = dof('U',k,nLayers);
    iu2 = dof('U',k+1,nLayers);
    iw1 = dof('W',k,nLayers);
    iw2 = dof('W',k+1,nLayers);
    ib1 = dof('B1',k,nLayers);
    ib2 = dof('B2',k,nLayers);
    id1 = dof('D1',k,nLayers);
    id2 = dof('D2',k,nLayers);

    uDofs = [iu1 iu2 ib1 ib2];
    wDofs = [iw1 iw2 id1 id2];

    BxH(1,uDofs) = Psi;
    B0H(2,wDofs) = dPsidz;
    B0H(3,uDofs) = dPsidz;
    BxH(3,wDofs) = Psi;

    if k == 1
        s = -0.5;
        psiD = psi3;
        dPsiDdz = dpsi3dz;
    elseif k == 2
        s = 0.5;
        psiD = psi1;
        dPsiDdz = dpsi1dz;
    else
        error('This implementation currently assumes two sublaminates.');
    end

    % uD contributions:
    %   epsilon_x = s*PsiD*uD,x
    %   gamma_xz  = s*PsiD,z*uD
    B0UD = [0;0;s*dPsiDdz];
    BxUD = [s*psiD;0;0];

    % wD contributions:
    %   epsilon_z = s*PsiD,z*wD
    %   gamma_xz  = s*PsiD*wD,x
    B0WD = [0;s*dPsiDdz;0];
    BxWD = [0;0;s*psiD];
end
