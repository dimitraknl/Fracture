function C = transformedPlaneStrainApproximation(theta,mat)

    if abs(theta) < 1e-12
        Ex = mat.E11;
        Ez = mat.E33;
        Gxz = mat.G13;
        nuxz = mat.nu13;
    elseif abs(abs(theta)-90) < 1e-12
        Ex = mat.E22;
        Ez = mat.E33;
        Gxz = mat.G23;
        nuxz = mat.nu23;
    else
        m = cosd(theta);
        n = sind(theta);
        Ex = 1/(m^4/mat.E11+n^4/mat.E22+ ...
            m^2*n^2*(1/mat.G12-2*mat.nu12/mat.E11));
        Ez = mat.E33;
        Gxz = mat.G13*m^2+mat.G23*n^2;
        nuxz = mat.nu13*m^2+mat.nu23*n^2;
    end

    nuzx = nuxz*Ez/Ex;
    den = 1-nuxz*nuzx;

    C = [Ex/den, nuxz*Ez/den, 0; ...
         nuzx*Ex/den, Ez/den, 0; ...
         0, 0, Gxz];
end
