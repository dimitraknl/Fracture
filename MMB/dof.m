function id = dof(type,j,nLayers)
    switch type
        case 'U'
            id = j;
        case 'W'
            id = nLayers+1+j;
        case 'B1'
            id = 2*(nLayers+1)+2*(j-1)+1;
        case 'B2'
            id = 2*(nLayers+1)+2*(j-1)+2;
        case 'D1'
            id = 2*(nLayers+1)+2*nLayers+2*(j-1)+1;
        case 'D2'
            id = 2*(nLayers+1)+2*nLayers+2*(j-1)+2;
        otherwise
            error('Unknown DOF type.');
    end
end
