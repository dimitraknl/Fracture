function [uInt,wInt,b,d] = unpackCSLTNode(qLocal,nLayers)

    nInterface = nLayers+1;
    idx = 0;

    uInt = qLocal(idx+(1:nInterface));
    idx = idx+nInterface;

    wInt = qLocal(idx+(1:nInterface));
    idx = idx+nInterface;

    b = zeros(2,nLayers);
    for k = 1:nLayers
        b(1,k) = qLocal(idx+1);
        b(2,k) = qLocal(idx+2);
        idx = idx+2;
    end

    d = zeros(2,nLayers);
    for k = 1:nLayers
        d(1,k) = qLocal(idx+1);
        d(2,k) = qLocal(idx+2);
        idx = idx+2;
    end
end
