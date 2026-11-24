function post = recoverModeIISolution( ...
    solution,model,W,nQH,Nx,Nz,nLayers)

    nQf = model.nQf;
    nSigma = model.nSigma;
    nLt = model.nLt;
    nCq = model.nCtip;

    qFree = solution(1:nQf);
    sigma = solution(nQf+(1:nSigma));
    lambdaTraction = solution(nQf+nSigma+(1:nLt)); %#ok<NASGU>
    muQ = solution(nQf+nSigma+nLt+(1:nCq));

    qact = zeros(model.nQact,1);
    qact(model.freeQ) = qFree;
    qact(model.prescribedDofs) = model.prescribedVals;

    qH = qact(1:nQH);
    aUD = qact(model.iUD);
    aWD = qact(model.iWD);

    uD = zeros(Nx,1);
    wD = zeros(Nx,1);
    uD(model.enrichedNodes) = aUD;
    wD(model.enrichedNodes) = aWD;

    nQnodeH = nQH/Nx;
    qNodeH = zeros(Nx,nQnodeH);
    for d = 1:nQnodeH
        qNodeH(:,d) = qH((d-1)*Nx+(1:Nx));
    end

    sigmaX = cell(nLayers,1);
    sigmaZ = cell(nLayers,1);
    tauXZ = cell(nLayers,1);

    for k = 1:nLayers
        sigmaX{k} = zeros(Nx,Nz);
        sigmaZ{k} = zeros(Nx,Nz);
        tauXZ{k} = zeros(Nx,Nz);

        for ix = 1:Nx
            for iz = 1:Nz
                sigmaX{k}(ix,iz) = ...
                    sigma(stressIndex(k,ix,iz,1,Nx,Nz));
                sigmaZ{k}(ix,iz) = ...
                    sigma(stressIndex(k,ix,iz,2,Nx,Nz));
                tauXZ{k}(ix,iz) = ...
                    sigma(stressIndex(k,ix,iz,3,Nx,Nz));
            end
        end
    end

    qReaction = model.Bact.'*W*sigma + model.Ctip.'*muQ;

    leftSupportReaction  = qReaction(model.leftSupportDof);
    rightSupportReaction = qReaction(model.rightSupportDof);
    loadReaction         = qReaction(model.loadDof);
    reactionForce        = abs(loadReaction);
    supportImbalance = abs(leftSupportReaction+rightSupportReaction+ ...
        loadReaction)/max([abs(leftSupportReaction), ...
        abs(rightSupportReaction),abs(loadReaction),eps]);

    tipSlidingReaction  = muQ(1);
    tipOpeningReaction  = muQ(2);
    tipRotationReaction = muQ(3);

    if model.nContact > 0
        contactReactions = muQ(4:end);
    else
        contactReactions = zeros(0,1);
    end

    uInterfaceHealthy = zeros(Nx,1);
    wInterfaceHealthy = zeros(Nx,1);
    uInterfaceUpper = zeros(Nx,1);
    uInterfaceLower = zeros(Nx,1);
    wInterfaceUpper = zeros(Nx,1);
    wInterfaceLower = zeros(Nx,1);
    wOuterUpper = zeros(Nx,1);
    wOuterLower = zeros(Nx,1);
    slopeJump = zeros(Nx,1);

    for ix = 1:Nx
        [uInt,wInt,b,~] = unpackCSLTNode(qNodeH(ix,:).',nLayers);

        uInterfaceHealthy(ix) = uInt(2);
        wInterfaceHealthy(ix) = wInt(2);
        uInterfaceUpper(ix) = uInt(2)+0.5*uD(ix);
        uInterfaceLower(ix) = uInt(2)-0.5*uD(ix);
        wInterfaceUpper(ix) = wInt(2)+0.5*wD(ix);
        wInterfaceLower(ix) = wInt(2)-0.5*wD(ix);
        wOuterLower(ix) = wInt(1);
        wOuterUpper(ix) = wInt(end);
        slopeJump(ix) = b(1,2)-b(2,1);
    end

    midpointDisplacement = abs(wOuterUpper(model.loadNode));

    prescribedRecovered = qact(model.prescribedDofs);
    boundaryConditionError = max(abs( ...
        prescribedRecovered-model.prescribedVals));

    post.qact = qact;
    post.qNodeH = qNodeH;
    post.uD = uD;
    post.wD = wD;
    post.slopeJump = slopeJump;
    post.sigmaX = sigmaX;
    post.sigmaZ = sigmaZ;
    post.tauXZ = tauXZ;
    post.reactionForce = reactionForce;
    post.leftSupportReaction = leftSupportReaction;
    post.rightSupportReaction = rightSupportReaction;
    post.loadReaction = loadReaction;
    post.supportImbalance = supportImbalance;
    post.tipSlidingReaction = tipSlidingReaction;
    post.tipOpeningReaction = tipOpeningReaction;
    post.tipRotationReaction = tipRotationReaction;
    post.contactNodes = model.contactNodes;
    post.contactReactions = contactReactions;
    post.uInterfaceHealthy = uInterfaceHealthy;
    post.wInterfaceHealthy = wInterfaceHealthy;
    post.uInterfaceUpper = uInterfaceUpper;
    post.uInterfaceLower = uInterfaceLower;
    post.wInterfaceUpper = wInterfaceUpper;
    post.wInterfaceLower = wInterfaceLower;
    post.wOuterUpper = wOuterUpper;
    post.wOuterLower = wOuterLower;
    post.midpointDisplacement = midpointDisplacement;
    post.boundaryConditionError = boundaryConditionError;
    post.interfaceSigmaZLower = sigmaZ{1}(:,end);
    post.interfaceSigmaZUpper = sigmaZ{2}(:,1);
    post.interfaceTauLower = tauXZ{1}(:,end);
    post.interfaceTauUpper = tauXZ{2}(:,1);
end
