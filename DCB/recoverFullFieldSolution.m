function post = recoverFullFieldSolution( ...
    solution,model,W,nQH,Nx,Nz,nLayers)

    nQf = model.nQf;
    nSigma = model.nSigma;
    nLt = model.nLt;
    nCtip = model.nCtip;

    qFree = solution(1:nQf);
    sigma = solution(nQf+(1:nSigma));
    lambdaTraction = solution(nQf+nSigma+(1:nLt)); %#ok<NASGU>
    muTip = solution(nQf+nSigma+nLt+(1:nCtip));

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

    % Generalized reactions, including trial-tip constraint reactions.
    qReaction = model.Bact.'*W*sigma + model.Ctip.'*muTip;

    upperDof = model.prescribedDofs(1);
    lowerDof = model.prescribedDofs(2);
    upperReaction = qReaction(upperDof);
    lowerReaction = qReaction(lowerDof);
    switch model.armConfiguration
        case "single_arm"
            % Report the reaction conjugate to the prescribed upper-arm motion.
            reactionForce = abs(upperReaction);
        case "two_arm"
            % Equal and opposite reactions; report the force carried by one arm.
            reactionForce = 0.5*(abs(upperReaction)+abs(lowerReaction));
        otherwise
            error('Unknown armConfiguration in post-processing.');
    end
    reactionImbalance = abs(upperReaction+lowerReaction)/ ...
        max([abs(upperReaction),abs(lowerReaction),eps]);

    tipSlidingReaction  = muTip(1);
    tipOpeningReaction  = muTip(2);
    tipRotationReaction = muTip(3);

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

    actualCMOD = wOuterUpper(1)-wOuterLower(1);
    meanMouthMotion = 0.5*(wOuterUpper(1)+wOuterLower(1));

    post.qact = qact;
    post.qNodeH = qNodeH;
    post.uD = uD;
    post.wD = wD;
    post.slopeJump = slopeJump;
    post.sigmaX = sigmaX;
    post.sigmaZ = sigmaZ;
    post.tauXZ = tauXZ;
    post.reactionForce = reactionForce;
    post.upperReaction = upperReaction;
    post.lowerReaction = lowerReaction;
    post.reactionImbalance = reactionImbalance;
    post.KmixSize = size(model.Kmix);
    post.tipSlidingReaction = tipSlidingReaction;
    post.tipOpeningReaction = tipOpeningReaction;
    post.tipRotationReaction = tipRotationReaction;
    post.uInterfaceHealthy = uInterfaceHealthy;
    post.wInterfaceHealthy = wInterfaceHealthy;
    post.uInterfaceUpper = uInterfaceUpper;
    post.uInterfaceLower = uInterfaceLower;
    post.wInterfaceUpper = wInterfaceUpper;
    post.wInterfaceLower = wInterfaceLower;
    post.wOuterUpper = wOuterUpper;
    post.wOuterLower = wOuterLower;
    post.actualCMOD = actualCMOD;
    post.meanMouthMotion = meanMouthMotion;
    post.interfaceSigmaZLower = sigmaZ{1}(:,end);
    post.interfaceSigmaZUpper = sigmaZ{2}(:,1);
end
