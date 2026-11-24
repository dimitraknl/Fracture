function beta = evaluateInitialModeRatioLever(leverLength,BH,BUD,BWD,W,Sblock, ...
    nQH,nSigma,Nx,Nz,nLayers,lastReleased,localWBottom,localWTop, ...
    deltaRef,x,xLoad,L,width,newton,contact,vcct)

    contactNodes = zeros(0,1);
    contactStable = false;
    solution = [];
    post = [];

    while ~contactStable
        model = buildMMBModel(BH,BUD,BWD,W,Sblock,nQH,nSigma, ...
            Nx,Nz,nLayers,lastReleased,localWBottom,localWTop, ...
            deltaRef,x,xLoad,L,contactNodes,width,leverLength);

        [solution,~] = solveMixedSystemNewton(model,solution,newton);
        post = recoverMMBSolution(solution,model,W,nQH,Nx,Nz,nLayers);

        if contact.mode == "active_set"
            releasedNodes = (1:lastReleased).';
            penetration = releasedNodes( ...
                post.wD(releasedNodes) < -contact.penetrationTol);
            newContact = setdiff(penetration,contactNodes);
            if isempty(newContact)
                contactStable = true;
            else
                contactNodes = unique([contactNodes;newContact]);
                solution = [];
            end
        elseif contact.mode == "closed_crack"
            contactNodes = (1:lastReleased).';
            if model.nContact == numel(contactNodes)
                contactStable = true;
            else
                solution = [];
            end
        else
            contactStable = true;
        end
    end

    tipNode = lastReleased+1;
    behindNode = tipNode-1;
    da1 = x(tipNode)-x(behindNode);
    if tipNode < Nx
        da2 = x(tipNode+1)-x(tipNode);
    else
        da2 = da1;
    end

    workI = post.tipOpeningReaction*post.wD(behindNode) + ...
        post.tipRotationReaction*post.slopeJump(behindNode);
    workII = post.tipSlidingReaction*post.uD(behindNode);

    if vcct.signMode == "absolute_work"
        GI = abs(workI)/(2*width*da2)*(da2/da1);
        GII = abs(workII)/(2*width*da2)*(da2/da1);
    else
        GI = max(-(workI)/(2*width*da2)*(da2/da1),0);
        GII = max(-(workII)/(2*width*da2)*(da2/da1),0);
    end

    beta = GII/max(GI+GII,eps);
end
