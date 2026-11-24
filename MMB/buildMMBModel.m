function model = buildMMBModel( ...
    BH,BUD,BWD,W,Sblock,nQH,nSigma,Nx,Nz,nLayers, ...
    lastReleased,localWBottom,localWTop,deltaControl, ...
    x,xLoad,L,contactNodes,width,leverLength)

    %#ok<INUSD>

    tipNode = lastReleased+1;
    enrichedNodes = (1:tipNode).';
    nA = numel(enrichedNodes);

    BUDact = BUD(:,enrichedNodes);
    BWDact = BWD(:,enrichedNodes);
    Bact = [BH,BUDact,BWDact];

    nQact = nQH+2*nA;
    iUD = nQH+(1:nA);
    iWD = nQH+nA+(1:nA);

    bondedNodes = tipNode:Nx;
    nBonded = numel(bondedNodes);
    nLt = 2*nBonded;
    G = spalloc(nLt,nSigma,2*nLt);

    row = 0;
    for ix = bondedNodes
        row = row+1;
        G(row,stressIndex(1,ix,Nz,2,Nx,Nz)) =  1;
        G(row,stressIndex(2,ix,1,2,Nx,Nz)) = -1;

        row = row+1;
        G(row,stressIndex(1,ix,Nz,3,Nx,Nz)) =  1;
        G(row,stressIndex(2,ix,1,3,Nx,Nz)) = -1;
    end

    contactNodes = unique(contactNodes(:));
    contactNodes = contactNodes(contactNodes >= 1 & ...
        contactNodes <= lastReleased);
    nContact = numel(contactNodes);

    % Constraint rows:
    %   1: crack-tip sliding uD = 0
    %   2: crack-tip opening wD = 0
    %   3: crack-tip relative rotation = 0
    %   4: one rigid-lever actuator displacement constraint
    %   5...: unilateral contact constraints
    Cq = sparse(4+nContact,nQact);
    gq = zeros(4+nContact,1);

    Cq(1,iUD(end)) = 1;
    Cq(2,iWD(end)) = 1;

    idBUpperBottom = qHIndex(tipNode,dof('B1',2,nLayers),Nx);
    idBLowerTop    = qHIndex(tipNode,dof('B2',1,nLayers),Nx);
    Cq(3,idBUpperBottom) =  1;
    Cq(3,idBLowerTop)    = -1;

    [~,loadNode] = min(abs(x-xLoad));
    endNode = 1;

    alphaEnd = -2*leverLength/L;
    alphaMid = 1 + 2*leverLength/L;

    endLeverDof = qHIndex(endNode,localWTop,Nx);
    midLeverDof = qHIndex(loadNode,localWTop,Nx);

    Cq(4,endLeverDof) = alphaEnd;
    Cq(4,midLeverDof) = alphaMid;
    gq(4) = -deltaControl;

    for c = 1:nContact
        localIndex = find(enrichedNodes == contactNodes(c),1);
        Cq(4+c,iWD(localIndex)) = 1;
    end

    % Only physical supports and one axial anchor are prescribed directly.
    % The end and midpoint vertical motions remain unknown and are coupled
    % only by the single rigid-lever equation above.
    prescribedDofs = [ ...
        qHIndex(1,localWBottom,Nx); ...
        qHIndex(Nx,localWBottom,Nx); ...
        qHIndex(1,dof('U',1,nLayers),Nx)];

    prescribedVals = [0;0;0];

    allQ = (1:nQact).';
    freeQ = setdiff(allQ,prescribedDofs);

    Bf = Bact(:,freeQ);
    Bp = Bact(:,prescribedDofs);
    Cqf = Cq(:,freeQ);
    Cqp = Cq(:,prescribedDofs);

    nQf = numel(freeQ);
    nCq = size(Cq,1);

    Kmix = [ ...
        sparse(nQf,nQf), Bf.'*W,               sparse(nQf,nLt), Cqf.';
        W*Bf,              -W*Sblock,          G.',              sparse(nSigma,nCq);
        sparse(nLt,nQf),   G,                  sparse(nLt,nLt),  sparse(nLt,nCq);
        Cqf,                sparse(nCq,nSigma), sparse(nCq,nLt), sparse(nCq,nCq)];

    rhs = [ ...
        zeros(nQf,1); ...
        -W*Bp*prescribedVals; ...
        zeros(nLt,1); ...
        gq-Cqp*prescribedVals];

    model.Kmix = Kmix;
    model.rhs = rhs;
    model.rows.q = (1:nQf).';
    model.rows.constitutive = nQf+(1:nSigma).';
    model.rows.interface = nQf+nSigma+(1:nLt).';
    model.rows.tip = nQf+nSigma+nLt+(1:nCq).';

    model.Bact = Bact;
    model.freeQ = freeQ;
    model.prescribedDofs = prescribedDofs;
    model.prescribedVals = prescribedVals;
    model.nQf = nQf;
    model.nSigma = nSigma;
    model.nLt = nLt;
    model.nCtip = nCq;
    model.nQact = nQact;
    model.nA = nA;
    model.iUD = iUD;
    model.iWD = iWD;
    model.enrichedNodes = enrichedNodes;
    model.tipNode = tipNode;
    model.Ctip = Cq;
    model.width = width;
    model.contactNodes = contactNodes;
    model.nContact = nContact;
    model.loadNode = loadNode;
    model.endNode = endNode;
    model.leftSupportDof = prescribedDofs(1);
    model.rightSupportDof = prescribedDofs(2);
    model.axialAnchorDof = prescribedDofs(3);
    model.endLeverDof = endLeverDof;
    model.midLeverDof = midLeverDof;
    model.leverConstraintRow = 4;
    model.alphaEnd = alphaEnd;
    model.alphaMid = alphaMid;
    model.leverLength = leverLength;
    model.constraintRhs = gq;
end
