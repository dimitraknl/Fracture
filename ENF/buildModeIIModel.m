function model = buildModeIIModel( ...
    BH,BUD,BWD,W,Sblock,nQH,nSigma,Nx,Nz,nLayers, ...
    lastReleased,localWBottom,localWTop,deltaControl, ...
    x,xLoad,contactNodes,width)

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

    Cq = sparse(3+nContact,nQact);
    Cq(1,iUD(end)) = 1;
    Cq(2,iWD(end)) = 1;

    idBUpperBottom = qHIndex(tipNode,dof('B1',2,nLayers),Nx);
    idBLowerTop    = qHIndex(tipNode,dof('B2',1,nLayers),Nx);
    Cq(3,idBUpperBottom) =  1;
    Cq(3,idBLowerTop)    = -1;

    for c = 1:nContact
        localIndex = find(enrichedNodes == contactNodes(c),1);
        Cq(3+c,iWD(localIndex)) = 1;
    end

    [~,loadNode] = min(abs(x-xLoad));

    % ENF three-point-bending boundary conditions:
    %   lower surface, x=0:   w=0 (left roller)
    %   lower surface, x=L:   w=0 (right roller)
    %   upper surface, x=L/2: w=-deltaControl (loading roller)
    %   lower surface, x=0:   u=0 (remove axial rigid-body mode)
    prescribedDofs = [ ...
        qHIndex(1,localWBottom,Nx); ...
        qHIndex(Nx,localWBottom,Nx); ...
        qHIndex(loadNode,localWTop,Nx); ...
        qHIndex(1,dof('U',1,nLayers),Nx)];

    prescribedVals = [0;0;-deltaControl;0];

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
        -Cqp*prescribedVals];

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
    model.leftSupportDof = prescribedDofs(1);
    model.rightSupportDof = prescribedDofs(2);
    model.loadDof = prescribedDofs(3);
    model.axialAnchorDof = prescribedDofs(4);
end
