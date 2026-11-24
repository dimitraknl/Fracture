function model = buildFullFieldModel( ...
    BH,BUD,BWD,W,Sblock,nQH,nSigma,Nx,Nz,nLayers, ...
    lastReleased,localWBottom,localWTop,deltaControl, ...
    armConfiguration,width)

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

    % Traction continuity over the bonded interface, including the tip.
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

    % Trial-tip constraints:
    %   uD_tip = 0
    %   wD_tip = 0
    %   b_upper,bottom - b_lower,top = 0
    Ctip = sparse(3,nQact);
    Ctip(1,iUD(end)) = 1;
    Ctip(2,iWD(end)) = 1;

    idBUpperBottom = qHIndex(tipNode,dof('B1',2,nLayers),Nx);
    idBLowerTop    = qHIndex(tipNode,dof('B2',1,nLayers),Nx);
    Ctip(3,idBUpperBottom) =  1;
    Ctip(3,idBLowerTop)    = -1;

    % Prescribed outer-surface DCB displacement control.
    % deltaControl is always the displacement prescribed to the moving arm.
    switch armConfiguration
        case "single_arm"
            % One-arm test: move the upper arm and hold the lower arm fixed.
            wTopValue = deltaControl;
            wBottomValue = 0.0;
        case "two_arm"
            % Symmetric two-arm test: equal and opposite arm displacements.
            wTopValue =  deltaControl;
            wBottomValue = -deltaControl;
        otherwise
            error('Unknown armConfiguration. Use "single_arm" or "two_arm".');
    end

    prescribedDofs = [ ...
        qHIndex(1,localWTop,Nx); ...
        qHIndex(1,localWBottom,Nx); ...
        qHIndex(Nx,dof('U',2,nLayers),Nx)];

    prescribedVals = [wTopValue; wBottomValue; 0];

    allQ = (1:nQact).';
    freeQ = setdiff(allQ,prescribedDofs);

    Bf = Bact(:,freeQ);
    Bp = Bact(:,prescribedDofs);
    Ctf = Ctip(:,freeQ);
    Ctp = Ctip(:,prescribedDofs);

    nQf = numel(freeQ);
    nCtip = size(Ctip,1);

    Kmix = [ ...
        sparse(nQf,nQf), Bf.'*W,             sparse(nQf,nLt), Ctf.';
        W*Bf,              -W*Sblock,        G.',              sparse(nSigma,nCtip);
        sparse(nLt,nQf),   G,                sparse(nLt,nLt),  sparse(nLt,nCtip);
        Ctf,                sparse(nCtip,nSigma), sparse(nCtip,nLt), sparse(nCtip,nCtip)];

    rhs = [ ...
        zeros(nQf,1); ...
        -W*Bp*prescribedVals; ...
        zeros(nLt,1); ...
        -Ctp*prescribedVals];

    model.Kmix = Kmix;
    model.rhs = rhs;

    % Residual block indices for ODIL diagnostics.
    model.rows.q = (1:nQf).';
    model.rows.constitutive = nQf+(1:nSigma).';
    model.rows.interface = nQf+nSigma+(1:nLt).';
    model.rows.tip = nQf+nSigma+nLt+(1:nCtip).';

    model.Bact = Bact;
    model.freeQ = freeQ;
    model.prescribedDofs = prescribedDofs;
    model.prescribedVals = prescribedVals;
    model.armConfiguration = armConfiguration;
    model.nQf = nQf;
    model.nSigma = nSigma;
    model.nLt = nLt;
    model.nCtip = nCtip;
    model.nQact = nQact;
    model.nA = nA;
    model.iUD = iUD;
    model.iWD = iWD;
    model.enrichedNodes = enrichedNodes;
    model.tipNode = tipNode;
    model.Ctip = Ctip;
    model.width = width;
end
