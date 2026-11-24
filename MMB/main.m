%% ODIL-XCSLT NODEWISE VCCT: MIXED-MODE BENDING (MMB)
% Mixed stress-displacement ODIL formulation for the MMB benchmark of
% Siorikis et al. (2021), using the validated single-domain
% Chebyshev-Lobatto discretization.
%
% Main features
%   - Cubic-Hermite layerwise kinematics through the thickness
%   - Single-domain Chebyshev-Lobatto collocation in x
%   - Independent sigma_x, sigma_z and tau_xz stress fields
%   - Full sliding uD(x) and opening wD(x) discontinuities
%   - Physical rigid MMB lever represented by one multipoint constraint
%   - One prescribed actuator displacement and one work-conjugate load
%   - Optional calibration of lever length c to target GII/(GI+GII)
%   - Active-set crack-face contact
%   - Nodewise VCCT on the nonuniform Lobatto grid
%   - Benzeggagh-Kenane mixed-mode propagation criterion
%   - Newton solution of the mixed stationarity residual
%   - ODIL loss 0.5*R''*R used only as convergence diagnostic
%
% Newton update:
%       Kmix*dy = -R,    R = Kmix*y-rhs
%
% No normal equations and no conventional FE stiffness assembly are used.

clear; clc; close all;

% Results folder
resultsDir = fullfile(pwd,'results3');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

%% ========================================================================
% 1. USER CONTROLS
% ========================================================================

analysisMode = "propagation";    % "fixed_crack" or "propagation"

grid.mode = "single_domain";      % validated branch
grid.single.Nx = 41;

% Multidomain data are retained only for future work. The present
% MMB validation deliberately uses the stable single-domain branch.
grid.multi.breaks = [];
grid.multi.nodes  = [];

% Positive actuator displacement magnitude at the external lever load point.
% The rigid-lever kinematic constraint is
%   delta_P = (1+2*c/L)*w_mid - (2*c/L)*w_end.
% Downward actuator motion is imposed as delta_P = -delta.
loadControl.fixedDisplacement = 0.10e-3;
loadControl.maxDisplacement   = 6.0e-3;
loadControl.nSteps            = 25;
loadControl.explicitHistory   = [];

% MMB case selection and rigid-lever calibration
% Available benchmark cases: "20_percent", "50_percent", "80_percent".
mixedMode.case = "80_percent";
mixedMode.calibrateLeverLength = true;
mixedMode.initialLeverLength = 45e-3;     % c [m]
mixedMode.leverLengthLower = 0.0;         % [m]
mixedMode.leverLengthUpper = 250e-3;      % [m]
mixedMode.ratioTolerance = 1.0e-4;
mixedMode.maxCalibrationIter = 35;
mixedMode.referenceDisplacement = 0.05e-3;
mixedMode.verboseCalibration = true;
mixedMode.scanPoints = 60;

% Newton controls
newton.stopMode       = "tolerance";
newton.maxIter        = 8;
newton.minIter        = 1;
newton.fixedIter      = 2;
newton.damping        = 1.0;
newton.tolRelResidual = 1.0e-4;
newton.tolRelStep     = 1.0e-12;
newton.verbose        = false;
newton.warmStart      = true;

% Matrix diagnostics
newton.checkMatrix                 = true;
newton.failOnStructuralSingularity = true;
newton.conditionWarning            = 1.0e14;
newton.symmetryTolerance           = 1.0e-10;
newton.printBlockResiduals         = false;

% Contact controls
%   "active_set"   : add wD=0 constraints only where penetration occurs
%   "closed_crack" : keep all released crack nodes closed
%   "none"         : diagnostic only
contact.mode           = "active_set";
contact.penetrationTol = 1.0e-10;
contact.maxUpdates     = 12;
contact.verbose        = false;

% VCCT controls
vcct.signMode = "absolute_work";   % "absolute_work" or "signed"
vcct.printDiagnostics = false;
vcct.modeRatioTolerance = 0.05;

% Plots
plotOptions.showStressProfiles = true;
plotOptions.showContactStatus  = true;

% Timing
timing.enabled       = true;
timing.printPerStep  = false;
timing.plotPerStep   = true;
timing.totalTimer    = tic;
timing.setupTimer    = tic;

%% ========================================================================
% 2. MMB GEOMETRY, MATERIAL AND FRACTURE DATA FROM THE PAPER
% ========================================================================

L  = 102.0e-3;
B  = 25.4e-3;
hT = 2*1.56e-3;
xLoad = L/2;

nLayers    = 2;
layerAngle = [0 0];
layerThick = [hT/2 hT/2];
zInterface = [-hT/2; 0; hT/2];

mat.E11  = 122.70e9;
mat.E22  = 10.10e9;
mat.E33  = 10.10e9;
mat.G12  = 5.50e9;
mat.G13  = 5.50e9;
mat.G23  = 3.70e9;
mat.nu12 = 0.25;
mat.nu13 = 0.25;
mat.nu23 = 0.25;

GIc  = 0.969e3;
GIIc = 1.719e3;
etaBK = 2.284;

switch mixedMode.case
    case "20_percent"
        targetModeRatio = 0.20;
        a0 = 22.7e-3;
        measuredGc = 1.103e3;
    case "50_percent"
        targetModeRatio = 0.50;
        a0 = 34.1e-3;
        measuredGc = 1.131e3;
    case "80_percent"
        targetModeRatio = 0.80;
        a0 = 31.4e-3;
        measuredGc = 1.376e3;
    otherwise
        error('Unknown mixedMode.case.');
end

fprintf('MMB case: %s\n',mixedMode.case);
fprintf('  target GII/GT = %.3f\n',targetModeRatio);
fprintf('  a0            = %.4f mm\n',1e3*a0);
fprintf('  measured Gc   = %.4f J/m^2\n\n',measuredGc);

%% ========================================================================
% 3. SINGLE- OR MULTIDOMAIN CHEBYSHEV-LOBATTO GRID
% ========================================================================
%
% Both choices use nodal physical fields and direct differentiation of the
% discrete fields. No conventional FE stiffness matrix is assembled.
%
% The validated MMB implementation uses one global domain.

Nz = 5;

switch grid.mode
    case "single_domain"
        Nx = grid.single.Nx;
        [xiX,DxiX,wXiX] = chebyshevLobattoClenshawCurtis(Nx);
        x  = 0.5*L*(xiX+1);
        Dx = (2/L)*DxiX;
        wx = 0.5*L*wXiX;
        domainMap = {(1:Nx).'};

    case "multidomain"
        [x,Dx,wx,domainMap] = multidomainChebyshevLobatto( ...
            grid.multi.breaks,grid.multi.nodes);
        Nx = numel(x);

        if min(abs(grid.multi.breaks-a0)) > 100*eps(max(L,1))
            error('For multidomain mode, a0 must be a domain boundary.');
        end
        if min(abs(grid.multi.breaks-xLoad)) > 100*eps(max(L,1))
            error('For multidomain mode, L/2 must be a domain boundary.');
        end

    otherwise
        error('Unknown grid.mode. Use "single_domain" or "multidomain".');
end

if any(wx <= 0)
    error('Non-positive axial quadrature weight detected.');
end
if any(diff(x) <= 0)
    error('The axial grid must be strictly increasing.');
end

fprintf('Chebyshev-Lobatto grid mode: %s\n',grid.mode);
fprintf('  Unique x nodes       = %d\n',Nx);
fprintf('  Number of domains    = %d\n',numel(domainMap));
fprintf('  Nearest a0-node error   = %.6e micrometres\n', ...
    1e6*min(abs(x-a0)));
fprintf('  Nearest midpoint error  = %.6e micrometres\n\n', ...
    1e6*min(abs(x-xLoad)));

[zeta,~,wZeta] = chebyshevLobattoClenshawCurtis(Nz);

zLayer  = cell(nLayers,1);
wzLayer = cell(nLayers,1);
for k = 1:nLayers
    hk = layerThick(k);
    zMid = 0.5*(zInterface(k)+zInterface(k+1));
    zLayer{k}  = zMid + 0.5*hk*zeta;
    wzLayer{k} = 0.5*hk*wZeta;
end

%% ========================================================================
% 4. MATERIAL MATRICES
% ========================================================================

Ccell = cell(nLayers,1);
Scell = cell(nLayers,1);
for k = 1:nLayers
    Ccell{k} = transformedPlaneStrainApproximation(layerAngle(k),mat);
    Scell{k} = inv(Ccell{k});
end

% Constitutive consistency checks
for k = 1:nLayers
    if norm(Ccell{k}-Ccell{k}.','fro')/max(norm(Ccell{k},'fro'),eps) > 1e-12
        error('Layer %d constitutive matrix is not symmetric.',k);
    end
    [~,pchol] = chol(Ccell{k});
    if pchol ~= 0
        error('Layer %d constitutive matrix is not positive definite.',k);
    end
end

%% ========================================================================
% 5. HEALTHY CSLT AND FULL DISCONTINUITY UNKNOWNS
% ========================================================================

nInterface = nLayers+1;
nQnodeH = 2*nInterface + 4*nLayers;
nQH = Nx*nQnodeH;

% Candidate discontinuity variables at every x node.
% Only nodes up to and including the constrained trial tip are activated.
nUDCandidate = Nx;
nWDCandidate = Nx;

localWBottom = dof('W',1,nLayers);
localWTop    = dof('W',nInterface,nLayers);

%% ========================================================================
% 6. INDEPENDENT STRESS VARIABLES
% ========================================================================

nPoint = nLayers*Nx*Nz;
nSigma = 3*nPoint;

fprintf('ODIL-XCSLT mixed-mode MMB nodewise formulation\n');
fprintf('-------------------------------------------------\n');
fprintf('Nx                           = %d\n',Nx);
fprintf('Nz per sublaminate           = %d\n',Nz);
fprintf('Healthy CSLT unknowns        = %d\n',nQH);
fprintf('Candidate axial jumps        = %d\n',nUDCandidate);
fprintf('Candidate opening jumps      = %d\n',nWDCandidate);
fprintf('Independent stresses         = %d\n\n',nSigma);

%% ========================================================================
% 7. BUILD EXTENDED STRAIN OPERATORS
% ========================================================================
%
% epsilon = BH*qH + BUD*uD + BWD*wD
%
% Centered discontinuity enrichment retained intentionally:
%   lower sublaminate sign = -1/2, PsiD = Psi3
%   upper sublaminate sign = +1/2, PsiD = Psi1
%
% At the interface:
%   uGamma+ - uGamma- = uD
%   wGamma+ - wGamma- = wD

BH  = spalloc(nSigma,nQH,nPoint*3*Nx*10);
BUD = spalloc(nSigma,Nx,nPoint*3*Nx*4);
BWD = spalloc(nSigma,Nx,nPoint*3*Nx*4);
pointWeight = zeros(nPoint,1);
Sblock = spalloc(nSigma,nSigma,9*nPoint);

for k = 1:nLayers
    hk = layerThick(k);

    for ix = 1:Nx
        for iz = 1:Nz
            p = pointIndex(k,ix,iz,Nx,Nz);
            rows = 3*(p-1)+(1:3);

            [B0H,BxH,B0UD,BxUD,B0WD,BxWD] = ...
                localExtendedStrainMatrices( ...
                zeta(iz),hk,k,nLayers,nQnodeH);

            % Healthy direct contribution
            for d = 1:nQnodeH
                v = B0H(:,d);
                if any(v ~= 0)
                    BH(rows,qHIndex(ix,d,Nx)) = ...
                        BH(rows,qHIndex(ix,d,Nx)) + v;
                end
            end

            % Healthy x derivatives
            for jx = 1:Nx
                Dij = Dx(ix,jx);
                if abs(Dij) > 0
                    for d = 1:nQnodeH
                        v = Dij*BxH(:,d);
                        if any(v ~= 0)
                            BH(rows,qHIndex(jx,d,Nx)) = ...
                                BH(rows,qHIndex(jx,d,Nx)) + v;
                        end
                    end
                end
            end

            % Axial discontinuity direct term
            if any(B0UD ~= 0)
                BUD(rows,ix) = BUD(rows,ix) + B0UD;
            end

            % Opening direct term
            if any(B0WD ~= 0)
                BWD(rows,ix) = BWD(rows,ix) + B0WD;
            end

            % Discontinuity x derivatives
            for jx = 1:Nx
                Dij = Dx(ix,jx);
                if abs(Dij) > 0
                    vUD = Dij*BxUD;
                    vWD = Dij*BxWD;
                    if any(vUD ~= 0)
                        BUD(rows,jx) = BUD(rows,jx) + vUD;
                    end
                    if any(vWD ~= 0)
                        BWD(rows,jx) = BWD(rows,jx) + vWD;
                    end
                end
            end

            pointWeight(p) = B*wx(ix)*wzLayer{k}(iz);
            Sblock(rows,rows) = Scell{k};
        end
    end
end

Wdiag = repelem(pointWeight,3);
W = spdiags(Wdiag,0,nSigma,nSigma);

%% ========================================================================
% 8. APPLIED MIDPOINT DISPLACEMENT HISTORY
% ========================================================================

if ~isempty(loadControl.explicitHistory)
    displacementSteps = loadControl.explicitHistory(:);
else
    switch analysisMode
        case "fixed_crack"
            displacementSteps = loadControl.fixedDisplacement;
        case "propagation"
            displacementSteps = linspace(0,loadControl.maxDisplacement, ...
                loadControl.nSteps).';
        otherwise
            error('Unknown analysisMode.');
    end
end

if any(displacementSteps < 0)
    error('Use positive displacement magnitudes.');
end

%% ========================================================================
% 9. INITIAL NODEWISE CRACK DESCRIPTION
% ========================================================================

% The crack-tip node remains constrained (closed), while nodes strictly
% behind it are released. Multidomain mode places it exactly at a0.
% Single-domain mode uses the nearest Lobatto node and reports the error.
if grid.mode == "multidomain"
    tipNodeInitial = find(abs(x-a0) <= 100*eps(max(L,1)),1,'first');
    if isempty(tipNodeInitial)
        error('The multidomain grid does not contain the exact initial crack tip.');
    end
else
    [~,tipNodeInitial] = min(abs(x-a0));
end
if tipNodeInitial <= 1 || tipNodeInitial >= Nx
    error('The initial crack tip must be an interior node.');
end

lastReleased = tipNodeInitial-1;
initialLastReleased = lastReleased;

fprintf('Requested initial crack tip = %.9f mm\n',1e3*a0);
fprintf('Discrete initial crack tip  = %.9f mm\n',1e3*x(tipNodeInitial));
fprintf('Crack-tip error             = %.6e micrometres\n\n', ...
    1e6*(x(tipNodeInitial)-a0));

if timing.enabled
    timing.setupSeconds = toc(timing.setupTimer);
else
    timing.setupSeconds = NaN;
end

%% ========================================================================
% 10. CALIBRATE THE PHYSICAL RIGID-LEVER LENGTH
% ========================================================================
%
% A single external actuator displacement is imposed through rigid-lever
% kinematics:
%
%   delta_P = alphaMid*w_mid + alphaEnd*w_end
%
%   alphaMid = 1 + 2*c/L
%   alphaEnd = -2*c/L
%
% where c is the distance from the midpoint loading link to the external
% actuator point. The multiplier of this one constraint is the actuator
% force. No independent end and midpoint displacements are prescribed.

leverLength = mixedMode.initialLeverLength;

if mixedMode.calibrateLeverLength
    fprintf('Calibrating rigid-lever length for target GII/GT = %.4f...\n', ...
        targetModeRatio);

    calibrationArgs = {BH,BUD,BWD,W,Sblock,nQH,nSigma,Nx,Nz,nLayers, ...
        initialLastReleased,localWBottom,localWTop, ...
        mixedMode.referenceDisplacement,x,xLoad,L,B,newton,contact,vcct};

    cGrid = linspace(mixedMode.leverLengthLower, ...
        mixedMode.leverLengthUpper,mixedMode.scanPoints);
    betaGrid = nan(size(cGrid));

    for jc = 1:numel(cGrid)
        betaGrid(jc) = evaluateInitialModeRatioLever( ...
            cGrid(jc),calibrationArgs{:});
    end

    residualGrid = betaGrid-targetModeRatio;
    bracket = find(residualGrid(1:end-1).*residualGrid(2:end) <= 0,1,'first');

    if isempty(bracket)
        [~,iBest] = min(abs(residualGrid));
        leverLength = cGrid(iBest);
        warning(['No sign-changing lever-length bracket was found. ' ...
            'Using c = %.6f mm, which gives GII/GT = %.6f.'], ...
            1e3*leverLength,betaGrid(iBest));
    else
        cLo = cGrid(bracket);
        cHi = cGrid(bracket+1);
        betaLo = betaGrid(bracket);

        for ic = 1:mixedMode.maxCalibrationIter
            cMid = 0.5*(cLo+cHi);
            betaMid = evaluateInitialModeRatioLever( ...
                cMid,calibrationArgs{:});

            if mixedMode.verboseCalibration
                fprintf('  calibration %2d: c = %.8f mm, GII/GT = %.8f\n', ...
                    ic,1e3*cMid,betaMid);
            end

            leverLength = cMid;

            if abs(betaMid-targetModeRatio) <= mixedMode.ratioTolerance
                break;
            end

            if (betaLo-targetModeRatio)*(betaMid-targetModeRatio) <= 0
                cHi = cMid;
            else
                cLo = cMid;
                betaLo = betaMid;
            end
        end
    end
end

alphaEnd = -2*leverLength/L;
alphaMid = 1 + 2*leverLength/L;

fprintf('Rigid-lever length c     = %.8f mm\n',1e3*leverLength);
fprintf('Lever coefficient end    = %.8f\n',alphaEnd);
fprintf('Lever coefficient middle = %.8f\n\n',alphaMid);

%% ========================================================================
% 11. STORAGE
% ========================================================================

nSteps = numel(displacementSteps);
reactionHistory = nan(nSteps,1);
crackHistory    = nan(nSteps,1);
GIHistory       = nan(nSteps,1);
GIIHistory      = nan(nSteps,1);
GTHistory       = nan(nSteps,1);
actualCMODHistory = nan(nSteps,1);
meanMotionHistory = nan(nSteps,1);
upperArmHistory = nan(nSteps,1);
lowerArmHistory = nan(nSteps,1);

newtonFinalLoss       = nan(nSteps,1);
newtonFinalResidual   = nan(nSteps,1);
newtonFinalStep       = nan(nSteps,1);
newtonIterations      = nan(nSteps,1);
matrixCondEstimate    = nan(nSteps,1);
matrixStructuralRank  = nan(nSteps,1);
matrixSymmetryError   = nan(nSteps,1);
reactionImbalanceHist = nan(nSteps,1);
newtonSolveHistories  = cell(nSteps,1);

% Timing storage
timing.stepSeconds       = nan(nSteps,1);
timing.modelBuildSeconds = nan(nSteps,1);
timing.newtonSeconds     = nan(nSteps,1);
timing.postSeconds       = nan(nSteps,1);
timing.solvesPerStep     = nan(nSteps,1);

lastSolution = [];
finalPost = [];

%% ========================================================================
% 11. LOAD-STEP, CONTACT, AND NODEWISE PROPAGATION LOOP
% ========================================================================

contactNodes = [];
contactUpdateHistory = nan(nSteps,1);
contactCountHistory  = nan(nSteps,1);
midpointDisplacementHistory = nan(nSteps,1);
actuatorDisplacementHistory = nan(nSteps,1);
GIHistory = nan(nSteps,1);
GIIHistory = nan(nSteps,1);
GTHistory = nan(nSteps,1);
modeIIFractionHistory = nan(nSteps,1);
GcHistory = nan(nSteps,1);
BKIndexHistory = nan(nSteps,1);

for iStep = 1:nSteps
    stepTimer = tic;
    stepBuildSeconds = 0.0;
    stepNewtonSeconds = 0.0;
    stepPostSeconds = 0.0;

    deltaControl = displacementSteps(iStep);
    crackStable = false;
    solveCounter = 0;
    historiesThisStep = {};
    contactUpdatesThisStep = 0;

    while ~crackStable
        if lastReleased >= Nx-1
            warning('Crack tip has reached the final usable node.');
            crackStable = true;
            break;
        end

        tipNode = lastReleased+1;

        switch contact.mode
            case "closed_crack"
                contactNodes = (1:lastReleased).';
            case "none"
                contactNodes = zeros(0,1);
            case "active_set"
                contactNodes = contactNodes(contactNodes <= lastReleased);
            otherwise
                error('Unknown contact.mode.');
        end

        contactStable = false;
        while ~contactStable
            modelTimer = tic;
            model = buildMMBModel( ...
                BH,BUD,BWD,W,Sblock,nQH,nSigma,Nx,Nz,nLayers, ...
                lastReleased,localWBottom,localWTop,deltaControl, ...
                x,xLoad,L,contactNodes,B,leverLength);
            stepBuildSeconds = stepBuildSeconds + toc(modelTimer);

            if ~newton.warmStart
                lastSolution = [];
            end

            newtonTimer = tic;
            [solution,nr] = solveMixedSystemNewton( ...
                model,lastSolution,newton);
            stepNewtonSeconds = stepNewtonSeconds + toc(newtonTimer);
            lastSolution = solution;

            solveCounter = solveCounter+1;
            historiesThisStep{solveCounter} = nr; %#ok<SAGROW>

            postTimer = tic;
            post = recoverMMBSolution( ...
                solution,model,W,nQH,Nx,Nz,nLayers);
            stepPostSeconds = stepPostSeconds + toc(postTimer);

            if contact.mode == "active_set"
                releasedNodes = (1:lastReleased).';
                penetrationNodes = releasedNodes( ...
                    post.wD(releasedNodes) < -contact.penetrationTol);
                newContact = setdiff(penetrationNodes,contactNodes);

                if isempty(newContact)
                    contactStable = true;
                else
                    contactNodes = unique([contactNodes;newContact]);
                    contactUpdatesThisStep = contactUpdatesThisStep+1;
                    lastSolution = [];

                    if contact.verbose
                        fprintf('  Added %d contact constraints.\n', ...
                            numel(newContact));
                    end
                    if contactUpdatesThisStep >= contact.maxUpdates
                        warning('Maximum contact active-set updates reached.');
                        contactStable = true;
                    end
                end
            else
                contactStable = true;
            end
        end

        behindNode = tipNode-1;
        da1 = x(tipNode)-x(behindNode);
        if tipNode < Nx
            da2 = x(tipNode+1)-x(tipNode);
        else
            da2 = da1;
        end

        uBehind = post.uD(behindNode);
        wBehind = post.wD(behindNode);
        bBehind = post.slopeJump(behindNode);

        FxTip = post.tipSlidingReaction;
        FzTip = post.tipOpeningReaction;
        MxTip = post.tipRotationReaction;

        workII = FxTip*uBehind;
        workI  = FzTip*wBehind + MxTip*bBehind;

        GIIraw = -(1/(2*da2*B))*workII*(da2/da1);
        GIraw  = -(1/(2*da2*B))*workI *(da2/da1);

        switch vcct.signMode
            case "absolute_work"
                GII = abs(workII)/(2*da2*B)*(da2/da1);
                GI  = abs(workI )/(2*da2*B)*(da2/da1);
            case "signed"
                GII = max(GIIraw,0);
                GI  = max(GIraw,0);
            otherwise
                error('Unknown vcct.signMode.');
        end

        GT = GI+GII;
        modeIIFraction = GII/max(GT,eps);
        GcBK = GIc + (GIIc-GIc)*modeIIFraction^etaBK;
        BKindex = GT/max(GcBK,eps);

        if vcct.printDiagnostics
            fprintf('\nMode-II VCCT diagnostics at step %d\n',iStep);
            fprintf('  x_tip       = %.9f mm\n',1e3*x(tipNode));
            fprintf('  da1         = %.6e m\n',da1);
            fprintf('  da2         = %.6e m\n',da2);
            fprintf('  Fx*uD       = %+ .6e J\n',workII);
            fprintf('  Fz*wD+Mx*b  = %+ .6e J\n',workI);
            fprintf('  GIIraw      = %+ .6e J/m^2\n',GIIraw);
            fprintf('  GIraw       = %+ .6e J/m^2\n',GIraw);
        end

        fprintf(['Step %3d | midpoint = %8.4f mm | tip = %8.4f mm | ' ...
                 'GII = %9.3f | GI = %9.3f | Gc = %9.3f | Q = %9.3f N | ' ...
                 'contact = %d | NR = %d | relR = %.3e\n'], ...
                 iStep,1e3*post.midpointDisplacement,1e3*x(tipNode), ...
                 GII,GI,GcBK,post.reactionForce,numel(contactNodes), ...
                 nr.nIterations,nr.finalRelResidual);

        if GT > 0 && abs(modeIIFraction-targetModeRatio) > vcct.modeRatioTolerance
            warning('Mode ratio %.3f differs from target %.3f at step %d.', ...
                modeIIFraction,targetModeRatio,iStep);
        end

        if analysisMode == "propagation" && ...
                GT >= GcBK && tipNode < Nx-1
            lastReleased = tipNode;
            lastSolution = [];
        else
            crackStable = true;

            reactionHistory(iStep) = post.reactionForce;
            crackHistory(iStep) = x(tipNode);
            GIHistory(iStep) = GI;
            GIIHistory(iStep) = GII;
            GTHistory(iStep) = GT;
            modeIIFractionHistory(iStep) = modeIIFraction;
            GcHistory(iStep) = GcBK;
            BKIndexHistory(iStep) = BKindex;
            midpointDisplacementHistory(iStep) = ...
                post.midpointDisplacement;
            actuatorDisplacementHistory(iStep) = ...
                post.actuatorDisplacement;
            actualCMODHistory(iStep) = max(post.wD(1:lastReleased));
            meanMotionHistory(iStep) = 0;
            upperArmHistory(iStep) = post.wOuterUpper(1);
            lowerArmHistory(iStep) = post.wOuterLower(1);

            newtonFinalLoss(iStep) = nr.finalLoss;
            newtonFinalResidual(iStep) = nr.finalRelResidual;
            newtonFinalStep(iStep) = nr.finalRelStep;
            newtonIterations(iStep) = nr.nIterations;
            matrixCondEstimate(iStep) = nr.condEstimate;
            matrixStructuralRank(iStep) = nr.structuralRank;
            matrixSymmetryError(iStep) = nr.symmetryError;
            reactionImbalanceHist(iStep) = post.supportImbalance;
            newtonSolveHistories{iStep} = historiesThisStep;
            contactUpdateHistory(iStep) = contactUpdatesThisStep;
            contactCountHistory(iStep) = numel(contactNodes);

            timing.stepSeconds(iStep)       = toc(stepTimer);
            timing.modelBuildSeconds(iStep) = stepBuildSeconds;
            timing.newtonSeconds(iStep)     = stepNewtonSeconds;
            timing.postSeconds(iStep)       = stepPostSeconds;
            timing.solvesPerStep(iStep)     = solveCounter;

            finalPost = post;
        end
    end
end

reportedDisplacementHistory = actuatorDisplacementHistory;

if timing.enabled
    timing.totalSeconds = toc(timing.totalTimer);
else
    timing.totalSeconds = NaN;
end

%% ========================================================================
% 12. NUMERICAL CHECKS
% ========================================================================

fprintf('\nFinal MMB checks\n');
fprintf('--------------------\n');
fprintf('Final crack tip                 = %.8f mm\n', ...
    1e3*crackHistory(end));
fprintf('Final actuator displacement     = %.8f mm\n', ...
    1e3*actuatorDisplacementHistory(end));
fprintf('Final reaction force            = %.8f N\n', ...
    reactionHistory(end));
fprintf('Rigid-lever length c            = %.8f mm\n', ...
    1e3*leverLength);
fprintf('Final end-link force            = %.8f N\n', ...
    finalPost.endLeverForce);
fprintf('Final midpoint-link force       = %.8f N\n', ...
    finalPost.midLeverForce);
fprintf('Final support-force imbalance   = %.8e\n', ...
    finalPost.supportImbalance);
fprintf('Final lever kinematic error      = %.8e m\n', ...
    abs(finalPost.actuatorDisplacement-displacementSteps(end)));
fprintf('Final GII                       = %.8f J/m^2\n', ...
    GIIHistory(end));
fprintf('Final GI                        = %.8f J/m^2\n', ...
    GIHistory(end));
fprintf('Final GII/GT                    = %.8f\n', ...
    modeIIFractionHistory(end));
fprintf('Final BK toughness              = %.8f J/m^2\n', ...
    GcHistory(end));
fprintf('Final GT/Gc                     = %.8f\n', ...
    BKIndexHistory(end));
fprintf('Final active contact nodes      = %d\n', ...
    contactCountHistory(end));
fprintf('Final Newton relative residual  = %.8e\n', ...
    newtonFinalResidual(end));
fprintf('Final condest(Kmix)             = %.8e\n', ...
    matrixCondEstimate(end));
fprintf('Final structural rank           = %.0f\n', ...
    matrixStructuralRank(end));
fprintf('Total runtime                   = %.6f s\n', ...
    timing.totalSeconds);

if newtonFinalResidual(end) > 10*newton.tolRelResidual
    warning('Final Newton residual exceeds the requested tolerance.');
end
if abs(modeIIFractionHistory(end)-targetModeRatio) > vcct.modeRatioTolerance
    warning('The final mode ratio differs from the target value.');
end
if any(finalPost.wD(1:max(lastReleased,1)) < -10*contact.penetrationTol)
    warning('Residual crack-face penetration remains.');
end

%% ========================================================================
% 13. MIXED-MODE RESULTS
% ========================================================================

figure('Color','w');
plot(1e3*reportedDisplacementHistory,reactionHistory,'LineWidth',1.8);
xlabel('Actuator displacement [mm]');
ylabel('Actuator force [N]');
title(sprintf('MMB response: target G_{II}/G_T = %.0f%%',100*targetModeRatio));
grid on; box on;

figure('Color','w');
plot(1e3*reportedDisplacementHistory,1e3*crackHistory,'LineWidth',1.7);
xlabel('Actuator displacement [mm]');
ylabel('Crack-tip coordinate [mm]');
title('Nodewise mixed-mode crack propagation');
grid on; box on;

figure('Color','w');
plot(1e3*reportedDisplacementHistory,GIIHistory,'LineWidth',1.7); hold on;
plot(1e3*reportedDisplacementHistory,GIHistory,'--','LineWidth',1.4);
plot(1e3*reportedDisplacementHistory,GTHistory,':','LineWidth',1.5);
plot(1e3*reportedDisplacementHistory,GcHistory,'-.','LineWidth',1.5);
xlabel('Actuator displacement [mm]');
ylabel('Energy-release rate [J/m^2]');
legend('G_{II}','G_I','G_T','G_c^{BK}','Location','best');
title('Mixed-mode VCCT and B-K toughness');
grid on; box on;

figure('Color','w');
plot(1e3*reportedDisplacementHistory,modeIIFractionHistory,'LineWidth',1.7);
xlabel('Actuator displacement [mm]');
ylabel('G_{II}/G_T');
yline(targetModeRatio,'--','Target');
ylim([0 1.05]);
title('Mixed-mode ratio');
grid on; box on;

%% ========================================================================
% 14. FINAL FIELDS
% ========================================================================

figure('Color','w');
plot(1e3*x,1e3*finalPost.wOuterUpper,'LineWidth',1.6); hold on;
plot(1e3*x,1e3*finalPost.wOuterLower,'--','LineWidth',1.5);
xline(1e3*crackHistory(end),'-.','Crack tip');
xline(1e3*xLoad,':','Loading point');
xlabel('x [mm]');
ylabel('Transverse displacement [mm]');
legend('Upper surface','Lower surface','Location','best');
title('Final MMB transverse displacement');
grid on; box on;

figure('Color','w');
plot(1e3*x,1e3*finalPost.uD,'LineWidth',1.7);
xline(1e3*crackHistory(end),'-.','Crack tip');
xlabel('x [mm]');
ylabel('Axial sliding jump u_D [mm]');
title('MMB interfacial sliding');
grid on; box on;

figure('Color','w');
plot(1e3*x,1e6*finalPost.wD,'LineWidth',1.7); hold on;
if ~isempty(finalPost.contactNodes)
    plot(1e3*x(finalPost.contactNodes), ...
        1e6*finalPost.wD(finalPost.contactNodes),'ko','MarkerSize',5);
end
xline(1e3*crackHistory(end),'-.','Crack tip');
xlabel('x [mm]');
ylabel('Opening jump w_D [micrometres]');
legend('w_D','Active contact','Location','best');
title('Crack-face opening and contact');
grid on; box on;

%% ========================================================================
% 15. NEWTON, CONTACT, AND TIMING
% ========================================================================

figure('Color','w');
semilogy(1:nSteps,newtonFinalResidual,'-o','LineWidth',1.4);
xlabel('Applied displacement step');
ylabel('Final relative residual');
title('Newton convergence per applied step');
grid on; box on;

figure('Color','w');
plot(1:nSteps,newtonIterations,'-o','LineWidth',1.4);
xlabel('Applied displacement step');
ylabel('Newton corrections');
title('Newton corrections per applied step');
grid on; box on;

figure('Color','w');
plot(1:nSteps,contactCountHistory,'LineWidth',1.5); hold on;
plot(1:nSteps,contactUpdateHistory,'--','LineWidth',1.4);
xlabel('Applied displacement step');
ylabel('Count');
legend('Active contact nodes','Contact-set updates','Location','best');
title('Contact active-set history');
grid on; box on;

if timing.enabled && timing.plotPerStep
    figure('Color','w');
    plot(1:nSteps,timing.stepSeconds,'LineWidth',1.5); hold on;
    plot(1:nSteps,timing.newtonSeconds,'--','LineWidth',1.4);
    xlabel('Applied displacement step');
    ylabel('Wall-clock time [s]');
    legend('Total step','Newton solves','Location','best');
    title('MMB solution timing');
    grid on; box on;

    figure('Color','w');
    plot(1:nSteps,cumsum(timing.stepSeconds,'omitnan'),'LineWidth',1.5);
    xlabel('Applied displacement step');
    ylabel('Cumulative time [s]');
    title('Cumulative MMB runtime');
    grid on; box on;
end

%% ========================================================================
% 16. FINAL INTERFACIAL STRESSES
% ========================================================================

if plotOptions.showStressProfiles
    figure('Color','w');
    plot(1e3*x,finalPost.interfaceTauLower/1e6,'LineWidth',1.5); hold on;
    plot(1e3*x,finalPost.interfaceTauUpper/1e6,'--','LineWidth',1.4);
    xline(1e3*crackHistory(end),'-.','Crack tip');
    xlabel('x [mm]');
    ylabel('tau_{xz} [MPa]');
    legend('Lower face','Upper face','Location','best');
    title('Interfacial shear stress');
    grid on; box on;

    figure('Color','w');
    plot(1e3*x,finalPost.interfaceSigmaZLower/1e6,'LineWidth',1.5); hold on;
    plot(1e3*x,finalPost.interfaceSigmaZUpper/1e6,'--','LineWidth',1.4);
    xline(1e3*crackHistory(end),'-.','Crack tip');
    xlabel('x [mm]');
    ylabel('sigma_z [MPa]');
    legend('Lower face','Upper face','Location','best');
    title('Interfacial normal stress');
    grid on; box on;
end

%% ========================================================================
% SAVE ALL FIGURES AS PNG
% ========================================================================

figs = findall(groot,'Type','figure');
[~,idxFig] = sort([figs.Number]);
figs = figs(idxFig);

for iFig = 1:numel(figs)
    ax = findall(figs(iFig),'Type','axes');
    figTitle = '';

    for iAx = 1:numel(ax)
        thisTitle = ax(iAx).Title.String;

        if iscell(thisTitle)
            thisTitle = strjoin(thisTitle,' ');
        end

        if ~isempty(thisTitle)
            figTitle = char(thisTitle);
            break;
        end
    end

    if isempty(figTitle)
        figTitle = sprintf('figure_%02d',iFig);
    end

    % Convert title to safe filename
    figTitle = regexprep(figTitle,'[\\/:*?"<>|]','');
    figTitle = regexprep(figTitle,'[{}_$^]','');
    figTitle = strrep(figTitle,'\','');
    figTitle = regexprep(strtrim(figTitle),'\s+','_');

    fileName = fullfile(resultsDir,[figTitle '.png']);
    exportgraphics(figs(iFig),fileName,'Resolution',300);
end

fprintf('\nSaved %d figures to:\n%s\n',numel(figs),resultsDir);
