%% ODIL-XCSLT NODEWISE VCCT: MODE II ENF DELAMINATION
% Mixed stress-displacement ODIL formulation for the End-Notch Flexure
% benchmark described in Siorikis et al. (2021).
%
% Main features
%   - Cubic-Hermite layerwise kinematics through the thickness
%   - Chebyshev-Lobatto collocation and Clenshaw-Curtis quadrature in x
%   - Independent stress unknowns: sigma_x, sigma_z, tau_xz
%   - Full discontinuity field: axial sliding uD(x), opening wD(x)
%   - Centered discontinuity enrichment, consistent with the Mode-I code
%   - Three-point-bending ENF displacement control
%   - Contact active set to prevent crack-face interpenetration
%   - Nodewise Mode-II VCCT using the actual nonuniform nodal distances
%   - Direct Newton solution of the mixed stationarity residual
%   - ODIL loss used as a convergence diagnostic only
%
% ODIL residual and monitored loss
%   R(y)       = Kmix*y - rhs
%   L_ODIL(y)  = 0.5*R(y)'*R(y)
%
% Newton update (no normal equations)
%   Kmix*dy = -R
%
% For a fixed crack and contact active set, the problem is linear.
% Nonlinearity enters only through crack-node release and contact updates.

clear; clc; close all;

% Results folder
resultsDir = fullfile(pwd,'results2');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

%% ========================================================================
% 1. USER CONTROLS
% ========================================================================

analysisMode = "propagation";    % "fixed_crack" or "propagation"

grid.mode = "single_domain";       % "single_domain" or "multidomain"
grid.single.Nx = 51;

% Include the initial crack tip and midpoint loading point exactly.
grid.multi.breaks = [0,29.21,40,50.8,65,85,101.6]*1e-3;
grid.multi.nodes  = [9,7,7,7,7,7];

% Positive downward midpoint displacement magnitudes.
loadControl.fixedDisplacement = 0.25e-3;
loadControl.maxDisplacement   = 3.0e-3;
loadControl.nSteps            = 15;
loadControl.explicitHistory   = [];

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
contact.mode           = "closed_crack";
contact.penetrationTol = 1.0e-10;
contact.maxUpdates     = 15;
contact.verbose        = false;

% VCCT controls
vcct.signMode = "absolute_work";   % "absolute_work" or "signed"
vcct.printDiagnostics = false;
vcct.modePurityWarning = 0.10;

% Crack-propagation controls
% The critical event is located using GII proportional to delta^2 for a
% fixed crack/contact state. The solution is then re-equilibrated at the
% same critical displacement after releasing one node.
propagation.useCriticalInterpolation = true;
propagation.relativeGTolerance       = 1.0e-5;
propagation.maxReleasesPerStep       = 40;
propagation.recordOnlyStableState    = true;

% Benchmark and mechanics checks
validation.checkExactGeometry        = true;
validation.checkBoundaryConditions   = true;
validation.bcTolerance               = 1.0e-10;
validation.forceBalanceTolerance     = 1.0e-8;
validation.minimumModeIIFraction     = 0.90;

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
% 2. ENF GEOMETRY AND MATERIAL DATA FROM THE PAPER
% ========================================================================

L  = 101.6e-3;
B  = 7.62e-3;
hT = 2*1.524e-3;
a0 = 29.21e-3;
xLoad = L/2;

nLayers    = 2;
layerAngle = [0 0];
layerThick = [hT/2 hT/2];
zInterface = [-hT/2; 0; hT/2];

mat.E11  = 148.00e9;
mat.E22  = 10.50e9;
mat.E33  = 10.50e9;
mat.G12  = 5.61e9;
mat.G13  = 5.61e9;
mat.G23  = 3.17e9;
mat.nu12 = 0.27;
mat.nu13 = 0.27;
mat.nu23 = 0.27;

GIIc = 0.55e3;

fprintf('\nENF benchmark used in this script\n');
fprintf('  L       = %.4f mm\n',1e3*L);
fprintf('  B       = %.4f mm\n',1e3*B);
fprintf('  2h      = %.4f mm\n',1e3*hT);
fprintf('  a0      = %.4f mm\n',1e3*a0);
fprintf('  load x  = %.4f mm\n',1e3*xLoad);
fprintf('  GIIc    = %.4f J/m^2\n\n',GIIc);

%% ========================================================================
% 3. SINGLE- OR MULTIDOMAIN CHEBYSHEV-LOBATTO GRID
% ========================================================================
%
% Both choices use nodal physical fields and direct differentiation of the
% discrete fields. No conventional FE stiffness matrix is assembled.
%
% In multidomain mode, both a0 and L/2 are shared domain endpoints.

Nz = 7;

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

if validation.checkExactGeometry && grid.mode == "multidomain"
    geomTol = 100*eps(max(L,1));
    if min(abs(x-a0)) > geomTol
        error('Exact initial crack location check failed.');
    end
    if min(abs(x-xLoad)) > geomTol
        error('Exact midpoint loading location check failed.');
    end
end

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

fprintf('ODIL-XCSLT Mode-II ENF nodewise formulation\n');
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
% 10. STORAGE
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
boundaryConditionErrorHist = nan(nSteps,1);

criticalEventDisplacement = [];
criticalEventForceBefore  = [];
criticalEventForceAfter   = [];
criticalEventCrackBefore  = [];
criticalEventCrackAfter   = [];
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
%
% For each requested midpoint displacement, the following nested sequence
% is used:
%   1. converge the contact active set;
%   2. evaluate nodewise VCCT;
%   3. if GII exceeds GIIc, return to the interpolated critical
%      displacement, release exactly one crack-tip node, and re-equilibrate;
%   4. repeat at the same displacement until the crack is stable;
%   5. return to the requested displacement and store only the final stable
%      state.

contactNodes = [];
contactUpdateHistory = nan(nSteps,1);
contactCountHistory  = nan(nSteps,1);
midpointDisplacementHistory = nan(nSteps,1);
GIHistory = nan(nSteps,1);
GIIHistory = nan(nSteps,1);
GTHistory = nan(nSteps,1);
modeIIFractionHistory = nan(nSteps,1);

for iStep = 1:nSteps
    stepTimer = tic;
    stepBuildSeconds = 0.0;
    stepNewtonSeconds = 0.0;
    stepPostSeconds = 0.0;

    deltaTarget = displacementSteps(iStep);
    deltaControl = deltaTarget;
    crackStableAtTarget = false;
    solveCounter = 0;
    historiesThisStep = {};
    contactUpdatesThisStep = 0;
    releasesThisStep = 0;
    criticalSubstepActive = false;
    pendingCriticalDelta = NaN;
    pendingCriticalForceBefore = NaN;
    pendingCriticalCrackBefore = NaN;

    while ~crackStableAtTarget
        if lastReleased >= Nx-1
            warning('Crack tip has reached the final usable node.');
            crackStableAtTarget = true;
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

        % ---------------------------------------------------------------
        % Converge contact for the current crack and displacement state
        % ---------------------------------------------------------------
        contactStable = false;
        while ~contactStable
            modelTimer = tic;
            model = buildModeIIModel( ...
                BH,BUD,BWD,W,Sblock,nQH,nSigma,Nx,Nz,nLayers, ...
                lastReleased,localWBottom,localWTop,deltaControl, ...
                x,xLoad,contactNodes,B);
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
            post = recoverModeIISolution( ...
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

        % ---------------------------------------------------------------
        % Nodewise VCCT using actual nonuniform nodal distances (Eq. 18)
        % ---------------------------------------------------------------
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

        fprintf(['Step %3d | target = %8.4f mm | solve = %8.4f mm | ' ...
                 'tip = %8.4f mm | GII = %9.3f | GI = %9.3f | ' ...
                 'F = %9.3f N | contact = %d | NR = %d | relR = %.3e\n'], ...
                 iStep,1e3*deltaTarget,1e3*post.midpointDisplacement, ...
                 1e3*x(tipNode),GII,GI,post.reactionForce, ...
                 numel(contactNodes),nr.nIterations,nr.finalRelResidual);

        if modeIIFraction < 1-vcct.modePurityWarning && GT > 0
            warning('Mode-II fraction is only %.3f at step %d.', ...
                modeIIFraction,iStep);
        end

        shouldPropagate = analysisMode == "propagation" && ...
            GII >= GIIc*(1-propagation.relativeGTolerance) && ...
            tipNode < Nx-1;

        if shouldPropagate
            releasesThisStep = releasesThisStep+1;
            if releasesThisStep > propagation.maxReleasesPerStep
                error('Exceeded propagation.maxReleasesPerStep at step %d.',iStep);
            end

            % First return to the critical displacement for this crack state.
            if propagation.useCriticalInterpolation && ...
                    ~criticalSubstepActive && deltaControl > 0 && GII > 0
                deltaCritical = deltaControl*sqrt(GIIc/GII);
                deltaCritical = max(0,min(deltaCritical,deltaControl));

                if deltaCritical < deltaControl*(1-1e-10)
                    pendingCriticalDelta = deltaCritical;
                    pendingCriticalForceBefore = ...
                        post.reactionForce*(deltaCritical/deltaControl);
                    pendingCriticalCrackBefore = x(tipNode);
                    deltaControl = deltaCritical;
                    criticalSubstepActive = true;
                    lastSolution = [];
                    continue;
                end
            end

            % Release one node and re-equilibrate at the same displacement.
            if ~isfinite(pendingCriticalDelta)
                pendingCriticalDelta = deltaControl;
                pendingCriticalForceBefore = post.reactionForce;
                pendingCriticalCrackBefore = x(tipNode);
            end

            lastReleased = tipNode;
            lastSolution = [];
            continue;
        end

        % ---------------------------------------------------------------
        % Crack is stable at the current substep displacement
        % ---------------------------------------------------------------
        if criticalSubstepActive
            % Save the stabilized post-release force at the same critical
            % displacement, then return to the requested load-step target.
            criticalEventDisplacement(end+1,1) = pendingCriticalDelta; %#ok<AGROW>
            criticalEventForceBefore(end+1,1)  = pendingCriticalForceBefore; %#ok<AGROW>
            criticalEventForceAfter(end+1,1)   = post.reactionForce; %#ok<AGROW>
            criticalEventCrackBefore(end+1,1)  = pendingCriticalCrackBefore; %#ok<AGROW>
            criticalEventCrackAfter(end+1,1)   = x(tipNode); %#ok<AGROW>

            if deltaControl < deltaTarget*(1-1e-12)
                deltaControl = deltaTarget;
                criticalSubstepActive = false;
                pendingCriticalDelta = NaN;
                pendingCriticalForceBefore = NaN;
                pendingCriticalCrackBefore = NaN;
                lastSolution = [];
                continue;
            end
        end

        % Store only the final stabilized state at the requested displacement.
        crackStableAtTarget = true;

        reactionHistory(iStep) = post.reactionForce;
        crackHistory(iStep) = x(tipNode);
        GIHistory(iStep) = GI;
        GIIHistory(iStep) = GII;
        GTHistory(iStep) = GT;
        modeIIFractionHistory(iStep) = modeIIFraction;
        midpointDisplacementHistory(iStep) = post.midpointDisplacement;
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
        boundaryConditionErrorHist(iStep) = post.boundaryConditionError;
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

reportedDisplacementHistory = midpointDisplacementHistory;

if timing.enabled
    timing.totalSeconds = toc(timing.totalTimer);
else
    timing.totalSeconds = NaN;
end

%% ========================================================================
% 12. NUMERICAL CHECKS
% ========================================================================

fprintf('\nFinal Mode-II checks\n');
fprintf('--------------------\n');
fprintf('Final crack tip                 = %.8f mm\n', ...
    1e3*crackHistory(end));
fprintf('Final midpoint displacement     = %.8f mm\n', ...
    1e3*midpointDisplacementHistory(end));
fprintf('Final reaction force            = %.8f N\n', ...
    reactionHistory(end));
fprintf('Final GII                       = %.8f J/m^2\n', ...
    GIIHistory(end));
fprintf('Final GI                        = %.8f J/m^2\n', ...
    GIHistory(end));
fprintf('Final GII/GT                    = %.8f\n', ...
    modeIIFractionHistory(end));
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
if modeIIFractionHistory(end) < validation.minimumModeIIFraction
    warning('The final response is not predominantly Mode II.');
end
if any(finalPost.wD(1:max(lastReleased,1)) < -10*contact.penetrationTol)
    warning('Residual crack-face penetration remains.');
end
if validation.checkBoundaryConditions && ...
        boundaryConditionErrorHist(end) > validation.bcTolerance
    warning('Boundary-condition error %.3e exceeds tolerance %.3e.', ...
        boundaryConditionErrorHist(end),validation.bcTolerance);
end
if reactionImbalanceHist(end) > validation.forceBalanceTolerance
    warning('Force-balance error %.3e exceeds tolerance %.3e.', ...
        reactionImbalanceHist(end),validation.forceBalanceTolerance);
end

%% ========================================================================
% 13. MODE-II RESULTS
% ========================================================================

figure('Color','w');
plot(1e3*reportedDisplacementHistory,reactionHistory,'LineWidth',1.8);
xlabel('Midpoint displacement [mm]');
ylabel('Reaction force [N]');
title('ENF Mode II: reaction force versus midpoint displacement');
grid on; box on;

if ~isempty(criticalEventDisplacement)
    figure('Color','w');
    plot(1e3*criticalEventDisplacement,criticalEventForceBefore,'-o', ...
        'LineWidth',1.4); hold on;
    plot(1e3*criticalEventDisplacement,criticalEventForceAfter,'--s', ...
        'LineWidth',1.4);
    xlabel('Critical midpoint displacement [mm]');
    ylabel('Reaction force [N]');
    legend('Before node release','After node release','Location','best');
    title('Interpolated nodewise Mode-II release events');
    grid on; box on;
end

figure('Color','w');
plot(1e3*reportedDisplacementHistory,1e3*crackHistory,'LineWidth',1.7);
xlabel('Midpoint displacement [mm]');
ylabel('Crack-tip coordinate [mm]');
title('Nodewise Mode-II crack propagation');
grid on; box on;

figure('Color','w');
plot(1e3*reportedDisplacementHistory,GIIHistory,'LineWidth',1.7); hold on;
plot(1e3*reportedDisplacementHistory,GIHistory,'--','LineWidth',1.4);
yline(GIIc,'-.','G_{IIc}');
xlabel('Midpoint displacement [mm]');
ylabel('Energy-release rate [J/m^2]');
legend('G_{II}','G_I','G_{IIc}','Location','best');
title('Nodewise VCCT components');
grid on; box on;

figure('Color','w');
plot(1e3*reportedDisplacementHistory,modeIIFractionHistory,'LineWidth',1.7);
xlabel('Midpoint displacement [mm]');
ylabel('G_{II}/G_T');
ylim([0 1.05]);
title('Mode-II fraction');
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
title('Final ENF transverse displacement');
grid on; box on;

figure('Color','w');
plot(1e3*x,1e3*finalPost.uD,'LineWidth',1.7);
xline(1e3*crackHistory(end),'-.','Crack tip');
xlabel('x [mm]');
ylabel('Axial sliding jump u_D [mm]');
title('Mode-II interfacial sliding');
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
semilogy(1:nSteps,max(boundaryConditionErrorHist,eps),'-o', ...
    'LineWidth',1.3); hold on;
semilogy(1:nSteps,max(reactionImbalanceHist,eps),'--', ...
    'LineWidth',1.3);
xlabel('Applied displacement step');
ylabel('Error');
legend('BC error [m]','Force-balance error','Location','best');
title('ENF mechanics checks');
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
    title('Mode-II solution timing');
    grid on; box on;

    figure('Color','w');
    plot(1:nSteps,cumsum(timing.stepSeconds,'omitnan'),'LineWidth',1.5);
    xlabel('Applied displacement step');
    ylabel('Cumulative time [s]');
    title('Cumulative Mode-II runtime');
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

    % Convert the plot title into a safe filename.
    figTitle = regexprep(figTitle,'[\\/:*?"<>|]','');
    figTitle = regexprep(figTitle,'[{}_$^]','');
    figTitle = strrep(figTitle,'\\','');
    figTitle = regexprep(strtrim(figTitle),'\s+','_');

    fileName = fullfile(resultsDir,[figTitle '.png']);
    exportgraphics(figs(iFig),fileName,'Resolution',300);
end

fprintf('\nSaved %d figures to:\n%s\n',numel(figs),resultsDir);

%% ========================================================================
