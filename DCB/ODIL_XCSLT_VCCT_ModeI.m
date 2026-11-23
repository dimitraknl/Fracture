%% ODIL-XCSLT NODEWISE VCCT: SINGLE- OR MULTIDOMAIN CHEBYSHEV-LOBATTO
% Mixed stress-displacement ODIL formulation for a DCB specimen.
% BAAAAAAA
%
% Main features
%   - Cubic-Hermite layerwise kinematics through the thickness
%   - Chebyshev-Lobatto collocation and Clenshaw-Curtis quadrature in x
%   - Independent stress unknowns: sigma_x, sigma_z, tau_xz
%   - Full discontinuity field: axial sliding uD(x), opening wD(x)
%   - Centered discontinuity enrichment (retained intentionally)
%   - Trial crack-tip constraints for uD, wD and relative rotation
%   - Nodewise VCCT using the actual nonuniform nodal distances
%   - Selectable single-arm or symmetric two-arm displacement control
%   - Direct Newton solution of the mixed residual equations
%   - ODIL loss used as a convergence diagnostic only
%
% ODIL residual and monitored loss
%   R(y)       = Kmix*y - rhs
%   L_ODIL(y)  = 0.5*R(y)'*R(y)
%
% Newton update
%   Kmix*dy = -R
%
% The fixed-crack problem is linear. With damping = 1, Newton should
% normally require one full correction after the initial residual.
%
% ODIL consistency:
%   - y contains discrete physical fields and constraint multipliers.
%   - Kmix*y-rhs is the stationarity residual of the discrete mixed
%     Hellinger-Reissner functional.
%   - L_ODIL = 0.5*R''*R is monitored, but no normal equations are formed.
%   - Newton solves the residual root directly: Kmix*dy = -R.
%   - Crack growth changes only the active discontinuity/constraint set.

clear; clc; close all;

% Results folder
resultsDir = fullfile(pwd,'results_mode1');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

%% ========================================================================
% 1. USER CONTROLS
% ========================================================================

analysisMode = "propagation";    % "fixed_crack" or "propagation"

% Axial collocation choice
%   "single_domain" : one global Chebyshev-Lobatto grid
%   "multidomain"   : fixed local Chebyshev-Lobatto subdomains with exact a0
grid.mode = "multidomain";

% Single-domain controls
grid.single.Nx = 51;

% Multidomain controls. a0 must be one of the domain boundaries.
grid.multi.breaks = [0,30,40,50,60,70,80,100]*1e-3;
grid.multi.nodes  = [7,5,5,5,5,5,7];

% Applied displacement control
% Loading configuration:
%   "single_arm" : upper arm prescribed by delta, lower arm held at zero.
%                  Plotted displacement = upper-arm displacement.
%   "two_arm"    : upper/lower arms prescribed by +delta/-delta.
%                  Plotted displacement = total CMOD = 2*delta.
loadControl.armConfiguration = "two_arm";
loadControl.fixedOpening = 1.0e-3;                  % [m], fixed-crack case
loadControl.maxOpening   = 6.0e-3;                  % [m], propagation case
loadControl.nSteps       = 45;
loadControl.explicitHistory = []; % e.g. [0;0.1;0.2;...]*1e-3

% Newton controls
newton.stopMode       = "tolerance"; % "tolerance" or "fixed"
newton.maxIter        = 8;           % safety cap
newton.minIter        = 1;
newton.fixedIter      = 2;           % used only if stopMode = "fixed"
newton.damping        = 1.0;
newton.tolRelResidual = 1.0e-4;
newton.tolRelStep     = 1.0e-12;
newton.verbose        = false;
newton.warmStart      = true;

% Matrix and residual diagnostics
newton.checkMatrix                 = true;
newton.failOnStructuralSingularity = true;
newton.conditionWarning            = 1.0e14;
newton.symmetryTolerance           = 1.0e-10;
newton.printBlockResiduals         = false;

% VCCT controls
vcct.useRotationalTerm = true;
vcct.signMode = "absolute_work"; % "absolute_work" or "signed"
vcct.printDiagnostics = false;

% Plot controls
plotOptions.figure7LikeStress = true;
plotOptions.maxNormalizedDistance = 3.0;

% Timing controls
timing.enabled       = true;
timing.printPerStep  = false;
timing.plotPerStep   = true;
timing.totalTimer    = tic;
timing.setupTimer    = tic;

%% ========================================================================
% 2. DCB GEOMETRY AND MATERIAL DATA FROM THE PAPER
% ========================================================================

L  = 100e-3;       % specimen length [m]
B  = 20e-3;        % specimen width [m]
hT = 3e-3;         % total thickness = 2h [m]
a0 = 30e-3;        % requested initial crack length [m]

nLayers    = 2;
layerAngle = [0 0];
layerThick = [hT/2 hT/2];
zInterface = [-hT/2; 0; hT/2];

mat.E11  = 135.30e9;
mat.E22  = 9.00e9;
mat.E33  = 9.00e9;
mat.G12  = 5.20e9;
mat.G13  = 5.20e9;
mat.G23  = 3.08e9;
mat.nu12 = 0.24;
mat.nu13 = 0.24;
mat.nu23 = 0.24;

GIc = 0.28e3;      % [J/m^2]

%% ========================================================================
% 3. SINGLE- OR MULTIDOMAIN CHEBYSHEV-LOBATTO GRID
% ========================================================================
%
% Both choices use nodal physical fields and direct differentiation of the
% discrete fields. No conventional FE stiffness matrix is assembled.
%
% In multidomain mode, a0 is a shared domain endpoint and is therefore
% represented exactly. In single-domain mode, the nearest Lobatto node is
% used and the geometric error is reported explicitly.

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
fprintf('  Nearest a0-node error= %.6e micrometres\n\n', ...
    1e6*min(abs(x-a0)));

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

localWBottom = dofW(1,nLayers);
localWTop    = dofW(nInterface,nLayers);

%% ========================================================================
% 6. INDEPENDENT STRESS VARIABLES
% ========================================================================

nPoint = nLayers*Nx*Nz;
nSigma = 3*nPoint;

fprintf('ODIL-XCSLT Chebyshev-Lobatto nodewise formulation\n');
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
% 8. APPLIED DISPLACEMENT HISTORY
% ========================================================================

if ~isempty(loadControl.explicitHistory)
    openingSteps = loadControl.explicitHistory(:);
else
    switch analysisMode
        case "fixed_crack"
            openingSteps = loadControl.fixedOpening;
        case "propagation"
            openingSteps = linspace(0,loadControl.maxOpening, ...
                loadControl.nSteps).';
        otherwise
            error('Unknown analysisMode.');
    end
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

nSteps = numel(openingSteps);
reactionHistory = nan(nSteps,1);
crackHistory    = nan(nSteps,1);
GIHistory       = nan(nSteps,1);
GIIHistory      = nan(nSteps,1);
GTHistory       = nan(nSteps,1);
actualCMODHistory = nan(nSteps,1);
meanMotionHistory = nan(nSteps,1);

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
% 11. LOAD-STEP AND NODEWISE PROPAGATION LOOP
% ========================================================================

for iStep = 1:nSteps
    stepTimer = tic;
    stepBuildSeconds = 0.0;
    stepNewtonSeconds = 0.0;
    stepPostSeconds = 0.0;

    deltaControl = openingSteps(iStep);
    crackStable = false;
    solveCounter = 0;
    historiesThisStep = {};

    while ~crackStable
        if lastReleased >= Nx-1
            warning('Crack tip has reached the final usable node.');
            crackStable = true;
            break;
        end

        tipNode = lastReleased+1;

        modelTimer = tic;
        model = buildFullFieldModel( ...
            BH,BUD,BWD,W,Sblock,nQH,nSigma,Nx,Nz,nLayers, ...
            lastReleased,localWBottom,localWTop,deltaControl, ...
            loadControl.armConfiguration,B);
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
        post = recoverFullFieldSolution( ...
            solution,model,W,nQH,Nx,Nz,nLayers);
        stepPostSeconds = stepPostSeconds + toc(postTimer);

        % ---------------------------------------------------------------
        % Nodewise VCCT using the actual nonuniform Chebyshev spacing.
        % Nodewise Eq. (18) of the paper for nonuniform Lobatto nodes.
        % Eq. (17) is the elementwise summation and is not used here.
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

        workItrans = FzTip*wBehind;
        workIrot   = MxTip*bBehind;
        workII     = FxTip*uBehind;

        if vcct.useRotationalTerm
            workI = workItrans + workIrot;
        else
            workI = workItrans;
        end

        GIraw = -(1/(2*da2*B))*workI*(da2/da1);
        GIIraw = -(1/(2*da2*B))*workII*(da2/da1);

        switch vcct.signMode
            case "absolute_work"
                GI  = abs(workI)/(2*da2*B)*(da2/da1);
                GII = abs(workII)/(2*da2*B)*(da2/da1);
            case "signed"
                GI  = max(GIraw,0);
                GII = max(GIIraw,0);
            otherwise
                error('Unknown vcct.signMode.');
        end
        GT = GI+GII;

        if vcct.printDiagnostics
            fprintf('\nVCCT diagnostics at step %d\n',iStep);
            fprintf('  x_tip       = %.9f mm\n',1e3*x(tipNode));
            fprintf('  da1         = %.6e m\n',da1);
            fprintf('  da2         = %.6e m\n',da2);
            fprintf('  Fz*wD       = %+ .6e J\n',workItrans);
            fprintf('  Mx*b        = %+ .6e J\n',workIrot);
            fprintf('  Fx*uD       = %+ .6e J\n',workII);
            fprintf('  GIraw       = %+ .6e J/m^2\n',GIraw);
            fprintf('  GIIraw      = %+ .6e J/m^2\n',GIIraw);
        end

        fprintf(['Step %3d | control = %8.4f mm | tip = %8.4f mm | ' ...
                 'GI = %9.3f | GII = %9.3f | F = %9.3f N | ' ...
                 'NR = %d | relR = %.3e\n'], ...
                 iStep,1e3*deltaControl,1e3*x(tipNode),GI,GII, ...
                 post.reactionForce,nr.nIterations,nr.finalRelResidual);

        if analysisMode == "propagation" && GI >= GIc && tipNode < Nx-1
            % Release exactly one node and resolve at the same displacement.
            lastReleased = tipNode;
            lastSolution = [];
        else
            crackStable = true;

            reactionHistory(iStep) = post.reactionForce;
            crackHistory(iStep) = x(tipNode);
            GIHistory(iStep) = GI;
            GIIHistory(iStep) = GII;
            GTHistory(iStep) = GT;
            actualCMODHistory(iStep) = post.actualCMOD;
            meanMotionHistory(iStep) = post.meanMouthMotion;

            newtonFinalLoss(iStep) = nr.finalLoss;
            newtonFinalResidual(iStep) = nr.finalRelResidual;
            newtonFinalStep(iStep) = nr.finalRelStep;
            newtonIterations(iStep) = nr.nIterations;
            matrixCondEstimate(iStep) = nr.condEstimate;
            matrixStructuralRank(iStep) = nr.structuralRank;
            matrixSymmetryError(iStep) = nr.symmetryError;
            reactionImbalanceHist(iStep) = post.reactionImbalance;
            newtonSolveHistories{iStep} = historiesThisStep;

            timing.stepSeconds(iStep)       = toc(stepTimer);
            timing.modelBuildSeconds(iStep) = stepBuildSeconds;
            timing.newtonSeconds(iStep)     = stepNewtonSeconds;
            timing.postSeconds(iStep)       = stepPostSeconds;
            timing.solvesPerStep(iStep)     = solveCounter;

            if timing.enabled && timing.printPerStep
                fprintf(['  Timing: total %.4f s | build %.4f s | ' ...
                         'Newton %.4f s | post %.4f s | solves %d\n'], ...
                    timing.stepSeconds(iStep),stepBuildSeconds, ...
                    stepNewtonSeconds,stepPostSeconds,solveCounter);
            end

            finalPost = post;
        end
    end
end

% Displacement reported on the horizontal axes:
%   single_arm -> actual displacement of the prescribed upper arm
%   two_arm    -> total relative opening between the two arms
reportedDisplacementHistory = reportedDisplacement( ...
    actualCMODHistory,meanMotionHistory, ...
    loadControl.armConfiguration);

if timing.enabled
    timing.totalSeconds = toc(timing.totalTimer);
else
    timing.totalSeconds = NaN;
end

%% ========================================================================
% 12. NUMERICAL CHECKS
% ========================================================================

fprintf('\nFinal checks\n');
fprintf('------------\n');
fprintf('Final crack tip             = %.8f mm\n',1e3*crackHistory(end));
fprintf('Final prescribed arm motion = %.8f mm\n',1e3*openingSteps(end));
fprintf('Final actual CMOD           = %.8f mm\n',1e3*actualCMODHistory(end));
fprintf('Final upper-arm displacement= %.8f mm\n', ...
    1e3*(meanMotionHistory(end)+0.5*actualCMODHistory(end)));
fprintf('Final lower-arm displacement= %.8f mm\n', ...
    1e3*(meanMotionHistory(end)-0.5*actualCMODHistory(end)));
fprintf('Final mean mouth motion     = %.6e mm\n',1e3*meanMotionHistory(end));
fprintf('Final reaction force        = %.8f N\n',reactionHistory(end));
fprintf('Final GI                    = %.8f J/m^2\n',GIHistory(end));
fprintf('Final GII                   = %.8f J/m^2\n',GIIHistory(end));
fprintf('Final ODIL relative residual= %.6e\n',newtonFinalResidual(end));
fprintf('Final matrix condest         = %.6e\n',matrixCondEstimate(end));
fprintf('Final structural rank        = %.0f / %.0f\n', ...
    matrixStructuralRank(end),finalPost.KmixSize(1));
fprintf('Final symmetry error         = %.6e\n',matrixSymmetryError(end));
fprintf('Final reaction imbalance     = %.6e\n',reactionImbalanceHist(end));

if timing.enabled
    fprintf('\nTiming summary\n');
    fprintf('--------------\n');
    fprintf('Setup/operator construction = %.6f s\n',timing.setupSeconds);
    fprintf('Load-step solution total    = %.6f s\n',sum(timing.stepSeconds,'omitnan'));
    fprintf('  Model assembly total      = %.6f s\n',sum(timing.modelBuildSeconds,'omitnan'));
    fprintf('  Newton solution total     = %.6f s\n',sum(timing.newtonSeconds,'omitnan'));
    fprintf('  Post-processing total     = %.6f s\n',sum(timing.postSeconds,'omitnan'));
    fprintf('Total script time           = %.6f s\n',timing.totalSeconds);
    fprintf('Mean accepted-step time     = %.6f s\n',mean(timing.stepSeconds,'omitnan'));
    fprintf('Mean Newton-solve time      = %.6f s\n', ...
        sum(timing.newtonSeconds,'omitnan')/max(sum(timing.solvesPerStep,'omitnan'),1));
end

%% ========================================================================
% 13. LOAD-DISPLACEMENT AND PROPAGATION PLOTS
% ========================================================================

figure('Color','w');
plot(1e3*reportedDisplacementHistory,reactionHistory,'LineWidth',1.6);
xlabel(displacementAxisLabel(loadControl.armConfiguration));
ylabel('Reaction force [N]');
title('ODIL-XCSLT Mode-I load-displacement response');
grid on; box on;

figure('Color','w');
plot(1e3*reportedDisplacementHistory,1e3*crackHistory,'LineWidth',1.6);
xlabel(displacementAxisLabel(loadControl.armConfiguration));
ylabel('Crack-tip coordinate [mm]');
title('Nodewise Mode-I crack propagation');
grid on; box on;

figure('Color','w');
plot(1e3*reportedDisplacementHistory,GIHistory,'LineWidth',1.6); hold on;
plot(1e3*reportedDisplacementHistory,GIIHistory,'--','LineWidth',1.4);
yline(GIc,'-.','G_{Ic}');
xlabel(displacementAxisLabel(loadControl.armConfiguration));
ylabel('Energy release rate [J/m^2]');
legend('G_I','G_{II}','G_{Ic}','Location','best');
title('Nodewise VCCT energy release rates');
grid on; box on;

%% ========================================================================
% 14. NEWTON CONVERGENCE PER APPLIED DISPLACEMENT STEP
% ========================================================================

figure('Color','w');
semilogy(1:nSteps,newtonFinalResidual,'-o','LineWidth',1.3);
xlabel('Applied displacement step');
ylabel('Final relative residual');
title('Newton convergence per applied displacement step');
grid on; box on;

figure('Color','w');
semilogy(1:nSteps,newtonFinalLoss,'-o','LineWidth',1.3);
xlabel('Applied displacement step');
ylabel('Final ODIL loss, 0.5 R^T R');
title('Final ODIL loss per applied displacement step');
grid on; box on;

figure('Color','w');
plot(1:nSteps,newtonIterations,'-o','LineWidth',1.3);
xlabel('Applied displacement step');
ylabel('Newton corrections');
title('Newton corrections per applied displacement step');
grid on; box on;


figure('Color','w');
semilogy(1:nSteps,matrixCondEstimate,'-o','LineWidth',1.2);
xlabel('Applied displacement step');
ylabel('condest(K_{mix})');
title('Mixed-system conditioning per accepted step');
grid on; box on;

figure('Color','w');
semilogy(1:nSteps,max(reactionImbalanceHist,eps),'-o','LineWidth',1.2);
xlabel('Applied displacement step');
ylabel('Relative reaction imbalance');
title('Symmetric DCB reaction-equilibrium check');
grid on; box on;

% Internal Newton history of the final accepted equilibrium solve
lastHistories = newtonSolveHistories{end};
if ~isempty(lastHistories)
    nrLast = lastHistories{end};

    figure('Color','w');
    semilogy(0:nrLast.nIterations-1,nrLast.loss,'-o','LineWidth',1.4);
    xlabel('Newton iteration');
    ylabel('ODIL loss, 0.5 R^T R');
    title('ODIL loss in the final equilibrium solve');
    grid on; box on;

    figure('Color','w');
    semilogy(0:nrLast.nIterations-1,nrLast.relResidual,'-o','LineWidth',1.4);
    hold on;
    semilogy(0:nrLast.nIterations-1,nrLast.relStep,'--s','LineWidth',1.3);
    xlabel('Newton iteration');
    ylabel('Relative norm');
    legend('Residual','Correction','Location','best');
    title('Newton convergence in the final equilibrium solve');
    grid on; box on;
end

%% ========================================================================
% 15. TIMING DIAGNOSTICS
% ========================================================================

if timing.enabled && timing.plotPerStep
    figure('Color','w');
    plot(1:nSteps,timing.stepSeconds,'LineWidth',1.4); hold on;
    plot(1:nSteps,timing.newtonSeconds,'--','LineWidth',1.3);
    plot(1:nSteps,timing.modelBuildSeconds,':','LineWidth',1.3);
    xlabel('Applied displacement step');
    ylabel('Wall-clock time [s]');
    legend('Total accepted step','Newton solves','Model assembly', ...
           'Location','best');
    title('Computational time per applied displacement step');
    grid on; box on;

    figure('Color','w');
    plot(1:nSteps,cumsum(timing.stepSeconds,'omitnan'),'LineWidth',1.5);
    xlabel('Applied displacement step');
    ylabel('Cumulative solution time [s]');
    title('Cumulative ODIL-XCSLT solution time');
    grid on; box on;

    figure('Color','w');
    plot(1:nSteps,timing.solvesPerStep,'-o','LineWidth',1.3);
    xlabel('Applied displacement step');
    ylabel('Equilibrium solves');
    title('Equilibrium solves per step, including crack releases');
    grid on; box on;
end

%% ========================================================================
% 16. TRANSVERSE DISPLACEMENT AND DISCONTINUITY FIELDS
% ========================================================================

figure('Color','w');
plot(1e3*x,1e3*finalPost.wOuterUpper,'LineWidth',1.7); hold on;
plot(1e3*x,1e3*finalPost.wOuterLower,'LineWidth',1.7);
plot(1e3*x,1e3*finalPost.wInterfaceHealthy,'--','LineWidth',1.3);
xline(1e3*crackHistory(end),'-.','Crack tip');
xlabel('Length [mm]');
ylabel('Transverse displacement [mm]');
legend('Upper external surface','Lower external surface', ...
       'Healthy interface field','Location','best');
title('Final transverse displacement along the specimen');
grid on; box on;

figure('Color','w');
plot(1e3*x,1e3*finalPost.wInterfaceUpper,'LineWidth',1.6); hold on;
plot(1e3*x,1e3*finalPost.wInterfaceLower,'LineWidth',1.6);
xline(1e3*crackHistory(end),'-.','Crack tip');
xlabel('Length [mm]');
ylabel('Interface displacement [mm]');
legend('Upper crack face','Lower crack face','Location','best');
title('Delamination-interface transverse displacement');
grid on; box on;

figure('Color','w');
plot(1e3*x,1e3*finalPost.wD,'LineWidth',1.6); hold on;
plot(1e3*x,1e3*finalPost.uD,'--','LineWidth',1.5);
xline(1e3*crackHistory(end),'-.','Crack tip');
xlabel('Length [mm]');
ylabel('Displacement jump [mm]');
legend('Opening w_D','Sliding u_D','Location','best');
title('Nodewise discontinuity fields');
grid on; box on;

%% ========================================================================
% 17. FIGURE-7-LIKE NORMALIZED INTERLAMINAR NORMAL STRESS
% ========================================================================
%
% The paper's Fig. 7 is a separate isotropic DCB stress benchmark. The plot
% below has the same normalization and axes style, but uses the present
% composite DCB state. A direct Fig. 7 validation requires a separate run
% with E = 70 GPa, nu = 0.3, B = 1 m and a0/h = 4.

if plotOptions.figure7LikeStress
    currentTip = crackHistory(end);
    hArm = hT/2;
    distanceFromTip = x-currentTip;
    ahead = distanceFromTip >= -100*eps(L);

    sigmaInterface = 0.5*( ...
        finalPost.interfaceSigmaZUpper + ...
        finalPost.interfaceSigmaZLower);

    Fnorm = max(abs(reactionHistory(end)),eps);
    normalizedStress = B*hArm/Fnorm*sigmaInterface;
    normalizedDistance = distanceFromTip/hArm;

    keep = ahead & normalizedDistance <= plotOptions.maxNormalizedDistance;

    figure('Color','w');
    plot(normalizedDistance(keep),normalizedStress(keep), ...
        'o-','LineWidth',1.4,'MarkerSize',4);
    xlabel('Distance from crack tip, (x-a)/h');
    ylabel('Normalized normal stress, Bh\sigma_3/F');
    title('Interlaminar normal stress ahead of the crack tip');
    xlim([0 plotOptions.maxNormalizedDistance]);
    grid on; box on;
end

%% ========================================================================
% 18. INTERFACE STRESS DIAGNOSTICS
% ========================================================================

figure('Color','w');
plot(1e3*x,finalPost.interfaceSigmaZUpper,'LineWidth',1.5); hold on;
plot(1e3*x,finalPost.interfaceSigmaZLower,'--','LineWidth',1.5);
xline(1e3*crackHistory(end),'-.','Crack tip');
xlabel('Length [mm]');
ylabel('Interface \sigma_z [Pa]');
legend('Upper face','Lower face','Location','best');
title('Interfacial normal stress');
grid on; box on;

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
    figTitle = regexprep(strtrim(figTitle),'\\s+','_');

    fileName = fullfile(resultsDir,[figTitle '.png']);
    exportgraphics(figs(iFig),fileName,'Resolution',300);
end

fprintf('\nSaved %d figures to:\n%s\n',numel(figs),resultsDir);

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

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

    idBUpperBottom = qHIndex(tipNode,dofB1(2,nLayers),Nx);
    idBLowerTop    = qHIndex(tipNode,dofB2(1,nLayers),Nx);
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
        qHIndex(Nx,dofU(2,nLayers),Nx)];

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

function [solution,nr] = solveMixedSystemNewton(model,initialGuess,newton)

    K = model.Kmix;
    rhs = model.rhs;
    n = size(K,1);

    if size(K,2) ~= n || numel(rhs) ~= n
        error('The mixed system must be square and compatible with rhs.');
    end
    if any(~isfinite(nonzeros(K))) || any(~isfinite(rhs))
        error('Non-finite entries detected in Kmix or rhs.');
    end

    % Check each distinct crack-state matrix only once. The matrix does not
    % depend on the prescribed displacement magnitude for a fixed crack.
    persistent checkedKeys checkedValues
    if isempty(checkedKeys)
        checkedKeys = strings(0,1);
        checkedValues = struct('rank',{},'cond',{},'sym',{});
    end
    key = sprintf('%d_%d_%d',n,nnz(K),model.tipNode);
    hit = find(checkedKeys == string(key),1);

    if isempty(hit)
        structuralRank = sprank(K);
        symmetryError = norm(K-K.','fro')/max(norm(K,'fro'),eps);
        if newton.checkMatrix
            try
                condEstimate = condest(K);
            catch
                condEstimate = Inf;
            end
        else
            condEstimate = NaN;
        end

        if structuralRank < n
            msg = sprintf('Kmix is structurally rank deficient: %d/%d.', ...
                structuralRank,n);
            if newton.failOnStructuralSingularity
                error(msg);
            else
                warning(msg);
            end
        end
        if symmetryError > newton.symmetryTolerance
            warning('Kmix symmetry error %.3e exceeds tolerance %.3e.', ...
                symmetryError,newton.symmetryTolerance);
        end
        if isfinite(condEstimate) && condEstimate > newton.conditionWarning
            warning('Kmix is poorly conditioned: condest = %.3e.',condEstimate);
        elseif ~isfinite(condEstimate)
            warning('Kmix condition estimate is non-finite.');
        end

        checkedKeys(end+1,1) = string(key);
        checkedValues(end+1) = struct('rank',structuralRank, ...
            'cond',condEstimate,'sym',symmetryError);
    else
        structuralRank = checkedValues(hit).rank;
        condEstimate = checkedValues(hit).cond;
        symmetryError = checkedValues(hit).sym;
    end

    if isempty(initialGuess) || numel(initialGuess) ~= numel(rhs)
        solution = zeros(size(rhs));
    else
        solution = initialGuess;
    end

    loss = nan(newton.maxIter,1);
    relResidual = nan(newton.maxIter,1);
    relStep = nan(newton.maxIter,1);
    nUsed = 0;

    for it = 1:newton.maxIter
        R = K*solution-rhs;
        loss(it) = 0.5*real(R.'*R);
        relResidual(it) = norm(R)/max(norm(rhs),eps);

        lastwarn('');
        dY = -K\R;
        [warnMsg,warnId] = lastwarn; %#ok<ASGLU>
        if any(~isfinite(dY))
            error('Newton correction contains non-finite values.');
        end
        if ~isempty(warnMsg) && contains(lower(warnMsg),'singular')
            warning('Linear solve warning: %s',warnMsg);
        end

        relStep(it) = norm(dY)/max(norm(solution),1);

        if newton.verbose
            fprintf(['    Newton %2d | loss = %.6e | relR = %.6e | ' ...
                     'relStep = %.6e\n'], ...
                     it,loss(it),relResidual(it),relStep(it));
        end

        solution = solution+newton.damping*dY;
        nUsed = it;

        switch newton.stopMode
            case "tolerance"
                % relResidual is evaluated before the current correction.
                % The final post-update residual is checked after the loop.
                if it >= newton.minIter && ...
                   relResidual(it) <= newton.tolRelResidual && ...
                   relStep(it) <= newton.tolRelStep
                    break;
                end
            case "fixed"
                if it >= newton.fixedIter
                    break;
                end
            otherwise
                error('Unknown newton.stopMode.');
        end
    end

    Rfinal = K*solution-rhs;
    finalLoss = 0.5*real(Rfinal.'*Rfinal);
    finalRelResidual = norm(Rfinal)/max(norm(rhs),eps);

    blockNorm.q = norm(Rfinal(model.rows.q));
    blockNorm.constitutive = norm(Rfinal(model.rows.constitutive));
    blockNorm.interface = norm(Rfinal(model.rows.interface));
    blockNorm.tip = norm(Rfinal(model.rows.tip));

    if newton.printBlockResiduals
        fprintf('      ||R_q||       = %.6e\n',blockNorm.q);
        fprintf('      ||R_const||   = %.6e\n',blockNorm.constitutive);
        fprintf('      ||R_interface||= %.6e\n',blockNorm.interface);
        fprintf('      ||R_tip||     = %.6e\n',blockNorm.tip);
    end

    if newton.stopMode == "tolerance" && ...
       finalRelResidual > newton.tolRelResidual
        warning('Newton solve ended with relative residual %.3e.', ...
            finalRelResidual);
    end

    nr.loss = loss(1:nUsed);
    nr.relResidual = relResidual(1:nUsed);
    nr.relStep = relStep(1:nUsed);
    nr.nIterations = nUsed;
    nr.finalLoss = finalLoss;
    nr.finalRelResidual = finalRelResidual;
    nr.finalRelStep = relStep(nUsed);
    nr.blockNorm = blockNorm;
    nr.structuralRank = structuralRank;
    nr.condEstimate = condEstimate;
    nr.symmetryError = symmetryError;
end

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

function [B0H,BxH,B0UD,BxUD,B0WD,BxWD] = ...
    localExtendedStrainMatrices(zeta,hk,k,nLayers,nQnodeH)

    [Psi,dPsidz,psi1,dpsi1dz,psi3,dpsi3dz] = ...
        hermiteExtended(zeta,hk);

    B0H = zeros(3,nQnodeH);
    BxH = zeros(3,nQnodeH);

    iu1 = dofU(k,nLayers);
    iu2 = dofU(k+1,nLayers);
    iw1 = dofW(k,nLayers);
    iw2 = dofW(k+1,nLayers);
    ib1 = dofB1(k,nLayers);
    ib2 = dofB2(k,nLayers);
    id1 = dofD1(k,nLayers);
    id2 = dofD2(k,nLayers);

    uDofs = [iu1 iu2 ib1 ib2];
    wDofs = [iw1 iw2 id1 id2];

    BxH(1,uDofs) = Psi;
    B0H(2,wDofs) = dPsidz;
    B0H(3,uDofs) = dPsidz;
    BxH(3,wDofs) = Psi;

    if k == 1
        s = -0.5;
        psiD = psi3;
        dPsiDdz = dpsi3dz;
    elseif k == 2
        s = 0.5;
        psiD = psi1;
        dPsiDdz = dpsi1dz;
    else
        error('This implementation currently assumes two sublaminates.');
    end

    % uD contributions:
    %   epsilon_x = s*PsiD*uD,x
    %   gamma_xz  = s*PsiD,z*uD
    B0UD = [0;0;s*dPsiDdz];
    BxUD = [s*psiD;0;0];

    % wD contributions:
    %   epsilon_z = s*PsiD,z*wD
    %   gamma_xz  = s*PsiD*wD,x
    B0WD = [0;s*dPsiDdz;0];
    BxWD = [0;0;s*psiD];
end

function [Psi,dPsidz,psi1,dpsi1dz,psi3,dpsi3dz] = ...
    hermiteExtended(zeta,hk)

    psi1 = 0.25*(1-zeta)^2*(2+zeta);
    psi2 = 0.25*(1-zeta)^2*(1+zeta);
    psi3 = 0.25*(1+zeta)^2*(2-zeta);
    psi4 = 0.25*(1+zeta)^2*(zeta-1);

    dpsi1 = 0.25*(-2*(1-zeta)*(2+zeta)+(1-zeta)^2);
    dpsi2 = 0.25*(-2*(1-zeta)*(1+zeta)+(1-zeta)^2);
    dpsi3 = 0.25*( 2*(1+zeta)*(2-zeta)-(1+zeta)^2);
    dpsi4 = 0.25*( 2*(1+zeta)*(zeta-1)+(1+zeta)^2);

    Psi = [psi1 psi3 hk/2*psi2 hk/2*psi4];
    dPsidz = (2/hk)*[dpsi1 dpsi3 hk/2*dpsi2 hk/2*dpsi4];

    dpsi1dz = (2/hk)*dpsi1;
    dpsi3dz = (2/hk)*dpsi3;
end

function reported = reportedDisplacement(actualCMOD,meanMotion,armConfiguration)

    % Recover the two actual loaded-arm displacements.
    wTop = meanMotion + 0.5*actualCMOD;
    wBottom = meanMotion - 0.5*actualCMOD; %#ok<NASGU>

    switch armConfiguration
        case "single_arm"
            % One-arm loading: plot only the displacement of the moving arm.
            reported = abs(wTop);
        case "two_arm"
            % Two-arm loading: plot the complete relative opening.
            reported = abs(actualCMOD);
        otherwise
            error('Unknown armConfiguration.');
    end
end

function label = displacementAxisLabel(armConfiguration)
    switch armConfiguration
        case "single_arm"
            label = 'Upper-arm displacement [mm]';
        case "two_arm"
            label = 'Total crack-mouth opening [mm]';
        otherwise
            label = 'Applied displacement [mm]';
    end
end

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

function C = transformedPlaneStrainApproximation(theta,mat)

    if abs(theta) < 1e-12
        Ex = mat.E11;
        Ez = mat.E33;
        Gxz = mat.G13;
        nuxz = mat.nu13;
    elseif abs(abs(theta)-90) < 1e-12
        Ex = mat.E22;
        Ez = mat.E33;
        Gxz = mat.G23;
        nuxz = mat.nu23;
    else
        m = cosd(theta);
        n = sind(theta);
        Ex = 1/(m^4/mat.E11+n^4/mat.E22+ ...
            m^2*n^2*(1/mat.G12-2*mat.nu12/mat.E11));
        Ez = mat.E33;
        Gxz = mat.G13*m^2+mat.G23*n^2;
        nuxz = mat.nu13*m^2+mat.nu23*n^2;
    end

    nuzx = nuxz*Ez/Ex;
    den = 1-nuxz*nuzx;

    C = [Ex/den, nuxz*Ez/den, 0; ...
         nuzx*Ex/den, Ez/den, 0; ...
         0, 0, Gxz];
end

function [x,D,w] = chebyshevLobattoClenshawCurtis(N)

    if N < 2
        error('N must be at least 2.');
    end

    n = N-1;
    theta = pi*(0:n)'/n;
    xDesc = cos(theta);

    c = [2;ones(n-1,1);2].*(-1).^(0:n)';
    X = repmat(xDesc,1,N);
    dX = X-X.';

    Ddesc = (c*(1./c).')./(dX+eye(N));
    Ddesc = Ddesc-diag(sum(Ddesc,2));

    wDesc = zeros(N,1);
    ii = 2:n;
    v = ones(n-1,1);

    if mod(n,2)==0
        wDesc(1) = 1/(n^2-1);
        wDesc(N) = wDesc(1);
        for k = 1:n/2-1
            v = v-2*cos(2*k*theta(ii))/(4*k^2-1);
        end
        v = v-cos(n*theta(ii))/(n^2-1);
    else
        wDesc(1) = 1/n^2;
        wDesc(N) = wDesc(1);
        for k = 1:(n-1)/2
            v = v-2*cos(2*k*theta(ii))/(4*k^2-1);
        end
    end

    wDesc(ii) = 2*v/n;

    perm = N:-1:1;
    x = xDesc(perm);
    D = Ddesc(perm,perm);
    w = wDesc(perm);
end

function id = pointIndex(k,ix,iz,Nx,Nz)
    id = (k-1)*Nx*Nz+(iz-1)*Nx+ix;
end

function id = stressIndex(k,ix,iz,comp,Nx,Nz)
    p = pointIndex(k,ix,iz,Nx,Nz);
    id = 3*(p-1)+comp;
end

function id = qHIndex(ix,localDof,Nx)
    id = ix+(localDof-1)*Nx;
end

function id = dofU(j,~)
    id = j;
end

function id = dofW(j,nLayers)
    id = nLayers+1+j;
end

function id = dofB1(k,nLayers)
    id = 2*(nLayers+1)+2*(k-1)+1;
end

function id = dofB2(k,nLayers)
    id = 2*(nLayers+1)+2*(k-1)+2;
end

function id = dofD1(k,nLayers)
    id = 2*(nLayers+1)+2*nLayers+2*(k-1)+1;
end

function id = dofD2(k,nLayers)
    id = 2*(nLayers+1)+2*nLayers+2*(k-1)+2;
end


function [x,Dx,wx,domainMap] = multidomainChebyshevLobatto(breaks,nodesPerDomain)
% Build a fixed multidomain Chebyshev-Lobatto grid with shared interface
% nodes. Differentiation is local to each subdomain. At a shared interface,
% the left and right derivative rows are averaged for the smooth healthy
% fields. Clenshaw-Curtis weights from neighboring domains are summed.

    nDom = numel(nodesPerDomain);
    if numel(breaks) ~= nDom+1
        error('breaks must contain one more entry than nodesPerDomain.');
    end

    x = [];
    domainMap = cell(nDom,1);
    localData = cell(nDom,1);

    for e = 1:nDom
        Ne = nodesPerDomain(e);
        [xi,Di,wi] = chebyshevLobattoClenshawCurtis(Ne);
        xa = breaks(e); xb = breaks(e+1);
        xe = 0.5*(xb-xa)*(xi+1)+xa;
        De = (2/(xb-xa))*Di;
        we = 0.5*(xb-xa)*wi;

        if e == 1
            idx = (1:Ne).';
            x = xe;
        else
            idx = [(numel(x)); (numel(x)+(1:Ne-1)).'];
            x = [x; xe(2:end)]; %#ok<AGROW>
        end
        domainMap{e} = idx;
        localData{e} = struct('D',De,'w',we);
    end

    Nx = numel(x);
    Dx = sparse(Nx,Nx);
    wx = zeros(Nx,1);
    rowCount = zeros(Nx,1);

    for e = 1:nDom
        idx = domainMap{e};
        De = localData{e}.D;
        we = localData{e}.w;
        Dx(idx,idx) = Dx(idx,idx) + De;
        wx(idx) = wx(idx) + we;
        rowCount(idx) = rowCount(idx) + 1;
    end

    shared = find(rowCount > 1);
    for r = shared.'
        Dx(r,:) = Dx(r,:)/rowCount(r);
    end

    if any(diff(x) <= 0)
        error('The assembled multidomain grid is not strictly increasing.');
    end
end
