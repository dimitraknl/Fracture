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

localWBottom = dof('W',1,nLayers);
localWTop    = dof('W',nInterface,nLayers);

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

