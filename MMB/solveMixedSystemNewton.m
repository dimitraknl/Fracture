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
