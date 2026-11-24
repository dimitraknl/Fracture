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
