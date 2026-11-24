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
