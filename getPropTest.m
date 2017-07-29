% getPropTest selects the correct psw2 column based on a 5% bleaching
% target.  Choose the correct match according to the E and RCP values.
function [propTest] = getPropTest(E, RCP)
    propTest = 0;
    if E == 0
        if strcmp(RCP, 'rcp26') || strcmp(RCP, 'control400')
            propTest = 20;
        elseif strcmp(RCP, 'rcp85')
            propTest = 21;
        elseif strcmp(RCP, 'rcp45')
            propTest = 24;
        elseif strcmp(RCP, 'rcp60')
            propTest = 25;
        end
    elseif E == 1
        if strcmp(RCP, 'rcp26') || strcmp(RCP, 'control400')
            propTest = 22;
        elseif strcmp(RCP, 'rcp85')
            propTest = 23;
        elseif strcmp(RCP, 'rcp45')
            propTest = 26;
        elseif strcmp(RCP, 'rcp60')
            propTest = 27;
        end
    end
    if propTest == 0
        error('No propTest value is defined for E = %d and RCP = %s\n', E, RCP);
    end
    return;
end
