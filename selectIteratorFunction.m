
%% Select the appropriate version of the time iteration function.  Different
%  compiled MEX versions exist for common array sizes, and an uncompiled version
%  is used for anything else.  For now, the Mac computers get the uncompiled
%  code, but if we create compiled versions they can be inserted here.
function handle =  selectIteratorFunction(select, Computer)
    % "select" is the length of the time array, which is a proxy for
    % several other array sizes passed into the function.
    % Comments show the original case the code was compiled for, but any
    % case with the same array sizes is a match.
    
    
    sourceFile1 = dir('timeIteration.m');
    sourceDate1 = sourceFile1.datenum;
    sourceFile2 = dir('Runge_Kutta_2_min0_160914.m');
    sourceDate2 = sourceFile2.datenum;
    if Computer < 3; select = 0; end;  % Kludge so Macs work without compiled code.
    switch select
        case 23040
            % 1861 to 2100, dt = 0.125
            handle = @timeIteration_23040_mex;
            mexName = 'timeIteration_23040_mex';
        case 38400
            % control400, dt = 0.125
            handle = @timeIteration_38400_mex;
            mexName = 'timeIteration_38400_mex';
        case 19200
            % conrol400, dt = 0.25.
            handle = @timeIteration_400yr_mex;
            mexName = 'timeIteration_400yr_mex';
        case 11520
            % 1861 to 2100, dt = 0.25
            handle = @timeIteration_2100_mex;
            mexName = 'timeIteration_2100_mex';
        case 46080
            % 1861 to 2100, dt = 0.0625
            handle = @timeIteration_46080_mex;
            mexName = 'timeIteration_45080_mex';
        otherwise
            % No MEX file for this case, run the slow way.
            handle = @timeIteration;
            mexName = '';
    end
    if ~isempty(mexName)
        mexFile = dir(strcat(mexName,'.mexw64'));
        if isempty(mexFile)
            disp('Warning: expected compiled version not found.  Using uncompiled timeIteration.m');
            handle = @timeIteration;
        else
            mexDate = mexFile.datenum; 
            if (mexDate < sourceDate1) || (mexDate < sourceDate2)
                disp('Warning: timeIteration has not been recompiled.  Using slower uncompiled version.');
                handle = @timeIteration;
            end
        end
    end
    
end