function [ ] = saveExcelHistory( basePath, now, RCP, E, everyx, queueMax, ...
        elapsed, Bleaching_85_10_By_Event, bleachParams, pswInputs)
    % Append the history of this run to an Excel file.
    filename = strcat(basePath, 'RunSummaries.xlsx');
    sheet = 'Run Info';
    oldLines = 0;
    if exist(filename, 'file')
        [~, txt, ~] = xlsread(filename, sheet, 'A:A');
        oldLines = length(txt);
        unsure = 1;
        while unsure
            fid=fopen(filename,'a');
            if fid < 0
                % Construct a questdlg with three options
                choice = questdlg('RunSummaries.xlsx may be open in Excel.  Close it and try again to save results from this run.', ...
                    'File Open Conflict', ...
                    'Skip saving','Try again','Try again');
                % Handle response
                switch choice
                    case 'Skip saving'
                        disp([choice ' not saving.'])
                        unsure = 0;
                        return;
                    case 'Try again'
                        disp([choice ' re-checking file.'])
                end
            else
                fclose(fid);
                unsure = 0;
            end
        end
    end
    mat = {datestr(now), RCP, E, everyx, queueMax, elapsed, ...
            Bleaching_85_10_By_Event, ...
            bleachParams.sBleach(1), bleachParams.cBleach(1), ...
            bleachParams.sRecoverySeedMult(1), bleachParams.cRecoverySeedMult(1), ...
            bleachParams.cSeedThresholdMult(1), ...
            pswInputs(1), pswInputs(2), pswInputs(3), pswInputs(4)
        };
    if ~oldLines
        matHeader = ...
        {'Date', 'RCP', 'E', 'everyx', 'workers', 'run time', ...
            '85-2010 bleaching', ...
            'sBleach 1', 'cBleach 1', ...
            'sRecSeedMult 1', 'cRecSeedMult 1', ...
            'cSeedThreshMult 1', ...
            'pMin', 'pMax', 'exponent', 'divisor'
            };
        mat = [matHeader; mat];
    end
    range = strcat('A', num2str(oldLines+1));
    % In the absence of proper write-locking, try twice.  That should be
    % enough to get around most race conditions when the amount of time
    % spent writing these is tiny compared to the computation time.
    
    [status, message] = xlswrite(filename, mat, sheet, range);
    % success = 1.
    if status == 0
        fprintf('xls file could not be written.  Trying once more.\nError: %s', message);
        xlswrite(filename, mat, sheet, range);
    end
end

