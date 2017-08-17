%#codegen
function [S, C, tOut, ri, gi, vgi, origEvolved] = tryDormandPrince(months, S0, C0, tMonths, ...
        temp, OA, omegaFactor, vgi, gi, MutVx, SelVx, C_seed, S_seed, suppressSuperIndex, ...
        superSeedFraction, oneShot, con)
    
    origEvolved = 0.0;
        
    % Note that S, C, gi, vgi, and others(?) are passed in a full arrays
    % but with only the first row initialized.  Since these can't be pre-computed,
    % they should now be treated
    % as additional ODEs with values at the time steps used by ode45,
    % rather than at fixed intervals.
    % 
    % ri is just zeros, and should be computed from vgi , gi, temp, and
    % constants.

        % Set up variables for using ode45.
        ri0 = zeros(1, con.Sn*con.Cn); % A byproduct of the calculation, but not an actual ODE of its own.
        inVar = [S0 C0 ri0 vgi(1, :) gi(1, :)]';
        % Time is in months in the equations, so the tspan input should be
        % in those units.
        % Output is a single column of time and multiple columns of
        % matching computed values.
        [tOut, yOut] = ode45(@(t, y) ...
            odeFunction(t, y, tMonths, temp, SelVx, C_seed, S_seed, OA, omegaFactor', con, MutVx), ...
            [0 months], inVar);
        
        cols = con.Sn*con.Cn;

        S = yOut(:, 1:cols);
        C = yOut(:, cols+1:cols*2);

        % XXX
        %    PROBLEM: super symbionts can't be introduced here - must be
        %    done in ode45 as a step function!
        % XXX

end