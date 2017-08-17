%#codegen
% This function conforms to the rules for functions which are arguments to
% ode45.  The first argument is a scalar time value, and the second is a
% column vector of the starting values for the odes.  The return value is a
% column vector of the first derivatives of those values at the given time.
% 
% In order to keep things readable, the input vector will be broken into
% familiar variable names inside the function.  A good compiler will be
% able to remove the associated overhead, or if the result is slow the code
% can be rewritten after it is trusted.
%
% The docs recommend interpolating any time-dependent input values in the
% function.  tMonths is a zero-based array of months for finding values in
% any required array.
function [dydt] = odeFunction(t, startVals, tMonths, ...
                            temp, SelVx, C_seed, S_seed, OA, omegaFactor, con, MutVx)
                       
    Sn    = con.Sn;
    Cn    = con.Cn;
    KSx   = con.KSx;
    G     = con.G;
    KS    = con.KS;
    KC    = con.KC;
    A     = con.A;
    Mu    = con.Mu;
    um    = con.um;
    a     = con.a;
    b     = con.b;
    EnvVx = con.EnvVx;
    %%
    cols = Sn*Cn;
    Sold = startVals(1:cols)';
    Cold = startVals(cols+1:cols*2)';
    ri0 = startVals(cols*2+1:cols*3)'; % not used, but here so its place in the inputs is not forgotten.
    vgi0 = startVals(cols*3+1:cols*4)';
    gi0 = startVals(cols*4+1:end)'; 

    assert(length(Sold) == Sn*Cn, 'There should be one symbiont entry for each coral type * each symbiont type.');
    assert(length(Cold) == Sn*Cn, 'There should be one coral entry for each coral type * each symbiont type.');
    
    %%
    % Note that storage order with two symbionts per coral time is
    % massive1, branching1, massive2, branching2
    % where we'll typically treat "1" as the wild or common symbiont type
    % and "2" as an introduced or other alternate type.
    %%
    
    Sm = sum(Sold(1:2:end));  % Sum symbionts in massives
    Sb = sum(Sold(2:2:end));  % Sum symbionts in branching
    SA = [Sm Sb];
    % SAx = Sold(1, :);  % replaced 1/30/2017 with line below
    SAx = repmat(SA, 1, 2);
    
    C1 = [Cold(2) Cold(1)];   % [branch mass]
    C2 = [Cold(1) Cold(2)];   % [mass branch]
    
    % Start getting interpolated values not needed in the original RK
    % approach.
    T = interp1(tMonths, temp, t);
    rm  = a*exp(b*T);
    oFactor = interp1(tMonths, omegaFactor, t);
    G = [con.Gm con.Gb] * oFactor;
    
    ri = (1- (vgi0 + EnvVx + (min(0, gi0 - T)).^2) ./ (2*SelVx)) .* exp(b*min(0, T - gi0)) * rm;


    
    % Baskett 2009 equations 4 and 5.  k1 indicates the derivative at t sub i
    dSdT = Sold ./ (KSx .* Cold) .* (ri .* KSx .* Cold - rm .* SAx ) ;  %Change in symbiont pops %OK
    dCdT = (C2 .* (G .* SA./ (KS .* C2).* (KC-A .* C1-C2)./KC - Mu ./(1+um .* SA./(KS .* C2))) ); %Change in coral pops %OK
    
    % We can't set a seed value rigidly, but we can refuse to drop if below
    % the seed.
    dSdT = dSdT .* max(0, sign(Sold-S_seed)); % max part is one when above seed, zero otherwise
    dCdT = repmat(dCdT, 1, Sn) .* max(0, sign(Cold-repmat(C_seed, 1, Sn))); % max part is one when above seed, zero otherwise
    
    
    dgi  = ((vgi0 .* (T - gi0) ) ./ SelVx) * rm; % Change in Symbiont Mean Genotype
    dvgi = MutVx - (vgi0) .^2 ./ SelVx * rm ; % Change in Symbiont Mean Variance

    dridt = ri - ri0;
    dvgidt = dvgi;
    dgidt = dgi;
    dydt = [dSdT dCdT dridt dvgidt dgidt]';

end
