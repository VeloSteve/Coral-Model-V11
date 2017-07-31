function growthRateFigure(fullDir, suffix, yearStr, k, temp, gi, vgi, ssi, con, SelVx)
    % To compare native and enhanced symbionts, take gi and vgi
    % values just past the start point ssi and the plot above and below
    % the two genotypes for massive only.
        % vgi - 2D array, the same size as S and C, with symbiont variance.
        % gi  - Symbiont mean genotype over time
    
    vg(1) = vgi(ssi+1, 1);
    vg(2) = vgi(ssi+1, 3);
    g(1) = gi(ssi+1, 1);
    g(2) = gi(ssi+1, 3);
    t0 = temp(ssi+1);
    %{
  Fancy T range selection - no longer used.
    tD = 2.0;
    tU = 1.0;
    tMinInner = min(t0, min(g));
    tMaxInner = max(t0, max(g));
    tMin = tMinInner - tD;
    tMax = tMaxInner + tU;
    temps = linspace(tMin, tMax, 50);
    %}
    tMin = 20;
    tMax = 32;
    
    % XXX - hardwired to 2035 min and max for paper - not automatic!
    if k == 610
        tMin = 20.796;
        tMax = 33.94;
    elseif k == 1463
        tMin = 29.097;
        tMax = 30.563;
    end
    points = 100;
    temps = linspace(tMin, tMax, points);
    %fprintf('tMin %5.1f g(1) %5.1f g(2) T %5.1f %5.1f tMax %5.1f\n', tMin, g(1), t0, g(2), tMax);
    % The equation in the loop is exactly the one in timeIteration, other
    % than variable naming.
    rates = NaN(points, 2);
    rates2009 = NaN(points, 2);
    for j = 1:points
        % Last term of Baskett 2009 eq. 3:
        rm  = con.a*exp(con.b*temps(j)) ; % maximum possible growth rate at optimal temp
        T = temps(j);
        % As used in Spring 2017 code:
        r = (1- (vg + con.EnvVx(1:2:3) + (min(0, g - T)).^2) ./ (2*SelVx(1:2:3))) .* exp(con.b*min(0, T - g)) * rm;
        %r2014 = (1- (vg + con.EnvVx(1:2:3) + (min(0, g - T)).^2) ./ (2*SelVx(1:2:3))) * rm ;% Prevents cold water bleaching
        % Baskett 2009 eq. 3
        r2009 = (1- (vg + con.EnvVx(1:2:3) + (g - T).^2) ./ (2*SelVx(1:2:3))) * rm ;% Prevents cold water bleaching
        %rEpp = rm;
        rates(j, :) = r;
        %rates2014(j, :) = r2014;
        rates2009(j, :) = r2009;
        %ratesEpp(j, :) = rEpp;
    end
    % yCutoff is the lower of the rate at tMinInner and tMaxInner (almost
    % always tMaxInner).
    %{
    yCutoff = 0.0;  % Capture the lowest growth rate inside the optima and actual T.
    for T = [tMinInner tMaxInner]
        rm  = con.a*exp(con.b*temps(j)) ; % maximum possible growth rate at optimal temp
        r = (1- (vg + con.EnvVx(1:2:3) + (min(0, g - T)).^2) ./ (2*SelVx(1:2:3))) .* exp(con.b*min(0, T - g)) * rm;
        yCutoff = min(yCutoff, min(r));
    end
    %}
    
    % The graphs are pretty flat near the optimum, so try to trim off 
    % uninteresting low negative parts of the growth (vertical) axis.  The
    % problem is to not trim parts that intersect the key temperatures.
    
    specs = {'-k', '-m', ':c', '--m', '-b', ':b', '-.k'};

    figHandle = figure(4000+k);
    set(figHandle, 'color', 'w', 'OuterPosition',[60 269 1000 783]);
    axes1 = axes;

    plot(temps, rates(:,1), specs{1}); %, gi(:,2)); %, gi(:,3), gi(:,4));
    hold on;
    plot(temps, rates2009(:,1), specs{5});
    plot([g(1) g(1)], [min(min(rates)) max(max(rates))], ':k');  % current optimum
    %plot([g(2) g(2)], [min(min(rates)) max(max(rates))], '-.k');
    plot([t0 t0], [min(min(rates)) max(max(rates))], '--k');  % current actual T
    t = sprintf('Growth rate vs T for Reef %d in %s', k, yearStr);
    title(t);
    xlabel('Temperature (C)');
    ylabel('Growth Rate');
    set(axes1,'FontSize',21);
    legend({'symbiont growth', 'Baskett 2009', 'Adapted T', 'Current T'}, ...
        'Location', 'best', 'FontSize',18);

    % Constrain the minimum y in some cases
    %{
    yl = ylim;
    if yl(1) < yCutoff
        yl(1) = min(-0.1, yCutoff);
    end
    ylim(yl);
    %}
    % Now just use a fixed y range
    ylim([-0.5 1]);
        
    hold off;
    print('-dpdf', '-r200', strcat(fullDir, 'GrowthCurve', suffix, '.pdf'));
    savefig(strcat(fullDir, 'GrowthCurve', suffix, '.fig'));
end

%{
To find T range for a time period
Find the numerical date:
datenum('2035-01-01')
datenum('2035-12-31')
then manually look up that value in TIME.  The indexes there will be the
indexes in SST for that date range.

To shade background for a T range:
ty = [-0.5 1 1 -0.5];
tx = [20.796 20.796 33.94 33.94];
patch('DisplayName', '2034 T Range', 'YData', ty, 'XData', tx, 'FaceAlpha', 0.2, 'LineStyle', 'none', 'FaceColor', [.8 .8 .8])
%}