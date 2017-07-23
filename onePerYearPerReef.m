%% Examine two lists of events and discard any event in one set which has a
%  matching event in the other set for the same year and reef number.  Then
%  combine the sets to produce on list with no year/reef duplicates.  Order
%  does not matter.
function [combined] = onePerYearPerReef(a, b)  
    % Split recursively if subsets are bigger than big.
    % Testing shows 32 is better than 16 or 64.
    big = 32;  
    la = length(a);
    lb = length(b);
    if la > big && lb > big 
        % Split both on the same reef boundary
        minK = min(min([a.k]), min([b.k]));
        maxK = max(max([a.k]), max([b.k]));
        if maxK > minK
            half = floor((minK + maxK)/2);
            c1 = onePerYearPerReef(a([a.k] <= half), b([b.k] <= half));
            c2 = onePerYearPerReef(a([a.k] >  half), b([b.k] >  half));
            combined = [c1 c2];
            return;
        end
    elseif la == 0
        combined = b;
        return;
    elseif lb == 0
        combined = a;
        return;
    end

    % This is O(n^2), and so very slow when a and b are large
    for i=1:la
        y = [a(i).year];
        k = [a(i).k];
        b([b.year] == y & [b.k] == k) = [];
    end
    combined = [a b];
end