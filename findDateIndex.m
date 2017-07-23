%% Find the index in "time" corresponding to the specified time range,
%  looking more widely if no index is an exact match and picking an index
%  near the center if more than one is found.
function I = findDateIndex(dateStr1, dateStr2, tArray)
    I = [];
    span = 0;
    while length(I) == 0
        I = find( datenum(dateStr1)-span < tArray & tArray < datenum(dateStr2)+span );
        span = span + 1;
        assert(span < 10, 'Not finding the requested date in a reasonable range!');
    end
    % For very small dt we may get more than one I value.  Pick one near
    % the middle.
    if length(I) > 1
        I = I(floor(length(I)/2));
    end
end