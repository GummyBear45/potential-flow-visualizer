function deleteTimer(~, ~, tmr)
%PF.DELETETIMER Figure DeleteFcn that stops and deletes the app's timer.
%   set(fig, 'DeleteFcn', {@pf.deleteTimer, tmr})
%   Lives in its own file on purpose: callbacks defined inside an app file
%   stop resolving once that file is edited, so a window opened from an
%   older version could not be closed cleanly.

    tmr = tmr(isvalid(tmr));
    if isempty(tmr), return; end
    stop(tmr);
    delete(tmr);
end
