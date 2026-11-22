function prepareFigure(fig, axList)
%PF.PREPAREFIGURE Light theme and no built-in axes interactions.
%   PF.PREPAREFIGURE(FIG, AXLIST) forces a light theme on FIG (the colours
%   of the apps assume one) and switches off the default pan/zoom
%   interactions and toolbars of the axes in AXLIST, which would otherwise
%   fight with the apps' own mouse handling.

    try
        theme(fig, 'light');
    catch
        % releases before R2025a have no figure themes
    end
    for ax = axList(:).'
        try
            disableDefaultInteractivity(ax);
            ax.Toolbar.Visible = 'off';
        catch
            % older releases have no default axes interactions
        end
    end
end
