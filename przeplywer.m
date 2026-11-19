function varargout = przeplywer
%PRZEPLYWER Interactive 2D potential-flow sandbox.
%   PRZEPLYWER opens a window with a uniform stream flowing left to right.
%   Sources, sinks, point vortices and doublets can be added, dragged,
%   strengthened and removed; the velocity vectors, streamlines and
%   stagnation points are recomputed live.
%
%   Mouse
%     left-drag          move an element
%     scroll             change strength of the element under the cursor
%     Shift + scroll     rotate a doublet axis
%     right-click        delete an element
%     right-click        panel: position, strength, doublet angle, blinking
%   Blinking: any element can jump between two positions every half of the
%   period "Blink T"; "Blinking vortex demo" sets up Aref's blinking vortex
%   with two-colour dye to show chaotic advection.
%   Keys (element appears under the cursor)
%     s source   k sink   v vortex CCW   c vortex CW   d doublet
%     Delete     remove the selected element
%
%   Everything is derived from the complex velocity
%       u - i v = dW/dz = U + sum A_k/(z - z_k) + sum B_k/(z - z_k)^2
%   with A_k = (m_k - i Gamma_k)/(2 pi) and B_k = -kappa_k e^{i alpha_k}/(2 pi).
%   The velocity is single valued, unlike the stream function (branch cut
%   behind every source) and the potential (branch cut behind every vortex).
%
%   FIG = PRZEPLYWER returns the figure handle.

    % ---- constants ------------------------------------------------------
    XL = [-4 4];  YL = [-2.5 2.5];      % domain
    DG = 0.25;                          % vector grid spacing
    H = 0.04;  NSTEP = 360;             % streamline step and step count
    NLEFT = 26;                         % streamlines entering from the left
    RSEED = 0.07;                       % seed radius around sources/sinks
    RKILL = 0.045;                      % streamlines end this close to an element
    PICK  = 0.22;                       % mouse pick radius
    SRC = 1;  VOR = 2;  DBL = 3;
    STEP = [0.25 0.25 0.25];            % scroll step per element type
    NPART = 700;                        % tracer particles
    DT = 0.03;                          % flow time advanced per animation frame
    MAXHOP = 0.1;                       % longest particle move per frame
    NDYE = 6000;                        % dye particles in the blinking demo

    % ---- state ----------------------------------------------------------
    U  = 1;                             % freestream speed
    ez = zeros(1,0);                    % element positions (complex)
    et = zeros(1,0);                    % element types
    es = zeros(1,0);                    % strengths: m, Gamma or kappa
    ea = zeros(1,0);                    % doublet axis angle
    eb = false(1,0);                    % element blinks between two positions
    ez2 = zeros(1,0);                   % the position it is NOT at right now
    eOn = zeros(1,0);                   % half period the element is in (0/1)
    TB = 4;                             % blink period (flow time)
    tFlow = 0;                          % flow time, advanced by the timer
    dyeMode = false;                    % particles coloured by origin, no respawn
    A  = zeros(1,0);  B = zeros(1,0);   % pole coefficients, see header
    sel = 0;  dragIdx = 0;  dragOff = 0;  busy = false;
    menuIdx = 0;                        % element shown in the properties panel
    pz = zeros(1,0);                    % particle positions (complex)
    pAge = zeros(1,0);  pLife = zeros(1,0);

    % ---- figure and axes ------------------------------------------------
    fig = figure('Name','Przeplywer - potential flow','NumberTitle','off', ...
        'Color','w','MenuBar','none','ToolBar','none', ...
        'Units','pixels','Position',[80 60 1180 800]);
    try
        theme(fig,'light');             % colours below assume a light figure
    catch
        % releases before R2025a have no figure themes
    end
    ax = axes('Parent',fig,'Units','normalized','Position',[0.05 0.21 0.86 0.77], ...
        'XLim',XL,'YLim',YL,'DataAspectRatio',[1 1 1],'Box','on','Layer','top');
    hold(ax,'on');
    try
        disableDefaultInteractivity(ax);
        ax.Toolbar.Visible = 'off';
    catch
        % older releases have no default axes interactions
    end

    [gx,gy] = meshgrid(XL(1):DG:XL(2), YL(1):DG:YL(2));
    gz = gx + 1i*gy;
    np = {'HitTest','off','PickableParts','none'};
    hQ = quiver(ax,gx,gy,0*gx,0*gy,0,'Color',[0.45 0.49 0.58],'MaxHeadSize',0.6,np{:});
    hS = line('Parent',ax,'XData',NaN,'YData',NaN,'Color',[0.13 0.33 0.62],'LineWidth',1,np{:});
    hDiv = line('Parent',ax,'XData',NaN,'YData',NaN,'Color',[0.90 0.45 0.05],'LineWidth',2,np{:});
    hP = scatter(ax,NaN,NaN,14,0,'filled','MarkerEdgeColor',[0.15 0.15 0.15], ...
        'LineWidth',0.25,np{:});
    colormap(ax,turbo(256));
    hCb = colorbar(ax,'Units','normalized','Position',[0.925 0.26 0.014 0.67]);
    ylabel(hCb,'speed |V|');
    hAxis = line('Parent',ax,'XData',NaN,'YData',NaN,'Color','k','LineWidth',2,np{:});
    hBlinkLn = line('Parent',ax,'XData',NaN,'YData',NaN,'Color',[0.3 0.3 0.3], ...
        'LineStyle',':','LineWidth',1.2,np{:});
    hGhost = line('Parent',ax,'XData',NaN,'YData',NaN,'LineStyle','none','Marker','s', ...
        'MarkerSize',11,'MarkerEdgeColor',[0.3 0.3 0.3],'LineWidth',1.2,np{:});
    mk = @(m,c) line('Parent',ax,'XData',NaN,'YData',NaN,'LineStyle','none', ...
        'Marker',m,'MarkerSize',11,'MarkerFaceColor',c,'MarkerEdgeColor','k', ...
        'LineWidth',1,np{:});
    hSrc = mk('o',[0.86 0.20 0.18]);
    hSnk = mk('o',[0.15 0.45 0.90]);
    hVp  = mk('s',[0.20 0.68 0.32]);
    hVn  = mk('s',[0.60 0.32 0.78]);
    hDb  = mk('d',[0.25 0.25 0.25]);
    hStag = line('Parent',ax,'XData',NaN,'YData',NaN,'LineStyle','none','Marker','x', ...
        'MarkerSize',11,'LineWidth',2.2,'Color','k',np{:});
    hSel = line('Parent',ax,'XData',NaN,'YData',NaN,'LineStyle','none','Marker','o', ...
        'MarkerSize',22,'LineWidth',1.5,'Color',[0.2 0.2 0.2],np{:});
    hT = gobjects(1,0);                 % label pool

    % ---- controls -------------------------------------------------------
    names = {'+ Source','+ Sink','+ Vortex CCW','+ Vortex CW','+ Doublet','Clear'};
    cbs = {@(~,~) addEl(SRC,2), @(~,~) addEl(SRC,-2), @(~,~) addEl(VOR,3), ...
           @(~,~) addEl(VOR,-3), @(~,~) addEl(DBL,1.5), @(~,~) clearAll()};
    for b = 1:numel(names)
        uicontrol(fig,'Style','pushbutton','String',names{b},'Units','normalized', ...
            'Position',[0.05+(b-1)*0.101 0.125 0.095 0.045],'FontSize',10,'Callback',cbs{b});
    end
    hUtxt = uicontrol(fig,'Style','text','String','U = 1.00','Units','normalized', ...
        'Position',[0.665 0.123 0.075 0.035],'BackgroundColor','w','FontSize',10, ...
        'HorizontalAlignment','right');
    hUsl = uicontrol(fig,'Style','slider','Min',0,'Max',2,'Value',U,'Units','normalized', ...
        'Position',[0.75 0.128 0.22 0.035],'Callback',@(~,~) setU());
    addlistener(hUsl,'ContinuousValueChange',@(~,~) setU());

    chk = @(str,x,w) uicontrol(fig,'Style','checkbox','String',str,'Value',1, ...
        'Units','normalized','Position',[x 0.075 w 0.035],'BackgroundColor','w', ...
        'FontSize',10,'Callback',@(~,~) redraw());
    cVec  = chk('Velocity vectors',0.05,0.13);
    cStr  = chk('Streamlines',0.19,0.10);
    cStag = chk('Stagnation points + dividing streamlines',0.30,0.30);
    cPart = chk('Particles',0.61,0.10);
    set(cPart,'Callback',@(~,~) toggleParticles());
    uicontrol(fig,'Style','text','String','Blink T','Units','normalized', ...
        'Position',[0.715 0.071 0.055 0.035],'BackgroundColor','w','FontSize',10, ...
        'HorizontalAlignment','right');
    hTB = uicontrol(fig,'Style','edit','String',num2str(TB),'Units','normalized', ...
        'Position',[0.775 0.075 0.045 0.035],'BackgroundColor','w','FontSize',10, ...
        'Callback',@(~,~) setPeriod());
    uicontrol(fig,'Style','pushbutton','String','Blinking vortex demo','Units','normalized', ...
        'Position',[0.83 0.07 0.14 0.045],'FontSize',10,'Callback',@(~,~) blinkDemo());

    uicontrol(fig,'Style','text','Units','normalized','Position',[0.05 0.012 0.92 0.05], ...
        'BackgroundColor','w','ForegroundColor',[0.35 0.35 0.35],'FontSize',9, ...
        'HorizontalAlignment','left','String', ...
        ['Drag: move  |  Scroll: strength  |  Shift+Scroll: rotate doublet  |  ' ...
         'Right-click: properties & blinking  |  Keys: s source, k sink, ' ...
         'v/c vortex CCW/CW, d doublet, Del remove']);

    % properties panel, shown by right-clicking an element
    MW = 200;  MH = 253;
    hMenu = uipanel(fig,'Units','pixels','Position',[0 0 MW MH],'Visible','off', ...
        'Title','','FontWeight','bold','BackgroundColor',[0.96 0.96 0.96]);
    rowLbl = {'x','y','strength','angle [deg]','other x','other y'};
    rowPos = [1 2 3 4 6 7];             % row 5 holds the blink checkbox
    hLbl = gobjects(1,6);  hEd = gobjects(1,6);
    hBlink = uicontrol(hMenu,'Style','checkbox','String','Blink between two positions', ...
        'Units','pixels','Position',[8 MH-22-5*27 185 22], ...
        'BackgroundColor',[0.96 0.96 0.96],'Callback',@(~,~) toggleBlink());
    for row = 1:6
        yr = MH - 22 - rowPos(row)*27;
        hLbl(row) = uicontrol(hMenu,'Style','text','String',rowLbl{row},'Units','pixels', ...
            'Position',[8 yr-3 78 20],'HorizontalAlignment','left', ...
            'BackgroundColor',[0.96 0.96 0.96]);
        hEd(row) = uicontrol(hMenu,'Style','edit','Units','pixels', ...
            'Position',[90 yr 100 22],'BackgroundColor','w','Callback',@(~,~) applyMenu());
    end
    uicontrol(hMenu,'Style','pushbutton','String','Delete','Units','pixels', ...
        'Position',[8 5 88 24],'Callback',@(~,~) deleteFromMenu());
    uicontrol(hMenu,'Style','pushbutton','String','Close','Units','pixels', ...
        'Position',[102 5 88 24],'Callback',@(~,~) closeMenu());

    set(fig,'WindowButtonDownFcn',@onDown,'WindowButtonMotionFcn',@onMove, ...
        'WindowButtonUpFcn',@onUp,'WindowScrollWheelFcn',@onScroll, ...
        'WindowKeyPressFcn',@onKey);

    % programmatic access, used for testing and scripting
    setappdata(fig,'api',struct('add',@addEl,'clear',@clearAll,'redraw',@redraw, ...
        'setU',@apiSetU,'stagnation',@stagnation,'velocity',@conjVel, ...
        'tick',@tick,'particles',@getParticles,'elements',@getElements, ...
        'openMenu',@openMenu,'blinkDemo',@blinkDemo,'setBlink',@apiSetBlink, ...
        'setPeriod',@apiSetPeriod,'time',@getTime));

    redraw();

    % particle animation runs on a timer so it keeps going between mouse events
    respawn(true(1,NPART));
    pAge = floor(rand(1,NPART).*pLife);         % stagger the first generation
    tmr = timer('ExecutionMode','fixedSpacing','Period',0.03,'BusyMode','drop', ...
        'TimerFcn',@tick);
    set(fig,'DeleteFcn',@(~,~) killTimer());
    start(tmr);
    if nargout, varargout{1} = fig; end

    % =====================================================================
    %  Flow field
    % =====================================================================
    function w = conjVel(z)
        % Complex conjugate velocity u - i*v at the points z.
        w = U + zeros(size(z));
        for k = 1:numel(ez)
            d = z - ez(k);
            w = w + A(k)./d + B(k)./d.^2;
        end
    end

    function v = dirn(z)
        % Unit vector along the flow, as a complex number.
        w = conjVel(z);
        v = conj(w)./abs(w);
    end

    function [x,y] = trace(z0,sg)
        % Integrate streamlines from seeds z0 (row) by arc length with the
        % midpoint rule; sg = +1 follows the flow, -1 goes upstream.
        n = numel(z0);
        if n == 0, x = NaN; y = NaN; return; end
        Z = nan(NSTEP+2,n);             % last row stays NaN: line separator
        Z(1,:) = z0;  z = z0;  alive = true(1,n);
        for i = 1:NSTEP
            zm = z + 0.5*H*sg.*dirn(z);
            z  = z + H*sg.*dirn(zm);
            z(~alive) = NaN;
            Z(i+1,:) = z;
            alive = alive & real(z) >= XL(1)-1e-9 & real(z) <= XL(2) ...
                          & imag(z) >= YL(1) & imag(z) <= YL(2) ...
                          & ~any(abs(z - ez.') < RKILL,1);
            if ~any(alive), break; end
        end
        x = real(Z(:));  y = imag(Z(:));
    end

    function zs = stagnation()
        % Stagnation points: dW/dz = 0 multiplied through by all the
        % denominators is a polynomial, so ROOTS finds every one of them.
        zs = zeros(0,1);
        K = numel(ez);
        if K == 0, return; end
        owner = repelem(1:K, 1 + (et == DBL));   % doublets are double poles
        r = ez(owner);
        N = numel(r);
        P = U*poly(r);
        for k = 1:K
            if et(k) == DBL
                t = B(k)*poly(r(owner ~= k));
            else
                i1 = find(owner == k,1);
                t = A(k)*poly(r([1:i1-1, i1+1:N]));
            end
            P = P + [zeros(1,N+1-numel(t)), t];
        end
        if ~any(P), return; end
        zs = roots(P);
        ok = isfinite(zs) & real(zs) > XL(1) & real(zs) < XL(2) ...
           & imag(zs) > YL(1) & imag(zs) < YL(2) ...
           & all(abs(zs - ez) > 1e-3, 2);
        zs = zs(ok);
    end

    % =====================================================================
    %  Particles
    % =====================================================================
    function respawn(idx,atInflow)
        % Give the particles IDX a new position and a fresh lifetime.
        % Particles that left the domain re-enter with the freestream at the
        % left edge, which keeps the density even; the others (too old, or
        % swallowed by a sink) reappear anywhere, so that fluid coming out
        % of sources gets particles too.
        n = nnz(idx);
        x = XL(1) + diff(XL)*rand(1,n);
        if nargin > 1 && U > 0
            in = atInflow(idx);
            x(in) = XL(1) + U*DT*rand(1,nnz(in));
        end
        pz(idx) = x + 1i*(YL(1) + diff(YL)*rand(1,n));
        pAge(idx) = 0;
        pLife(idx) = 300 + floor(400*rand(1,n));
    end

    function t = getTime()
        t = tFlow;
    end

    function p = getParticles()
        p = pz;
    end

    function s = getElements()
        s = struct('z',ez,'type',et,'strength',es,'angle',ea, ...
            'blink',eb,'z2',ez2);
    end

    function tick(~,~)
        % One animation frame: advance flow time, move blinking elements,
        % advect the particles with the true velocity (midpoint rule) and
        % colour them by local speed, or by origin in dye mode.
        if ~isgraphics(fig) || busy, return; end
        try
            tFlow = tFlow + DT;
            updateBlink();
            if ~get(cPart,'Value'), return; end
            hop = @(d) d.*min(1, MAXHOP./max(abs(d),eps));   % cap the step length
            zm = pz + 0.5*hop(DT*conj(conjVel(pz)));
            v  = conj(conjVel(zm));
            dz = hop(DT*v);
            pz = pz + dz;
            if dyeMode
                % dye is never recycled, otherwise the mixing pattern is lost;
                % it only vanishes into sinks and doublets (vortices just spin it)
                sing = reshape(ez(et ~= VOR),[],1);
                gone = ~isfinite(pz) | any(abs(pz - sing) < max(RKILL,abs(dz)),1);
                pz(gone) = NaN;
                set(hP,'XData',real(pz),'YData',imag(pz));
                drawnow limitrate nocallbacks
                return
            end
            pAge = pAge + 1;
            out = real(pz) < XL(1) | real(pz) > XL(2) ...
                | imag(pz) < YL(1) | imag(pz) > YL(2);
            % A particle closer to an element than its own step would jump
            % over it and bounce around the pole forever, so it ends here.
            dead = out | ~isfinite(pz) | pAge > pLife ...
                 | any(abs(pz - ez.') < max(RKILL,abs(dz)),1);
            sp = abs(v);
            if any(dead)
                respawn(dead,out);
                sp(dead) = abs(conjVel(pz(dead)));
            end
            set(hP,'XData',real(pz),'YData',imag(pz),'CData',sp);
            drawnow limitrate nocallbacks
        catch err
            if isgraphics(fig), rethrow(err); end   % ignore a closing window
        end
    end

    function toggleParticles()
        % the timer keeps running either way, because it also drives blinking
        if get(cPart,'Value')
            set(hP,'Visible','on');
            if ~dyeMode, set(hCb,'Visible','on'); end
        else
            set([hP hCb],'Visible','off');
        end
    end

    % =====================================================================
    %  Blinking elements (Aref's blinking vortex, generalised)
    % =====================================================================
    function p = phase()
        % 0 in the first half of each blink period, 1 in the second
        p = mod(floor(tFlow/(TB/2)),2);
    end

    function updateBlink()
        % Swap a blinking element to its other position at every half period.
        % The flow is steady in between, so streamlines are redrawn only then.
        k = find(eb & eOn ~= phase());
        k(k == dragIdx) = [];                   % do not yank it from the mouse
        if isempty(k), return; end
        tmp = ez(k);  ez(k) = ez2(k);  ez2(k) = tmp;
        eOn(k) = phase();
        redraw();
    end

    function setBlink(k,on,z2)
        if nargin > 2
            ez2(k) = clampDom(z2);
        elseif ~isfinite(ez2(k))
            % default partner position 1.2 to the right, or left near the edge
            d = 1.2;
            if real(ez(k)) + d > XL(2), d = -d; end
            ez2(k) = clampDom(ez(k) + d);
        end
        eb(k) = on;
        eOn(k) = phase();
        redraw();
    end

    function toggleBlink()
        if menuIdx == 0, return; end
        setBlink(menuIdx, logical(get(hBlink,'Value')));
    end

    function apiSetBlink(k,on,z2)
        if nargin > 2, setBlink(k,on,z2); else, setBlink(k,on); end
    end

    function setPeriod()
        v = str2double(get(hTB,'String'));
        if isfinite(v) && v > 0
            % keep every element in its current half period
            p0 = phase();  TB = v;
            tFlow = (floor(tFlow/(TB/2)) + mod(p0 - phase(),2) + 1e-6)*TB/2;
        end
        set(hTB,'String',num2str(TB));
    end

    function apiSetPeriod(v)
        set(hTB,'String',num2str(v));
        setPeriod();
    end

    function blinkDemo()
        % Aref (1984): one vortex, switched between two positions every half
        % period. Each half period alone is plain rotation, yet alternating
        % the two stretches and folds the dye: chaotic advection.
        clearAll();
        set(hUsl,'Value',0);  U = 0;  set(hUtxt,'String','U = 0.00');
        tFlow = 0;
        ez = -0.6;  et = VOR;  es = 3;  ea = 0;
        eb = true;  ez2 = 0.6;  eOn = 0;
        sel = 0;
        % dye: a disc of fluid around both positions, upper half orange,
        % lower half blue
        r = 1.5*sqrt(rand(1,NDYE));  t = 2*pi*rand(1,NDYE);
        pz = r.*exp(1i*t);
        col = repmat([0.15 0.45 0.90],NDYE,1);
        col(imag(pz) > 0,:) = repmat([0.95 0.50 0.10],nnz(imag(pz) > 0),1);
        dyeMode = true;
        set(cPart,'Value',1);
        set(hP,'XData',real(pz),'YData',imag(pz),'CData',col,'SizeData',5, ...
            'MarkerEdgeColor','none','Visible','on');
        set(hCb,'Visible','off');
        redraw();
    end

    function leaveDye()
        if ~dyeMode, return; end
        dyeMode = false;
        pz = zeros(1,NPART);  pAge = zeros(1,NPART);  pLife = zeros(1,NPART);
        respawn(true(1,NPART));
        set(hP,'XData',real(pz),'YData',imag(pz),'CData',zeros(1,NPART), ...
            'SizeData',14,'MarkerEdgeColor',[0.15 0.15 0.15]);
        if get(cPart,'Value'), set(hCb,'Visible','on'); end
    end

    function killTimer()
        try
            stop(tmr);  delete(tmr);
        catch
            % timer already gone
        end
    end

    % =====================================================================
    %  Drawing
    % =====================================================================
    function redraw()
        A = zeros(size(ez));  B = A;
        A(et == SRC) = es(et == SRC)/(2*pi);
        A(et == VOR) = -1i*es(et == VOR)/(2*pi);
        B(et == DBL) = -es(et == DBL).*exp(1i*ea(et == DBL))/(2*pi);
        set(ax,'CLim',[0 2*max(U,0.5)]);        % particle colour scale

        % --- velocity vectors: length follows speed, capped near singularities
        if get(cVec,'Value')
            w = conjVel(gz);
            w(~isfinite(w)) = 0;
            sp = abs(w);  cap = 2*max(U,0.5);
            sc = 0.9*DG*min(sp,cap)/cap./max(sp,eps);
            set(hQ,'UData',real(w).*sc,'VData',-imag(w).*sc,'Visible','on');
        else
            set(hQ,'Visible','off');
        end

        % --- streamlines
        if get(cStr,'Value')
            dy = diff(YL)/NLEFT;
            if U > 1e-6
                z0 = XL(1) + 1i*(YL(1) + dy*((1:NLEFT) - 0.5));
            else
                z0 = zeros(1,0);        % nothing enters from the left without U
            end
            sg = ones(size(z0));
            for k = 1:numel(ez)
                switch et(k)
                    case SRC
                        % equal volume flux between neighbouring lines
                        if es(k) == 0, continue; end
                        n = min(max(round(abs(es(k))/(max(U,1e-6)*dy)),6),40);
                        zk = ez(k) + RSEED*exp(2i*pi*((1:n) - 0.5)/n);
                        sk = sign(es(k))*ones(1,n);
                    case VOR
                        zk = ez(k) + 1i*[0.25 0.5 0.8 0.25 0.5 0.8];
                        sk = [1 1 1 -1 -1 -1];
                    otherwise
                        % closed loops inside the dividing circle of radius a
                        a = sqrt(abs(es(k))/(2*pi*max(U,0.05)));
                        a = min(max(a,0.15),1.5);
                        p = 1i*exp(1i*ea(k))*a*[0.35 0.7 -0.35 -0.7];
                        zk = ez(k) + [p p];
                        sk = [1 1 1 1 -1 -1 -1 -1];
                end
                z0 = [z0 zk]; sg = [sg sk]; %#ok<AGROW>
            end
            [x,y] = trace(z0,sg);
            set(hS,'XData',x,'YData',y,'Visible','on');
        else
            set(hS,'Visible','off');
        end

        % --- stagnation points and the streamlines through them
        if get(cStag,'Value')
            zs = stagnation();
            set(hStag,'XData',real(zs),'YData',imag(zs),'Visible','on');
            if isempty(zs)
                set(hDiv,'XData',NaN,'YData',NaN);
            else
                % near a stagnation point u - i v ~ c*(z - zs); the flow
                % leaves along angle -arg(c)/2 and arrives perpendicular to it
                d = zs - ez;
                c = sum(-A./d.^2 - 2*B./d.^3, 2);
                e = 0.02*exp(-0.5i*angle(c)).';
                zr = zs.';
                [x,y] = trace([zr+e, zr-e, zr+1i*e, zr-1i*e], ...
                    [ones(1,2*numel(zr)), -ones(1,2*numel(zr))]);
                set(hDiv,'XData',x,'YData',y);
            end
            set(hDiv,'Visible','on');
        else
            set([hStag hDiv],'Visible','off');
        end

        drawMarkers();
    end

    function drawMarkers()
        pts = @(h,z) set(h,'XData',real(z),'YData',imag(z));
        pts(hSrc, ez(et == SRC & es >= 0));
        pts(hSnk, ez(et == SRC & es <  0));
        pts(hVp,  ez(et == VOR & es >= 0));
        pts(hVn,  ez(et == VOR & es <  0));
        pts(hDb,  ez(et == DBL));

        zd = ez(et == DBL);  ed = 0.24*exp(1i*ea(et == DBL));
        seg = [zd - ed; zd + ed; nan(size(zd))];
        set(hAxis,'XData',real(seg(:)),'YData',imag(seg(:)));

        % where blinking elements will jump next
        zb = ez(eb);  zg = ez2(eb);
        seg = [zb; zg; nan(size(zb))];
        set(hBlinkLn,'XData',real(seg(:)),'YData',imag(seg(:)));
        pts(hGhost, zg);

        K = numel(ez);
        for k = numel(hT)+1:K
            hT(k) = text(ax,0,0,'','FontSize',9,'Color','k','BackgroundColor',[1 1 1 0.75], ...
                'Margin',1,'Clipping','on',np{:});
        end
        for k = 1:K
            switch et(k)
                case SRC, s = sprintf('m = %+.2f',es(k));
                case VOR, s = sprintf('\\Gamma = %+.2f',es(k));
                otherwise
                    s = sprintf('\\kappa = %+.2f,  %.0f^\\circ',es(k), ...
                        mod(rad2deg(ea(k)),360));
            end
            set(hT(k),'Position',[real(ez(k))+0.14, imag(ez(k))+0.17, 0], ...
                'String',s,'Visible','on');
        end
        set(hT(K+1:end),'Visible','off');

        if sel >= 1 && sel <= K
            pts(hSel, ez(sel));
        else
            sel = 0;
            set(hSel,'XData',NaN,'YData',NaN);
        end
        syncMenu();
    end

    % =====================================================================
    %  Properties panel
    % =====================================================================
    function openMenu(k)
        % Show the panel for element K next to the mouse pointer.
        menuIdx = k;  sel = k;
        fp = get(fig,'Position');  cp = get(fig,'CurrentPoint');
        x = min(max(cp(1) + 14, 1), fp(3) - MW);
        y = min(max(cp(2) - MH - 6, 1), fp(4) - MH);
        set(hMenu,'Position',[x y MW MH],'Visible','on');
        drawMarkers();
    end

    function deleteFromMenu()
        removeEl(menuIdx);
    end

    function closeMenu()
        menuIdx = 0;
        set(hMenu,'Visible','off');
    end

    function syncMenu()
        % Refresh the panel fields from the element, e.g. while dragging.
        if menuIdx == 0, return; end
        if menuIdx > numel(ez), closeMenu(); return; end
        k = menuIdx;
        switch et(k)
            case SRC
                if es(k) >= 0, ttl = 'Source'; else, ttl = 'Sink'; end
                lbl = 'm';
            case VOR
                if es(k) >= 0, ttl = 'Vortex CCW'; else, ttl = 'Vortex CW'; end
                lbl = 'Gamma';
            otherwise
                ttl = 'Doublet';  lbl = 'kappa';
        end
        set(hMenu,'Title',ttl);
        set(hLbl(3),'String',lbl);
        set(hEd(1),'String',sprintf('%.3f',real(ez(k))));
        set(hEd(2),'String',sprintf('%.3f',imag(ez(k))));
        set(hEd(3),'String',sprintf('%.2f',es(k)));
        set(hEd(4),'String',sprintf('%.1f',mod(rad2deg(ea(k)),360)));
        if et(k) == DBL, vis = 'on'; else, vis = 'off'; end
        set([hLbl(4) hEd(4)],'Enable',vis);
        set(hBlink,'Value',eb(k));
        if isfinite(ez2(k))
            set(hEd(5),'String',sprintf('%.3f',real(ez2(k))));
            set(hEd(6),'String',sprintf('%.3f',imag(ez2(k))));
        else
            set(hEd(5:6),'String','');
        end
        if eb(k), vis = 'on'; else, vis = 'off'; end
        set([hLbl(5:6) hEd(5:6)],'Enable',vis);
    end

    function applyMenu()
        % Take the typed values; anything that is not a number is ignored
        % and the field snaps back to the current value.
        k = menuIdx;
        if k == 0 || k > numel(ez), return; end
        val = str2double(get(hEd,'String'));
        if isfinite(val(1)) && isfinite(val(2))
            ez(k) = clampDom(val(1) + 1i*val(2));
        end
        if isfinite(val(3)), es(k) = val(3); end
        if isfinite(val(4)), ea(k) = deg2rad(val(4)); end
        if isfinite(val(5)) && isfinite(val(6))
            ez2(k) = clampDom(val(5) + 1i*val(6));
        end
        redraw();
    end

    function tf = overMenu()
        % True when the mouse pointer is on the open properties panel.
        tf = false;
        if menuIdx == 0, return; end
        cp = get(fig,'CurrentPoint');  mp = get(hMenu,'Position');
        tf = cp(1) >= mp(1) && cp(1) <= mp(1)+mp(3) ...
          && cp(2) >= mp(2) && cp(2) <= mp(2)+mp(4);
    end

    % =====================================================================
    %  Editing
    % =====================================================================
    function addEl(type,strength,z)
        if nargin < 3
            z = (2*rand - 1.5) + 1i*(1.6*rand - 0.8);
        end
        ez(end+1) = clampDom(z);  et(end+1) = type;
        es(end+1) = strength;     ea(end+1) = 0;
        eb(end+1) = false;  ez2(end+1) = NaN;  eOn(end+1) = 0;
        sel = numel(ez);
        redraw();
    end

    function removeEl(k)
        if k < 1 || k > numel(ez), return; end
        ez(k) = [];  et(k) = [];  es(k) = [];  ea(k) = [];
        eb(k) = [];  ez2(k) = [];  eOn(k) = [];
        sel = 0;  dragIdx = 0;
        if menuIdx == k
            closeMenu();
        elseif menuIdx > k
            menuIdx = menuIdx - 1;
        end
        redraw();
    end

    function clearAll()
        ez = zeros(1,0);  et = ez;  es = ez;  ea = ez;
        eb = false(1,0);  ez2 = ez;  eOn = ez;
        sel = 0;  dragIdx = 0;
        closeMenu();
        leaveDye();
        redraw();
    end

    function setU()
        U = get(hUsl,'Value');
        set(hUtxt,'String',sprintf('U = %.2f',U));
        redraw();
    end

    function apiSetU(val)
        set(hUsl,'Value',min(max(val,0),2));
        setU();
    end

    function z = clampDom(z)
        z = min(max(real(z),XL(1)),XL(2)) + 1i*min(max(imag(z),YL(1)),YL(2));
    end

    % =====================================================================
    %  Mouse and keyboard
    % =====================================================================
    function z = mousePos()
        cp = get(ax,'CurrentPoint');
        z = cp(1,1) + 1i*cp(1,2);
    end

    function k = pick(z)
        k = 0;
        if isempty(ez), return; end
        [d,i] = min(abs(ez - z));
        if d <= PICK, k = i; end
    end

    function onDown(~,~)
        if overMenu(), return; end
        z = mousePos();
        k = pick(z);
        if k == 0
            if real(z) >= XL(1) && real(z) <= XL(2) && imag(z) >= YL(1) && imag(z) <= YL(2)
                sel = 0;  closeMenu();  drawMarkers();
            end
            return
        end
        if strcmp(get(fig,'SelectionType'),'alt')
            openMenu(k);
        else
            sel = k;  dragIdx = k;  dragOff = ez(k) - z;
            if menuIdx > 0, menuIdx = k; end    % open panel follows the selection
            drawMarkers();
        end
    end

    function onMove(~,~)
        % BUSY drops mouse events that arrive while a frame is being drawn
        if dragIdx == 0 || busy, return; end
        busy = true;
        try
            ez(dragIdx) = clampDom(mousePos() + dragOff);
            redraw();
            drawnow limitrate
        catch err
            busy = false;
            rethrow(err);
        end
        busy = false;
    end

    function onUp(~,~)
        if dragIdx == 0, return; end
        ez(dragIdx) = clampDom(mousePos() + dragOff);
        dragIdx = 0;
        redraw();
    end

    function onScroll(~,evt)
        if overMenu(), return; end
        k = pick(mousePos());
        if k == 0, k = sel; end
        if k == 0, return; end
        n = -evt.VerticalScrollCount;           % wheel up increases
        if et(k) == DBL && any(strcmp(get(fig,'CurrentModifier'),'shift'))
            ea(k) = ea(k) + n*pi/12;
        else
            es(k) = round((es(k) + n*STEP(et(k)))*100)/100;
        end
        sel = k;
        redraw();
    end

    function onKey(~,evt)
        % typing in a panel field must not add or delete elements
        obj = get(fig,'CurrentObject');
        if ~isempty(obj) && isprop(obj,'Style') && strcmp(get(obj,'Style'),'edit')
            return
        end
        z = mousePos();
        inside = real(z) >= XL(1) && real(z) <= XL(2) && imag(z) >= YL(1) && imag(z) <= YL(2);
        new = [];
        switch evt.Key
            case 's', new = {SRC, 2};
            case 'k', new = {SRC,-2};
            case 'v', new = {VOR, 3};
            case 'c', new = {VOR,-3};
            case 'd', new = {DBL, 1.5};
            case {'delete','backspace'}
                if sel > 0, removeEl(sel); end
            case 'escape'
                sel = 0;  drawMarkers();
        end
        if ~isempty(new)
            if inside, addEl(new{:},z); else, addEl(new{:}); end
        end
    end
end
