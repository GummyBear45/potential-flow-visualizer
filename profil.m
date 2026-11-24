function varargout = profil
%PROFIL Airfoils in 2D potential flow: Joukowski mapping and thin airfoil theory.
%   PROFIL opens a window with two modes.
%
%   Joukowski (exact)
%     Flow around a circle in the auxiliary plane zeta is the PRZEPLYWER
%     cylinder: freestream + doublet + vortex at the circle centre mu. The
%     map z = zeta + b^2/zeta turns the circle into an airfoil. Moving mu
%     left gives thickness, moving it up gives camber; the Kutta condition
%     fixes the circulation so the flow leaves the trailing edge smoothly.
%
%   Thin airfoil (distributions)
%     Sources (density sigma) and vortices (density gamma, clockwise
%     positive) are spread along the chord. Drag their control points: the
%     thickness follows from sigma, the camber line from gamma
%         y_t(x) = 1/(2U) * int_0^x sigma,
%         dy_c/dx = v_gamma(x)/U,  v_gamma = -(1/2pi) PV int gamma/(x - xi).
%     The last gamma point is locked at 0 (Kutta condition). The theory is
%     linear, so the sliders scale sigma (thickness) and gamma (camber)
%     directly; the angle of attack adds the exact flat-plate sheet
%     gamma = 2 U alpha sqrt((1-x)/x), whose downwash is a uniform -U alpha.
%     "Keep body closed" shifts the other source points while one is
%     dragged, so that int sigma = 0. The small bottom slider adds an exact
%     round nose y_N = sqrt(2 r x)(1 - x) of radius r.
%
%   Boundary layer (both modes, optional)
%     Thwaites -> Michel transition -> Head, Squire-Young drag. Thin mode
%     uses the Lighthill-Van Dyke leading-edge correction of the surface
%     speed, which needs the nose radius.
%
%   Pressure coefficient
%     The right-hand plot shows -Cp along the chord. In thin airfoil theory
%     Cp is linear in the disturbance, so it splits exactly into thickness
%     (same on both surfaces), camber and angle of attack (opposite signs
%     on the two surfaces). In Joukowski mode the exact Cp is drawn over
%     the thin-theory curves for the same geometry.
%
%   FIG = PROFIL returns the figure handle.

    % ---- constants ------------------------------------------------------
    XL = [-0.8 1.8];  YL = [-0.75 0.75];    % physical window, chord about [0, 1]
    DG = 0.1;                               % arrow grid spacing
    b = 0.25;  X0 = 0.5;                    % Joukowski constant; shift puts TE at x = 1
    U = 1;                                  % freestream speed
    H = 0.02;  NSTEP = 450;  NLEFT = 26;    % streamlines
    NPAN = 120;  NFOU = 16;                 % thin-airfoil panels, Glauert terms
    NPART = 500;  DT = 0.012;  MAXHOP = 0.03;
    XC = [0 0.01 0.03 0.07 0.13 0.22 0.35 0.5 0.65 0.8 0.92 1];   % editor control points
    CPLIM = [-3 1.2];
    PICK = 0.05;

    % panels with cosine spacing; midpoints are uniform in theta
    xe = (1 - cos(pi*(0:NPAN)/NPAN))/2;
    thm = pi*((1:NPAN) - 0.5)/NPAN;  xm = (1 - cos(thm))/2;
    dxp = diff(xe);  dth = pi/NPAN;
    % principal-value kernel: int over panel j of 1/(x_i - xi) dxi
    Lmat = log(abs((xm.' - xe(1:end-1))./(xm.' - xe(2:end))));

    % viscous-inviscid interaction (thin mode): wake panels behind the TE
    NW = 40;  LW = 1.5;                     % wake panels and length [chords]
    VIT = 150;  OMG = 0.5;  VTOL = 5e-3;    % iterations, relaxation, tolerance
    HSW = 2.0;                              % shape factor where the march turns inverse
    NCRIT = 9;                              % e^N transition: critical amplification (XFOIL default)
    XCPL = 0.96;                            % coupling stops here (TE singularity of linear theory)
    xwE = 1 + LW*linspace(0,1,NW+1).^1.6;  xwM = (xwE(1:end-1) + xwE(2:end))/2;
    dxw = diff(xwE);
    kern = @(xq,ed) log(abs((xq.' - ed(1:end-1))./(xq.' - ed(2:end))));
    Law = kern(xm,xwE);  Lwa = kern(xwM,xe);  Lww = kern(xwM,xwE);

    % ---- state ----------------------------------------------------------
    mode = 1;                               % 1 Joukowski, 2 thin airfoil
    thk = 0.12;  cmb = 0.06;  aoa = deg2rad(5);
    tauTE = 0;  nKT = 2;                    % trailing-edge angle [deg], map exponent
    blOn = false;  ReN = 1e6;  bl = [];     % boundary layer: switch, Reynolds number, result
    mu = 0;  a = b;  betaK = 0;  Gam = 0;  Uc = U;  Ap = 0;  Bp = 0;
    sc = zeros(size(XC));  gc = zeros(size(XC));   % editor shapes (before scaling)
    kS = 1;  kG = 1;  aoaT = 0;             % thin mode: thickness/camber scale, AoA
    rN = 0.0159;                            % thin mode: nose radius at thickness scale 1
    sig = zeros(1,NPAN);  gam = sig;  gUsr = sig;   % densities; gUsr = editor part of gamma
    utC = sig;                              % thickness-induced surface speed u_t at xm
    sigUsr = sig;                           % editor part of sigma (without the nose term)
    sigV = sig;  gamV = sig;  sigW = zeros(1,NW);   % viscous transpiration / decambering
    vs = [];                                % VII state, kept as warm start
    ytE = zeros(1,NPAN+1);  ycE = ytE;  ytUsrE = ytE;
    info = {};
    pz = zeros(1,0);  pAge = pz;  pLife = pz;
    dragPt = [0 0];  busy = false;  quick = false;   % quick: fast BL solve while dragging

    % ---- figure ---------------------------------------------------------
    fig = figure('Name','Profil - airfoils in potential flow','NumberTitle','off', ...
        'Color','w','MenuBar','none','ToolBar','none', ...
        'Units','pixels','Position',[40 40 1420 860]);
    np = {'HitTest','off','PickableParts','none'};

    axF = axes('Parent',fig,'Units','normalized','Position',[0.04 0.40 0.56 0.56], ...
        'XLim',XL,'YLim',YL,'DataAspectRatio',[1 1 1],'Box','on','Layer','top');
    hold(axF,'on');  title(axF,'Flow around the airfoil');
    colormap(axF,turbo(256));  set(axF,'CLim',[0 2*U]);
    [gx,gy] = meshgrid(XL(1):DG:XL(2), YL(1):DG:YL(2));
    gz = gx + 1i*gy;
    hQ = quiver(axF,gx,gy,0*gx,0*gy,0,'Color',[0.45 0.49 0.58],'MaxHeadSize',0.6,np{:});
    hS = line('Parent',axF,'XData',NaN,'YData',NaN,'Color',[0.13 0.33 0.62],'LineWidth',1,np{:});
    hDiv = line('Parent',axF,'XData',NaN,'YData',NaN,'Color',[0.90 0.45 0.05],'LineWidth',2,np{:});
    hBody = patch('Parent',axF,'XData',NaN,'YData',NaN,'FaceColor',[0.80 0.80 0.80], ...
        'EdgeColor','k','LineWidth',1.2,np{:});
    hCamber = line('Parent',axF,'XData',NaN,'YData',NaN,'Color',[0.2 0.55 0.3], ...
        'LineStyle','--','LineWidth',1,np{:});
    hP = scatter(axF,NaN,NaN,12,0,'filled','MarkerEdgeColor',[0.15 0.15 0.15], ...
        'LineWidth',0.25,np{:});
    hStag = line('Parent',axF,'XData',NaN,'YData',NaN,'LineStyle','none','Marker','x', ...
        'MarkerSize',11,'LineWidth',2.2,'Color','k',np{:});

    % pressure coefficient along the chord (-Cp up)
    axC = axes('Parent',fig,'Units','normalized','Position',[0.66 0.56 0.31 0.38], ...
        'XLim',[0 1],'YLim',CPLIM,'YDir','reverse','Box','on');
    hold(axC,'on');  grid(axC,'on');
    title(axC,'Pressure coefficient along the chord');
    xlabel(axC,'x / c');  ylabel(axC,'C_p  (negative up)');
    cU = [0.75 0.15 0.12];  cL = [0.15 0.35 0.75];
    cT = [0.45 0.45 0.45];  cC = [0.20 0.60 0.30];  cA = [0.55 0.30 0.75];
    ln = @(c,s,w) line('Parent',axC,'XData',NaN,'YData',NaN,'Color',c,'LineStyle',s, ...
        'LineWidth',w,np{:});
    hExU = ln(cU,'-',2.2);   hExL = ln(cL,'-',2.2);
    hThU = ln(cU,'--',1.4);  hThL = ln(cL,'--',1.4);
    hCt  = ln(cT,'-',1.2);
    hCcU = ln(cC,'-',1.2);   hCcL = ln(cC,':',1.4);
    hCaU = ln(cA,'-',1.2);   hCaL = ln(cA,':',1.4);
    line('Parent',axC,'XData',[0 1],'YData',[0 0],'Color','k','LineWidth',0.5,np{:});
    hLeg = legend(axC,[hExU hExL hThU hThL hCt hCcU hCaU], ...
        {'exact, upper','exact, lower','thin theory, upper','thin theory, lower', ...
         'thickness part','camber part (dotted: lower)','angle-of-attack part (dotted: lower)'}, ...
        'Location','northeast','FontSize',8,'AutoUpdate','off');

    % bottom right: circle plane (Joukowski) or distribution editor (thin)
    pnlZ = uipanel(fig,'Units','normalized','Position',[0.62 0.03 0.37 0.45], ...
        'BorderType','none','BackgroundColor','w');
    axZ = axes('Parent',pnlZ,'Units','normalized','Position',[0.10 0.12 0.86 0.78], ...
        'XLim',[-0.85 0.85],'YLim',[-0.55 0.55],'DataAspectRatio',[1 1 1],'Box','on');
    hold(axZ,'on');
    title(axZ,['Circle plane ' char(950) ': freestream + doublet + vortex at ' char(956)]);
    hZs = line('Parent',axZ,'XData',NaN,'YData',NaN,'Color',[0.13 0.33 0.62],'LineWidth',0.8,np{:});
    hZd = line('Parent',axZ,'XData',NaN,'YData',NaN,'Color',[0.90 0.45 0.05],'LineWidth',1.8,np{:});
    hZc = patch('Parent',axZ,'XData',NaN,'YData',NaN,'FaceColor',[0.85 0.85 0.85], ...
        'EdgeColor','k',np{:});
    line('Parent',axZ,'XData',[-b b],'YData',[0 0],'LineStyle','none','Marker','+', ...
        'MarkerSize',10,'LineWidth',1.5,'Color',[0.8 0.1 0.1],np{:});
    hZmu = line('Parent',axZ,'XData',NaN,'YData',NaN,'LineStyle','none','Marker','d', ...
        'MarkerSize',10,'MarkerFaceColor',[0.25 0.25 0.25],'MarkerEdgeColor','k',np{:});
    hZstag = line('Parent',axZ,'XData',NaN,'YData',NaN,'LineStyle','none','Marker','x', ...
        'MarkerSize',10,'LineWidth',2,'Color','k',np{:});
    hZtxt = text(axZ,-0.82,-0.47,'','FontSize',9,'BackgroundColor',[1 1 1 0.8],np{:});

    pnlE = uipanel(fig,'Units','normalized','Position',[0.62 0.03 0.37 0.45], ...
        'BorderType','none','BackgroundColor','w','Visible','off');
    axE = axes('Parent',pnlE,'Units','normalized','Position',[0.10 0.12 0.86 0.78], ...
        'XLim',[0 1],'YLim',[-1 3],'Box','on');
    hold(axE,'on');  grid(axE,'on');
    title(axE,['Drag the points: sources ' char(963) '/U (red), vortices ' char(947) '/U (green)']);
    xlabel(axE,'x / c');
    line('Parent',axE,'XData',[0 1],'YData',[0 0],'Color','k','LineWidth',0.5,np{:});
    xf = linspace(0,1,300);
    hEs = line('Parent',axE,'XData',xf,'YData',0*xf,'Color',[0.84 0.22 0.18],'LineWidth',1.6,np{:});
    hEg = line('Parent',axE,'XData',xf,'YData',0*xf,'Color',[0.2 0.6 0.3],'LineWidth',1.6,np{:});
    hEa = line('Parent',axE,'XData',xf,'YData',NaN*xf,'Color',[0.55 0.30 0.75], ...
        'LineStyle','--','LineWidth',1.3,np{:});
    text(axE,0.45,2.45,'dashed: AoA part of \gamma (exact)','FontSize',9, ...
        'Color',[0.55 0.30 0.75],np{:});
    hEn = line('Parent',axE,'XData',xf,'YData',NaN*xf,'Color',[0.84 0.22 0.18], ...
        'LineStyle','--','LineWidth',1.3,np{:});
    text(axE,0.45,2.15,'dashed: round-nose part of \sigma (exact)','FontSize',9, ...
        'Color',[0.84 0.22 0.18],np{:});
    hEsp = line('Parent',axE,'XData',XC,'YData',sc,'LineStyle','none','Marker','o','MarkerSize',8, ...
        'MarkerFaceColor',[0.84 0.22 0.18],'MarkerEdgeColor','k',np{:});
    hEgp = line('Parent',axE,'XData',XC(1:end-1),'YData',gc(1:end-1),'LineStyle','none', ...
        'Marker','s','MarkerSize',8,'MarkerFaceColor',[0.2 0.6 0.3],'MarkerEdgeColor','k',np{:});
    line('Parent',axE,'XData',1,'YData',0,'LineStyle','none','Marker','s','MarkerSize',8, ...
        'MarkerEdgeColor',[0.2 0.6 0.3],'LineWidth',1.5,np{:});
    text(axE,0.70,2.75,[char(947) '(1) = 0: Kutta condition'],'FontSize',9,np{:});

    % boundary layer: displacement thickness (exaggerated), transition, separation
    BLX = 5;                                % drawing exaggeration of delta*
    hBL = line('Parent',axF,'XData',NaN,'YData',NaN,'Color',[0.85 0.33 0.10], ...
        'LineWidth',1.4,np{:});
    hTr = line('Parent',axF,'XData',NaN,'YData',NaN,'LineStyle','none','Marker','^', ...
        'MarkerSize',8,'MarkerFaceColor',[0.2 0.65 0.3],'MarkerEdgeColor','k',np{:});
    hSep = line('Parent',axF,'XData',NaN,'YData',NaN,'LineStyle','none','Marker','p', ...
        'MarkerSize',13,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k',np{:});

    pf.prepareFigure(fig,[axF axC axZ axE]);

    % ---- controls -------------------------------------------------------
    uicontrol(fig,'Style','text','String','Mode','Units','normalized', ...
        'Position',[0.04 0.297 0.04 0.03],'BackgroundColor','w','FontSize',10, ...
        'HorizontalAlignment','left');
    hMode = uicontrol(fig,'Style','popupmenu','Units','normalized', ...
        'Position',[0.08 0.30 0.20 0.035],'FontSize',10, ...
        'String',{'Joukowski airfoil (exact)','Thin airfoil (source/vortex distributions)'}, ...
        'Callback',@(src,~) setMode(get(src,'Value')));
    chk = @(s,pos) uicontrol(fig,'Style','checkbox','String',s,'Value',1,'Units','normalized', ...
        'Position',pos,'BackgroundColor','w','FontSize',10,'Callback',@(~,~) redraw());
    cVec = chk('Vectors',[0.30 0.30 0.07 0.035]);
    cStr = chk('Streamlines',[0.37 0.30 0.09 0.035]);
    cPart = chk('Particles',[0.46 0.30 0.08 0.035]);
    cComp = chk('Cp components',[0.04 0.06 0.12 0.035]);
    cClose = uicontrol(fig,'Style','checkbox','String','Keep body closed','Value',1, ...
        'Units','normalized','Position',[0.04 0.02 0.12 0.035],'BackgroundColor','w', ...
        'FontSize',10,'Callback',@(~,~) toggleClosed());

    lbl = @(s,y) uicontrol(fig,'Style','text','String',s,'Units','normalized', ...
        'Position',[0.04 y 0.13 0.03],'BackgroundColor','w','FontSize',10, ...
        'HorizontalAlignment','left');
    sld = @(lo,hi,v,y) uicontrol(fig,'Style','slider','Min',lo,'Max',hi,'Value',v, ...
        'Units','normalized','Position',[0.17 y+0.004 0.17 0.028]);
    val = @(y) uicontrol(fig,'Style','text','String','','Units','normalized', ...
        'Position',[0.345 y 0.06 0.03],'BackgroundColor','w','FontSize',10, ...
        'HorizontalAlignment','left');
    % the same three sliders serve both modes; setMode relabels them
    hL1 = lbl('',0.245);  hThk = sld(0,0.3,thk,0.245);  hThkV = val(0.245);
    hL2 = lbl('',0.200);  hCmb = sld(-0.15,0.3,cmb,0.200); hCmbV = val(0.200);
    hL3 = lbl('',0.155);  hAoa = sld(-10,15,rad2deg(aoa),0.155); hAoaV = val(0.155);
    % while a slider is dragged a quick (8-iteration) boundary-layer solve
    % keeps up; the full solve runs when it is released
    for hs_ = [hThk hCmb hAoa]
        set(hs_,'Callback',@(~,~) readSliders(false));
        addlistener(hs_,'ContinuousValueChange',@(~,~) readSliders(true));
    end

    btn = @(s,x,w,cb) uicontrol(fig,'Style','pushbutton','String',s,'Units','normalized', ...
        'Position',[x 0.105 w 0.04],'FontSize',10,'Callback',cb);
    btn('Symmetric',0.04,0.075,@(~,~) preset(1));
    btn('Cambered',0.12,0.075,@(~,~) preset(2));
    btn('Plate, 5 deg',0.20,0.075,@(~,~) preset(3));
    btn('Polar...',0.28,0.075,@(~,~) askPolar());

    % boundary layer (both modes); the slider below is the trailing-edge
    % angle in Joukowski mode and the nose radius in thin-airfoil mode
    cBL = uicontrol(fig,'Style','checkbox','String','Boundary layer','Value',0, ...
        'Units','normalized','Position',[0.17 0.06 0.11 0.035],'BackgroundColor','w', ...
        'FontSize',10,'Callback',@(~,~) setBL());
    uicontrol(fig,'Style','text','String','Re','Units','normalized', ...
        'Position',[0.17 0.018 0.025 0.03],'BackgroundColor','w','FontSize',10, ...
        'HorizontalAlignment','left');
    hRe = uicontrol(fig,'Style','edit','String','1e6','Units','normalized', ...
        'Position',[0.195 0.022 0.075 0.033],'BackgroundColor','w','FontSize',10, ...
        'Callback',@(~,~) setBL());
    hTEl = uicontrol(fig,'Style','text','String','','Units','normalized', ...
        'Position',[0.29 0.058 0.11 0.03],'BackgroundColor','w','FontSize',10, ...
        'HorizontalAlignment','left');
    hTE = uicontrol(fig,'Style','slider','Min',0,'Max',20,'Value',tauTE, ...
        'Units','normalized','Position',[0.29 0.026 0.11 0.028],'Callback',@(~,~) setBottom(false));
    addlistener(hTE,'ContinuousValueChange',@(~,~) setBottom(true));

    hInfo = uicontrol(fig,'Style','text','Units','normalized','Position',[0.41 0.10 0.19 0.19], ...
        'BackgroundColor',[0.97 0.97 0.97],'FontSize',9,'HorizontalAlignment','left', ...
        'FontName','Consolas');
    hHint = uicontrol(fig,'Style','text','Units','normalized','Position',[0.41 0.005 0.19 0.09], ...
        'BackgroundColor','w','ForegroundColor',[0.35 0.35 0.35],'FontSize',8, ...
        'HorizontalAlignment','left');

    set(fig,'WindowButtonDownFcn',@onDown,'WindowButtonMotionFcn',@onMove, ...
        'WindowButtonUpFcn',@onUp);

    setappdata(fig,'api',struct('setMode',@setMode,'setJoukowski',@apiSetJoukowski, ...
        'setThin',@apiSetThin,'setBL',@apiSetBL,'setControls',@apiSetControls,'preset',@preset, ...
        'polar',@makePolar, ...
        'velocity',@physVel,'state',@getState,'cp',@getCp,'tick',@tick));

    preset(2);
    setMode(1);
    respawn(true(1,NPART));
    pAge = floor(rand(1,NPART).*pLife);
    tmr = timer('ExecutionMode','fixedSpacing','Period',0.03,'BusyMode','drop', ...
        'TimerFcn',@tick);
    % the handler is in its own file, so closing works even after this file is edited
    set(fig,'DeleteFcn',{@pf.deleteTimer,tmr});
    start(tmr);
    if nargout, varargout{1} = fig; end

    % =====================================================================
    %  Joukowski airfoil
    % =====================================================================
    function setJoukowski()
        % circle through the trailing-edge point zeta = b, Kutta circulation.
        % The map exponent n sets the trailing-edge angle (Karman-Trefftz);
        % n = 2 is the cusped Joukowski airfoil. X0 keeps the TE at x = 1.
        nKT = 2 - deg2rad(tauTE)/pi;
        X0 = 1 - nKT*b;
        mu = b*(-thk + 1i*cmb);
        a = abs(b - mu);
        betaK = -angle(b - mu);
        Gam = -4*pi*a*U*sin(aoa + betaK);       % counter-clockwise positive
        Uc = U*exp(-1i*aoa);
        Ap = -1i*Gam/(2*pi);                     % vortex at mu
        Bp = -U*exp(1i*aoa)*a^2;                 % doublet making the circle a streamline
    end

    function w = wZeta(zeta)
        w = pf.poleVelocity(zeta,Uc,mu,Ap,Bp);
    end

    function z = fmap(zeta)
        t = ((zeta - b)./(zeta + b)).^nKT;
        z = nKT*b*(1 + t)./(1 - t) + X0;
    end

    function d = dmap(zeta)
        t = ((zeta - b)./(zeta + b)).^nKT;
        d = 4*nKT^2*b^2*t./((1 - t).^2.*(zeta.^2 - b^2));
    end

    function [zeta,inside] = invmap(z)
        % the preimage outside the circle; points inside the airfoil have none
        w = (z - X0 - nKT*b)./(z - X0 + nKT*b);
        best = nan(size(z));  dBest = -inf(size(z));
        for k = -1:1
            ph = (angle(w) + 2*pi*k)/nKT;        % t must have its angle in (-pi, pi]
            t = abs(w).^(1/nKT).*exp(1i*ph);
            zk = b*(1 + t)./(1 - t);
            ok = ph > -pi & ph <= pi & abs(zk - mu) > dBest;
            best(ok) = zk(ok);  dBest(ok) = abs(zk(ok) - mu);
        end
        zeta = best;
        inside = dBest < a;
    end

    function w = wParticleJ(zeta)
        % particles are moved in the circle plane: d zeta/dt = conj(w_zeta)/|f'|^2
        w = wZeta(zeta)./abs(dmap(zeta)).^2;
    end

    function [xu,cpu,xl,cpl,geo] = exactCp()
        % exact surface Cp, split at the leading edge (the point farthest from TE)
        s = linspace(1e-3,2*pi - 1e-3,6000);
        zeta = mu + a*exp(1i*(-betaK + s));     % starts at the TE, upper surface first
        z = fmap(zeta);
        V = abs(wZeta(zeta)./dmap(zeta));
        cp = 1 - (V/U).^2;
        zTE = fmap(b);
        [~,iLE] = max(abs(z - zTE));
        zLE = z(iLE);
        c = abs(zTE - zLE);  eC = (zTE - zLE)/c;
        loc = (z - zLE).*conj(eC)/c;             % chord frame, LE at 0, TE at 1
        xu = real(loc(1:iLE));    cpu = cp(1:iLE);
        xl = real(loc(iLE:end));  cpl = cp(iLE:end);
        % thickness and camber relative to the chord line, on the panel edges
        % (pchip, not linear: the camber slope is differentiated afterwards)
        [xsu,iu] = unique(xu);  etaU = imag(loc(1:iLE));  etaU = etaU(iu);
        [xsl,il] = unique(xl);  etaL = imag(loc(iLE:end));  etaL = etaL(il);
        yu = pchip(xsu,etaU,xe);
        yl = pchip(xsl,etaL,xe);
        yu([1 end]) = 0;  yl([1 end]) = 0;
        geo = struct('yt',(yu - yl)/2,'yc',(yu + yl)/2,'c',c, ...
            'alphaC',aoa - angle(eC),'zLE',zLE,'zTE',zTE, ...
            'zeta',zeta,'z',z,'V',V,'x',real(loc));
    end

    % =====================================================================
    %  Boundary layer (thin-airfoil mode): Thwaites (laminar), Michel
    %  (transition), Head (turbulent), Squire-Young (profile drag)
    % =====================================================================
    function [wU,wL] = wakeEdges(zEndU,zEndL,dU0,dL0,dsW,zc,zTE)
        % Outline of the wake's displacement thickness for the drawing.
        % zc is the wake centre line at the stations xwM (arc length
        % xwM - 1 from the trailing edge zTE). Each edge starts where its
        % surface line ends (the march stops at x/c = 0.99), the thickness
        % is split in the ratio of the two layers at the trailing edge and
        % relaxes to an even split as the wake becomes symmetric, and the
        % edge is blended out of the surface end point over 0.05 chord, so a
        % thick (separated) layer meets the thinner wake without a step.
        % Drawing only: no number depends on it.
        sw = xwM - 1;
        Tw = gradient(zc);  nw = 1i*Tw./abs(Tw);  t0 = Tw(1)/abs(Tw(1));
        dw = movmean(dsW,5);                     % first wake panels are very short
        ru = 0.5 + (dU0/max(dU0 + dL0,eps) - 0.5)*exp(-sw/0.15);
        wU = oneEdge(zEndU,BLX*dw.*ru);
        wL = oneEdge(zEndL,-BLX*dw.*(1 - ru));

        function zb = oneEdge(z0,off)
            loc = (z0 - zTE)*conj(t0);  s0 = real(loc);  o0 = imag(loc);
            behind = sw > s0;                    % no hook back past the end point
            if ~any(behind), zb = z0;  return, end
            iFirst = find(behind,1);
            o = off(behind) + (o0 - off(iFirst))*exp(-(sw(behind) - s0)/0.08);
            zb = [z0, zc(behind) + o.*nw(behind)];
        end
    end

    % =====================================================================
    %  Viscous-inviscid interaction, thin-airfoil mode
    %  The boundary layer displaces the outer flow by delta*. In thin
    %  theory that is (i) extra sources d(Ue delta*)/ds on the chord line
    %  and in the wake, which change the surface speed, and (ii) an extra
    %  camber (delta*_u - delta*_l)/2, re-solved with Glauert's series,
    %  which lowers the lift. The turbulent layer is marched in inverse
    %  mode (delta* given, Ue unknown), which stays regular through
    %  separation, and delta* is corrected with Carter's semi-inverse rule
    %      delta*_new = delta* (1 + w (Ue_BL/Ue_outer - 1))
    %  until both sides agree. Drag comes from Squire-Young at the far end
    %  of the wake, so it stays valid with separation on the airfoil.
    % =====================================================================
    function out = viscousThin(maxIt)
        % Warm start from the last solution. If that goes badly (a big jump
        % in the settings), step the angle of attack up from the last good
        % angle (or 0) in 2 degree steps: coupled solutions are found far
        % more reliably by continuation than from scratch.
        good = @(o) ~o.fail && (o.conv || o.res < 0.02);
        out = viscousSolve(maxIt);
        if good(out) || maxIt <= 8, return, end
        keep = vs;  k1 = out;  aT = aoaT;
        a0 = 0;
        if ~isempty(keep) && isfield(keep,'aoa') && isfinite(keep.aoa) && keep.aoa < aT
            a0 = keep.aoa;
        else
            vs = [];
        end
        n = max(1,ceil(abs(aT - a0)/deg2rad(2)));
        for kk = 1:n
            aoaT = a0 + (aT - a0)*kk/n;
            o = viscousSolve(maxIt);
            if o.fail, break, end
        end
        aoaT = aT;
        out = viscousSolve(maxIt);
        if out.fail || (~k1.fail && k1.res < out.res)
            out = k1;  vs = keep;
        end
    end

    function out = viscousSolve(maxIt)
        rA = kS^2*rN;  c = hypot(1,ycE(end));  nu = U*c/ReN;
        % Lighthill-Van Dyke leading-edge correction factors
        Fm = sqrt(xm./(xm + rA/2));  Gm = sqrt((1 - xm)./(xm + rA/2));
        % second order: the lifting part of the speed is carried by the
        % local thickness flow, q = U(1 + u_t/U)(1 + q_lift/U) to first order
        % in each; exact for an ellipse at any angle, and it gives the
        % thickness gain in lift slope (about 1 + 0.77 t/c)
        fac = min(max(1 + utC/U,0.5),2);
        quI = U*((1 + utC/U).*Fm + fac.*(gUsr/(2*U).*Fm + aoaT*Gm));
        qlI = U*((1 + utC/U).*Fm - fac.*(gUsr/(2*U).*Fm + aoaT*Gm));
        cpInt = @(qa,qb) sum((qa.^2 - qb.^2).*dxp)/(U^2*hypot(1,ycE(end)));
        lead = (1 + cos(thm))./sin(thm);  nn = (1:NFOU).';
        utW0 = (Lwa*sigUsr.').'/(2*pi) + kS*U*sqrt(2*rN)/pi*noseU(xwM);
        if isempty(vs) || ~all(isfinite([vs.dsU vs.dsL vs.dsW vs.UeU vs.UeL vs.UeW]))
            vs = struct('dsU',zeros(1,NPAN),'dsL',zeros(1,NPAN),'dsW',zeros(1,NW), ...
                'UeU',U*ones(1,NPAN),'UeL',U*ones(1,NPAN),'UeW',U*ones(1,NW),'cold',true, ...
                'turb',false(2,NPAN),'aoa',NaN);
        end
        prev = vs;  xTrFix = [NaN NaN];  xLSFix = [NaN NaN];  xSwFix = [inf inf];  jFix = NaN;
        om = OMG;  best = inf;  retries = 0;     % relaxation, cut only when diverging
        res = inf;  it = 0;  rW = 0;  hist = nan(1,maxIt);
        while it < maxIt
            it = it + 1;
            % --- outer flow with the displacement effect
            [sa,sw,dyc] = transpiration();
            dut  = (Lmat*sa.' + Law*sw.').'/(2*pi);
            dutW = (Lwa*sa.' + Lww*sw.').'/(2*pi);
            dS = diff(dyc)./dxp;
            dA0 = -(1/pi)*sum(dS)*dth;
            dAn = (2/pi)*sum(dS.*cos(nn*thm),2)*dth;
            dgc = 2*U*sum(dAn.*sin(nn*thm),1);
            qu = U*((1 + (utC + dut)/U).*Fm + fac.*((gUsr + dgc)/(2*U).*Fm + (aoaT + dA0)*Gm));
            ql = U*((1 + (utC + dut)/U).*Fm - fac.*((gUsr + dgc)/(2*U).*Fm + (aoaT + dA0)*Gm));
            UeWi = U + utW0 + dutW;

            % --- boundary layer from the stagnation point along both surfaces
            [P,zs,jS] = thinPaths(qu,ql,jFix);
            if it == 5, jFix = jS; end          % freeze the stagnation panel too
            dsS = [vs.dsU; vs.dsL];  UeS = [vs.UeU; vs.UeL];
            dsN = dsS;  r = cell(1,2);  res = 0;  turbN = false(2,NPAN);
            for k = 1:2
                p = P{k};  on = p.surf > 0;
                s = [0 cumsum(abs(diff(p.z)))];
                % prescribed delta* for the inverse part: last iteration's,
                % only where the layer was turbulent then
                dsP = nan(size(p.x));
                if ~vs.cold
                    ii = sub2ind([2 NPAN],p.surf(on),p.id(on));
                    v = dsS(ii);  v(~vs.turb(ii)) = NaN;
                    dsP(on) = v;
                end
                r{k} = marchVII(s,max(p.Ue,1e-9),dsP,nu,p.x,[xTrFix(k) xLSFix(k)],xSwFix(k));
                rr = r{k};
                % new delta*: direct stations (laminar, attached turbulent)
                % under-relaxed towards the march, inverse stations by Carter
                ii = sub2ind([2 NPAN],p.surf(on),p.id(on));
                dnew = rr.dstar;
                old = nan(size(p.x));  old(on) = dsS(ii);
                % local relaxation (Le Balleur): the outer response to a
                % delta* wiggle of spacing ds grows like pi/ds, so the step
                % must shrink like ds/(pi delta*) where the panels are dense
                dsl = gradient(s);
                omLoc = min(om,dsl./(pi*max(rr.dstar,1e-6)));
                dirK = 1:min(rr.ksw - 1,numel(s));
                inC = p.x <= XCPL;                   % stations taking part in the coupling
                okD = dirK(isfinite(old(dirK)) & dirK > 3 & ~vs.cold & inC(dirK));
                if ~isempty(okD)
                    rD = max(abs(dnew(okD) - old(okD)))/max(max(dnew(okD)),1e-9);
                    res = max(res,rD);
                    w_ = omLoc(okD);
                    dnew(okD) = (1 - w_).*old(okD) + w_.*dnew(okD);
                end
                if rr.ksw <= numel(s)
                    t = rr.ksw:numel(s);  t = t(inC(t));
                end
                if rr.ksw <= numel(s) && ~isempty(t)
                    ratio = rr.UeB(t)./max(p.Ue(t),0.05*U);
                    res = max(res,max(abs(ratio - 1)));
                    % Carter's step, limited to -20 % / +25 % per iteration
                    dnew(t) = rr.dsUsed(t).*min(max(1 + omLoc(t).*(ratio - 1),0.8),1.25);
                end
                % near the trailing edge linear theory stops the flow (finite
                % TE angle); keep that out of the iteration and extend delta*
                % linearly from just upstream
                tc = find(inC & (1:numel(s)) > max(rr.itr,3),2,'last');
                if numel(tc) == 2 && any(~inC)
                    slope = max((dnew(tc(2)) - dnew(tc(1)))/max(s(tc(2)) - s(tc(1)),1e-9),0);
                    tail = find(~inC & (1:numel(s)) > tc(2));
                    dnew(tail) = dnew(tc(2)) + slope*(s(tail) - s(tc(2)));
                end
                dsN(ii) = max(dnew(on),0);
                UeS(ii) = rr.UeB(on);
                tb = false(size(p.x));  tb(rr.itr:end) = rr.itr < numel(p.x);
                turbN(ii) = tb(on);
                if it == 5                          % freeze transition and laminar separation
                    xTrFix(k) = p.x(rr.itr);  xLSFix(k) = rr.xLS;
                end
                % freeze the direct/inverse switch as well: a station flipping
                % between the two modes keeps the iteration from settling
                if it == 8 && rr.ksw <= numel(s), xSwFix(k) = p.x(rr.ksw); end
                % stations past x/c = 0.99 are not marched: hold the last value
                last = find(on,1,'last');
                rest = p.id(last)+1:NPAN;
                dsN(p.surf(last),rest) = dsN(p.surf(last),p.id(last));
                UeS(p.surf(last),rest) = UeS(p.surf(last),p.id(last));
            end
            vs.dsU = dsN(1,:);  vs.dsL = dsN(2,:);
            vs.UeU = UeS(1,:);  vs.UeL = UeS(2,:);  vs.turb = turbN;

            % --- wake: the two layers merge at the trailing edge
            eU = find(isfinite(r{1}.theta),1,'last');  eL = find(isfinite(r{2}.theta),1,'last');
            th0 = r{1}.theta(eU) + r{2}.theta(eL);
            ds0 = r{1}.dstar(eU) + r{2}.dstar(eL);
            Ue0 = (r{1}.UeB(eU) + r{2}.UeB(eL))/2;
            [thW,HW] = wakeMarch(th0,ds0,Ue0,UeWi);
            dsWn = HW.*thW;
            if vs.cold
                vs.dsW = dsWn;
            else
                rW = mean(abs(dsWn - vs.dsW))/max(mean(dsWn),1e-9);
                wW = min(0.3,dxw./(pi*max(dsWn,1e-6)));   % local, as on the airfoil
                vs.dsW = (1 - wW).*vs.dsW + wW.*dsWn;
            end
            vs.UeW = UeWi;
            resBL = res;                             % boundary layer only: drives the relaxation
            if ~vs.cold, res = max(res,rW); end
            wasCold = vs.cold;  vs.cold = false;
            if wasCold, res = inf; end               % nothing to compare with yet
            hist(it) = res;
            bad = ~all(isfinite([vs.dsU vs.dsL vs.dsW vs.UeU vs.UeL])) ...
                || isnan(res) || max([vs.dsU vs.dsL]) > 0.2*c ...
                || ~all(isfinite(thW)) || ~all(isfinite([qu ql]));
            if bad
                % diverged: return to the last good state with a smaller step
                retries = retries + 1;
                vs = prev;  om = om/2;  best = inf;
                if retries > 3, break, end
                continue
            end
            % the last good iteration is what gets reported
            G = struct('r',{r},'P',{P},'zs',zs,'qu',qu,'ql',ql,'thW',thW,'HW',HW, ...
                'UeWi',UeWi,'ds0',ds0,'sa',sa,'sw',sw,'dA0',dA0,'dAn',dAn,'dgc',dgc);
            if wasCold || resBL < best
                best = min(resBL,best);  prev = vs;  % remember the best state so far
            elseif resBL > 3*best
                om = max(om*0.7,0.05);               % growing: relax harder
            end
            if res < VTOL && it > 5, break; end
        end
        if ~exist('G','var')
            % never got a good iteration: report nothing rather than garbage
            vs = [];  sigV(:) = 0;  gamV(:) = 0;  sigW(:) = 0;
            out = struct('fail',true,'conv',false,'it',it,'res',inf);
            return
        end
        r = G.r;  P = G.P;  zs = G.zs;  qu = G.qu;  ql = G.ql;  thW = G.thW;  HW = G.HW;
        UeWi = G.UeWi;  sa = G.sa;  sw = G.sw;
        dA0 = G.dA0;  dgc = G.dgc;

        % viscous corrections kept for the flow field
        sigV = sa;  sigW = sw;  gamV = 2*U*dA0*lead + dgc;

        % results for drawing and the info box
        if res < 0.02, vs.aoa = aoaT; else, vs.aoa = NaN; end   % for continuation
        out = struct('fail',false,'cd',0,'sep',false,'it',it,'res',res,'conv',res < VTOL, ...
            'qu',qu,'ql',ql,'zStag',zs,'CL',cpInt(qu,ql),'CLi',cpInt(quI,qlI));
        for k = 1:2
            p = P{k};  rr = r{k};
            T = gradient(p.z);  nrm = (3 - 2*k)*1i*T./abs(T);
            out.dsurf{k} = p.z + BLX*rr.dstar.*nrm;
            if rr.itr < numel(p.z)
                out.zTr(k) = p.z(rr.itr);  out.xTr(k) = p.x(rr.itr);
            else
                out.zTr(k) = NaN;  out.xTr(k) = NaN;
            end
            out.zSep(k) = NaN;  out.xSep(k) = NaN;  out.xRe(k) = NaN;
            if rr.isep > 0
                out.zSep(k) = p.z(rr.isep);  out.xSep(k) = p.x(rr.isep);  out.sep = true;
                if rr.ireat > 0, out.xRe(k) = p.x(rr.ireat); end
            end
        end
        % wake drawn along the trailing-edge level, as in the outer flow
        yTE = ycE(end);
        eU = find(isfinite(r{1}.dstar),1,'last');  eL = find(isfinite(r{2}.dstar),1,'last');
        [out.wake,out.wakeL] = wakeEdges(out.dsurf{1}(eU),out.dsurf{2}(eL), ...
            r{1}.dstar(eU),r{2}.dstar(eL),vs.dsW,xwM + 1i*yTE,1 + 1i*yTE);
        out.cd = 2*thW(end)/c*(UeWi(end)/U)^((HW(end) + 5)/2);
        out.thW = thW;  out.HW = HW;  out.P = P;  out.r = r;  out.hist = hist(1:it);

        function [sa_,sw_,dyc_] = transpiration()
            % sources d(Ue delta*)/dx from the mass defect Ue*delta*, and the
            % displacement camber (delta*_u - delta*_l)/2, on the panel edges
            % (smoothed over a few panels: the drop of delta* at transition
            % would otherwise be a point sink that drags transition upstream)
            % (to first order, as everywhere in thin theory, the mass defect
            % is U*delta*; using the iterated Ue here feeds back into itself
            % and locks the iteration in a two-cycle)
            % Smoothing is over a fixed length (3 % of chord) rather than a
            % number of panels: the cosine panels are very dense at both
            % edges, where differentiating delta* would amplify any noise.
            mU = smoothX(U*vs.dsU);  mL = smoothX(U*vs.dsL);
            mW = U*vs.dsW;
            toE = @(v) [0 interp1(xm,v,xe(2:end-1)) v(end)];
            mUe = toE(mU);  mLe = toE(mL);
            mWe = [mUe(end) + mLe(end), interp1(xwM,mW,xwE(2:end-1)), mW(end)];
            sa_ = diff(mUe + mLe)./dxp;
            sw_ = diff(mWe)./dxw;
            dyc_ = (toE(smoothX(vs.dsU)) - toE(smoothX(vs.dsL)))/2;
        end

        function v = smoothX(v)
            % moving average over 0.03 chord on a uniform grid in x
            xu = linspace(xm(1),xm(end),301);
            vu = movmean(interp1(xm,v,xu),9);
            v = interp1(xu,vu,xm);
        end
    end

    function [P,zs,j] = thinPaths(qu,ql,jFix)
        % one path round the body (lower TE -> nose -> upper TE) with the
        % speed signed along it; the stagnation point is where it changes
        % sign. Returns the two downstream paths, each point tagged with its
        % surface (1 upper, 2 lower, 0 the stagnation point) and panel index.
        yc = interp1(xe,ycE,xm);
        yt = interp1(xe,ytUsrE,xm) + kS*sqrt(2*rN*xm).*(1 - xm);
        zu = xm + 1i*(yc + yt);  zl = xm + 1i*(yc - yt);
        Pz = [fliplr(zl) zu];  Px = [fliplr(xm) xm];  Q = [-fliplr(ql) qu];
        Sf = [2*ones(1,NPAN) ones(1,NPAN)];  Id = [NPAN:-1:1, 1:NPAN];
        j = find(Q(1:end-1) <= 0 & Q(2:end) > 0);
        if isfinite(jFix)
            j = jFix;                            % held fixed while iterating
        elseif isempty(j)
            j = NPAN;                            % no crossing: start at the nose
        else
            [~,m] = min(Px(j) + Px(j+1));  j = j(m);
        end
        w = Q(j)/(Q(j) - Q(j+1));
        if ~isfinite(w), w = 0.5; end
        zs = Pz(j) + w*(Pz(j+1) - Pz(j));  xs = Px(j) + w*(Px(j+1) - Px(j));
        ix = {j+1:2*NPAN, j:-1:1};  sg = [1 -1];
        P = cell(1,2);
        for k = 1:2
            q = ix{k};
            % stop at x/c = 0.99: linear theory slows the flow to zero at a
            % finite-angle trailing edge (log term in u_t)
            q = q(1:find(Px(q) <= 0.99 | (1:numel(q)) <= 3,1,'last'));
            P{k} = struct('z',[zs Pz(q)],'x',[xs Px(q)],'Ue',[0 sg(k)*Q(q)], ...
                'surf',[0 Sf(q)],'id',[0 Id(q)]);
        end
    end

    function r = marchVII(s,Ue,dsP,nu,xp,xFix,xSwFix)
        % Laminar: Thwaites (direct) up to e^N transition (N = NCRIT) or
        % laminar separation; with xTrFix given, transition is held at that x.
        % Turbulent: Head's momentum and entrainment equations, inverse
        % (delta* = dsP given, Ue solved) once dsP is known, direct on the
        % first pass.
        n = numel(s);
        theta = nan(1,n);  Hf = theta;  UeB = Ue;
        dU = gradient(Ue,s);
        dUs = gradient(movmean(Ue,9),s);         % pressure gradient seen by lambda
        m3 = min(3,n);  dUs(1:m3) = dU(1:m3);
        I = cumtrapz(s,Ue.^5);
        th2 = 0.45*nu*I./Ue.^6;
        th2(1) = 0.075*nu/max(dU(1),1e-9);
        xTrFix = xFix(1);  xLSFix = xFix(2);  frozen = isfinite(xTrFix);
        itr = n;  Namp = 0;  dNprev = 0;  sLS = NaN;  xLS = inf;   % sLS: laminar separation
        for k = 1:n
            theta(k) = sqrt(th2(k));
            lam = th2(k)/nu*dUs(k);
            Hf(k) = thwaitesHL(lam);
            ReT = Ue(k)*theta(k)/nu;
            % e^N: amplification of the most unstable Tollmien-Schlichting
            % waves, integrated along the surface (Drela-Giles envelope).
            % Past laminar separation (lambda < -0.09) the layer goes on as a
            % free shear layer whose H climbs quickly (here from 3.55 to 7
            % over about 0.06 chord), and the waves grow much faster.
            if isnan(sLS) && k > 3 && (frozen && xp(k) >= xLSFix || ~frozen && lam < -0.09)
                sLS = s(k);  xLS = xp(k);
            end
            Hamp = Hf(k);
            if isfinite(sLS)
                Hamp = min(3.55 + 60*(s(k) - sLS),7);
                % no laminar recovery inside the bubble: H stays at its
                % separation value (Thwaites' lambda past separation would
                % make delta* follow every wiggle of the local gradient)
                Hf(k) = thwaitesHL(-0.09);
            end
            dN = enRate(Hamp,ReT,theta(k));
            if k > 1, Namp = Namp + 0.5*(dNprev + dN)*(s(k) - s(k-1)); end
            dNprev = dN;
            if frozen
                % frozen: no laminar-separation override, or transition
                % flickers between two stations and the iteration never settles
                go = k > 3 && xp(k) >= xTrFix;
            else
                % Transition at N = Ncrit; after laminar separation this
                % happens inside a short bubble, as in XFOIL. The bubble is
                % not allowed to grow beyond 0.1 chord.
                go = k > 3 && (Namp >= NCRIT || (isfinite(sLS) && (s(k) - sLS > 0.1 || xp(k) >= XCPL - 0.02)));
            end
            if go, itr = k;  break, end
        end
        trWhy = 0;                                % 1: e^N, 2: bubble limit, +10: after laminar separation
        if itr < n, trWhy = 1 + (Namp < NCRIT) + 10*isfinite(sLS); end
        isep = 0;  ireat = 0;  ksw = n + 1;  dsUsed = nan(1,n);
        if itr < n
            % --- turbulent, direct (outer speed prescribed) while attached
            th = theta(itr);  Hf(itr) = 1.4;  Q = Ue(itr)*th*headH1(1.4);
            for k = itr:n-1
                h = s(k+1) - s(k);
                [d1,d2] = headRHS(th,Q,Ue(k),dU(k),nu);
                [e1,e2] = headRHS(th + h*d1,Q + h*d2,Ue(k+1),dU(k+1),nu);
                thN = th + h/2*(d1 + e1);  QN = Q + h/2*(d2 + e2);
                HN = headH(max(QN/(Ue(k+1)*thN),3.31));
                if HN > HSW || ~isfinite(HN) || (xp(k+1) >= xSwFix && isfinite(dsP(k+1)))
                    ksw = k + 1;  break
                end
                th = thN;  Q = QN;  theta(k+1) = th;  Hf(k+1) = HN;
            end
            % --- turbulent, inverse (delta* prescribed) once H passes HSW:
            % regular through separation, where the direct form is singular
            if ksw <= n
                k0 = ksw - 1;                     % last direct station
                ds0_ = Hf(k0)*theta(k0);
                kk = ksw:n;
                dsUsed(kk) = dsP(kk);
                if ~all(isfinite(dsUsed(kk)))
                    % no previous inverse solution here: let delta* keep
                    % growing at the rate it had when the switch happened
                    g = max((ds0_ - Hf(max(k0-3,itr))*theta(max(k0-3,itr))) ...
                        /max(s(k0) - s(max(k0-3,itr)),1e-9),0);
                    guess = ds0_ + 1.5*g*(s(kk) - s(k0));
                    bad_ = ~isfinite(dsUsed(kk));
                    dsUsed(kk(bad_)) = guess(bad_);
                end
                dsUsed(kk) = max(dsUsed(kk),1.1*ds0_/HSW);
                Ub = Ue(k0);  dsPrev = ds0_;
                Umax = 1.5*max(Ue);
                for k = k0:n-1
                    h = s(k+1) - s(k);  dd = (dsUsed(k+1) - dsPrev)/h;
                    [a1,b1] = invRHS(th,Ub,dsPrev,dd,nu);
                    [a2,b2] = invRHS(th + h*a1,Ub + h*b1,dsUsed(k+1),dd,nu);
                    th = max(th + h/2*(a1 + a2),1e-8);
                    Ub = min(max(Ub + h/2*(b1 + b2),0.02*U),Umax);
                    theta(k+1) = th;  UeB(k+1) = Ub;  Hf(k+1) = dsUsed(k+1)/th;
                    dsPrev = dsUsed(k+1);
                end
            end
            t = find(Hf(itr+2:n) > 2.4,1);
            if ~isempty(t)
                isep = itr + 1 + t;
                t2 = find(Hf(isep:n) < 2.0,1);
                if ~isempty(t2), ireat = isep + t2 - 1; end
            end
        end
        r = struct('xLS',xLS,'theta',theta,'H',Hf,'dstar',Hf.*theta,'UeB',UeB, ...
            'itr',itr,'isep',isep,'ireat',ireat,'ksw',ksw,'dsUsed',dsUsed, ...
            'trWhy',trWhy,'Ntr',Namp);
    end

    function dNdx = enRate(Hk,ReT,th)
        % Growth rate dN/ds of the envelope e^N method (Drela and Giles,
        % AIAA J. 25, 1987, as in XFOIL). Waves grow only above the
        % critical Re_theta0(H), switched on smoothly over +-0.08 in log10.
        Hk = max(Hk,1.05);
        hm = 1/(Hk - 1);
        lgRe0 = (1.415*hm - 0.489)*tanh(20*hm - 12.9) + 3.295*hm + 0.44;
        rn = (log10(max(ReT,1)) - (lgRe0 - 0.08))/0.16;
        if rn <= 0, dNdx = 0;  return, end
        ramp = 1;
        if rn < 1, ramp = 3*rn^2 - 2*rn^3; end
        dNdRe = 0.01*sqrt((2.4*Hk - 3.7 + 2.5*tanh(1.5*Hk - 4.65))^2 + 0.25);
        el = (6.54*Hk - 14.07)/Hk^2;
        em = (0.058*(Hk - 4)^2/(Hk - 1) - 0.068)/el;
        dNdx = ramp*dNdRe*0.5*(em + 1)*el/th;
    end

    function [dth,dU] = invRHS(th,Ub,ds,dds,nu)
        % inverse Head method: with delta* given, momentum and entrainment
        % equations form a 2x2 system for theta' and Ue' that stays regular
        % through separation (the direct form is singular there).
        % (The shape factor is Hs, not H: nested functions share variables
        % with the app, and H there is the streamline step.)
        Hs = min(max(ds/th,1.15),8);
        H1 = headH1(Hs);  H1H = headH1d(Hs);
        Fe = 0.0306*(H1 - 3)^(-0.6169);
        cf = ludwiegTillmann(Hs,Ub*th/nu);
        % [1, (2+H) th/Ue; Ue (H1 - H H1'), H1 th] [th'; Ue'] = [cf/2; Ue (F - H1' dds)]
        % solved explicitly; det = th ((2+H) H H1' - (1+H) H1) never vanishes
        % because H1' < 0
        a11 = 1;  a12 = (2 + Hs)*th/Ub;  a21 = Ub*(H1 - Hs*H1H);  a22 = H1*th;
        b1 = cf/2;  b2 = Ub*(Fe - H1H*dds);
        dt = a11*a22 - a12*a21;
        dth = (b1*a22 - a12*b2)/dt;
        dU = (a11*b2 - a21*b1)/dt;
    end

    function d = headH1d(H)
        % derivative of headH1
        if H <= 1.6
            d = -1.287*0.8234*(H - 1.1)^(-2.287);
        else
            d = -3.064*1.5501*(H - 0.6778)^(-4.064);
        end
    end

    function [thW,HW] = wakeMarch(th0,ds0,Ue0,Uw)
        % Head's method with zero wall friction along the wake centre line
        xw = [1 xwM];  Uq = [Ue0 Uw];  n = numel(xw);
        dUq = gradient(Uq,xw);
        th = th0;  H0 = min(max(ds0/th0,1.12),6);  Q = Ue0*th*headH1(H0);
        thW = zeros(1,n-1);  HW = thW;
        for k = 1:n-1
            h = xw(k+1) - xw(k);
            [a1,b1] = wakeRHS(th,Q,Uq(k),dUq(k));
            [a2,b2] = wakeRHS(th + h*a1,Q + h*b1,Uq(k+1),dUq(k+1));
            th = th + h/2*(a1 + a2);  Q = Q + h/2*(b1 + b2);
            thW(k) = th;  HW(k) = headH(max(Q/(Uq(k+1)*th),3.31));
        end
    end

    function [dth,dQ] = wakeRHS(th,Q,Ue,dU)
        HH1 = max(Q/(Ue*th),3.31);  Hf = headH(HH1);
        dth = -(2 + Hf)*th/Ue*dU;               % no wall friction in the wake
        dQ = Ue*0.0306*(HH1 - 3)^(-0.6169);
    end

    function u = noseU(x)
        % thin-theory surface speed of the nose term sqrt(x)(1 - x), per
        % unit sqrt(2r)/pi: 3 + (1-3x)/(2 sqrt x) ln|(1+sqrt x)/(1-sqrt x)|
        q = sqrt(x);
        u = 3 + (1 - 3*x)./(2*q).*log(abs((1 + q)./(1 - q)));
    end

    function [Hf,ell] = thwaitesHL(lam)
        % Thwaites' correlations in the fits of Cebeci and Bradshaw
        lam = min(max(lam,-0.09),0.1);
        if lam >= 0
            ell = 0.22 + 1.57*lam - 1.8*lam^2;
            Hf = 2.61 - 3.75*lam + 5.24*lam^2;
        else
            ell = 0.22 + 1.402*lam + 0.018*lam/(lam + 0.107);
            Hf = 2.088 + 0.0731/(lam + 0.14);
        end
    end

    function [dth,dQ] = headRHS(th,Q,Ue,dU,nu)
        % momentum integral and entrainment equations, Q = Ue*theta*HH1
        HH1 = max(Q/(Ue*th),3.31);
        Hf = headH(HH1);
        cf = ludwiegTillmann(Hf,Ue*th/nu);
        dth = cf/2 - (2 + Hf)*th/Ue*dU;
        dQ = Ue*0.0306*(HH1 - 3)^(-0.6169);
    end

    function HH1 = headH1(Hf)
        if Hf <= 1.6
            HH1 = 3.3 + 0.8234*(Hf - 1.1)^(-1.287);
        else
            HH1 = 3.3 + 1.5501*(Hf - 0.6778)^(-3.064);
        end
    end

    function Hf = headH(HH1)
        % inverse of headH1
        if HH1 >= 5.3
            Hf = 1.1 + 0.86*(HH1 - 3.3)^(-0.777);
        else
            Hf = 0.6778 + 1.1536*(HH1 - 3.3)^(-0.326);
        end
    end

    function cf = ludwiegTillmann(H,ReT)
        cf = 0.246*10^(-0.678*H)*max(ReT,10)^(-0.268);
    end

    % =====================================================================
    %  Thin airfoil theory
    % =====================================================================
    function th = thinFromGeometry(yt,yc,alphaC)
        % sources from the thickness slope, vortices from Glauert's solution
        sig = 2*U*diff(yt)./dxp;
        utC = (Lmat*sig.').'/(2*pi);
        s = diff(yc)./dxp;                       % camber slope on the panels
        n = (1:NFOU).';
        An = (2/pi)*sum(s.*cos(n*thm),2)*dth;
        A0c = -(1/pi)*sum(s)*dth;
        lead = (1 + cos(thm))./sin(thm);
        gA = 2*U*alphaC*lead;
        gC = 2*U*(A0c*lead + sum(An.*sin(n*thm),1));
        gam = gA + gC;
        th = thinCp(gA,gC);
        th.CL = pi*(2*(alphaC + A0c) + An(1));
    end

    function solveDistributions()
        % shape from user distributions (freestream along x)
        % thin theory is linear: scaling sigma scales the thickness, scaling
        % gamma scales the camber; the AoA adds the exact flat-plate sheet,
        % leading-edge singularity included
        % The round nose is a separate exact thickness term
        %   y_N = sqrt(2 r x) (1 - x),  sigma_N = U sqrt(2 r) (1 - 3x)/sqrt(x),
        % closed at the trailing edge by itself; the editor shapes the rest.
        lead = (1 + cos(thm))./sin(thm);
        sigUsr = kS*U*pchip(XC,sc,xm);
        sigN = kS*U*sqrt(2*rN)*(1 - 3*xm)./sqrt(xm);
        sig = sigUsr + sigN;
        gUsr = kG*U*pchip(XC,gc,xm);
        gam = gUsr + 2*U*aoaT*lead;
        ytUsrE = [0 cumsum(sigUsr.*dxp)]/(2*U);
        ytE = ytUsrE + kS*sqrt(2*rN*xe).*(1 - xe);
        % u_t: panels for the editor part; the nose term in closed form, since
        % constant-strength panels cannot hold its 1/sqrt(x) density:
        %   (1/pi) PV int y_N'/(x - xi) dxi
        %     = sqrt(2r)/pi * (3 + (1 - 3x)/(2 sqrt x) ln((1 + sqrt x)/(1 - sqrt x)))
        utC = (Lmat*sigUsr.').'/(2*pi) + kS*U*sqrt(2*rN)/pi*noseU(xm);
        % downwash on the chord line: panels for the editor part, and the
        % exact uniform -U*alpha of the flat-plate sheet (panels would smear
        % its leading-edge singularity)
        v = -(Lmat*gUsr.').'/(2*pi) - U*aoaT;
        ycE = [0 cumsum(v.*dxp)]/U;
    end

    function th = thinCp(gA,gC)
        th.x = xm;
        th.t = -2*utC/U;
        th.au = -gA/U;  th.al = gA/U;
        th.cu = -gC/U;  th.cl = gC/U;
        th.u = th.t + th.au + th.cu;
        th.l = th.t + th.al + th.cl;
    end

    function w = wThin(z)
        % Source and vortex panels of constant density, placed along the
        % mean line (and the wake along the trailing-edge level) for the
        % flow picture. The strengths come from linear theory, which puts
        % them on y = 0; drawing them there would make the streamlines bend
        % around an invisible flat sheet instead of the airfoil drawn.
        % A straight panel from z1 to z2 at angle phi gives
        %   w = (sigma + i gamma) e^(-i phi)/(2 pi) ln((z - z1)/(z - z2))
        % (plus, with the boundary layer on, its displacement sources and
        % the decambering vortices)
        sz = size(z);  z = z(:);
        zE = xe + 1i*ycE;
        ph = angle(diff(zE));
        s = ((sig + sigV + 1i*(gam + gamV)).*exp(-1i*ph)/(2*pi)).';
        w = U + log((z - zE(1:end-1))./(z - zE(2:end)))*s;
        if any(sigW)
            zW = xwE + 1i*ycE(end);
            w = w + log((z - zW(1:end-1))./(z - zW(2:end)))*(sigW/(2*pi)).';
        end
        w = reshape(w,sz);
    end

    function in = insideThin(z)
        x = real(z);  y = imag(z);
        in = false(size(z));
        m = x >= 0 & x <= 1;
        if ~any(m), return; end
        yc = interp1(xe,ycE,x(m));  yt = interp1(xe,ytE,x(m));
        in(m) = abs(y(m) - yc) < max(yt,0) + 0.004;
    end

    % =====================================================================
    %  Drawing
    % =====================================================================
    function w = physVel(z)
        if mode == 1
            [zeta,inside] = invmap(z);
            w = wZeta(zeta)./dmap(zeta);
            w(inside) = NaN;
        else
            w = wThin(z);
            w(insideThin(z)) = NaN;
        end
    end

    function redraw()
        zLeft = XL(1) + 1i*(YL(1) + diff(YL)*((1:NLEFT) - 0.5)/NLEFT);
        boxOK = @(z) real(z) >= XL(1)-1e-9 & real(z) <= XL(2) ...
                   & imag(z) >= YL(1) & imag(z) <= YL(2);
        if mode == 1
            setJoukowski();
            bl = [];
            set([hBL hTr hSep],'XData',NaN,'YData',NaN);
            % body and circle
            zc = mu + a*exp(1i*linspace(0,2*pi,400));
            zb = fmap(zc);
            set(hBody,'XData',real(zb),'YData',imag(zb));
            set(hZc,'XData',real(zc),'YData',imag(zc));
            set(hZmu,'XData',real(mu),'YData',imag(mu));
            % streamlines are traced in the circle plane and mapped
            if get(cStr,'Value')
                keepS = @(zeta) abs(zeta - mu) > 1.002*a & boxOK(fmap(zeta));
                [xz,yz] = pf.traceStreamlines(@wZeta,invmap(zLeft),ones(1,NLEFT),H,NSTEP,keepS);
                zz = fmap(xz + 1i*yz);
                set(hS,'XData',real(zz),'YData',imag(zz),'Visible','on');
                set(hZs,'XData',xz,'YData',yz,'Visible','on');
            else
                set([hS hZs],'Visible','off');
            end
            % stagnation points: one at the front, one at the trailing edge
            zs = pf.stagnationPoints(Uc,mu,Ap,Bp,true);
            [~,iTE] = min(abs(zs - b));
            zF = zs(3 - iTE);  zT = zs(iTE);
            keep = @(zeta) abs(zeta - mu) > 0.999*a & boxOK(fmap(zeta));
            nrm = @(q) 0.004*(q - mu)/a;
            [xd,yd] = pf.traceStreamlines(@wZeta,[zF + nrm(zF), zT + nrm(zT)],[-1 1], ...
                H,NSTEP,keep);
            zd = fmap(xd + 1i*yd);
            set(hDiv,'XData',real(zd),'YData',imag(zd));
            set(hZd,'XData',xd,'YData',yd);
            set(hStag,'XData',real(fmap(zF)),'YData',imag(fmap(zF)));
            set(hZstag,'XData',real([zF zT]),'YData',imag([zF zT]));
            set(hZtxt,'String',sprintf(['doublet ' char(954) ' = 2' char(960) 'Ua^2 = %.3f\n' ...
                'vortex ' char(915) ' = %+.3f (Kutta)'],2*pi*U*a^2,Gam));
            set(hCamber,'XData',NaN,'YData',NaN);

            [xu,cpu,xl,cpl,geo] = exactCp();
            th = thinFromGeometry(geo.yt,geo.yc,geo.alphaC);
            set(hExU,'XData',xu,'YData',cpu);  set(hExL,'XData',xl,'YData',cpl);
            hLeg.String(1:2) = {'exact, upper','exact, lower'};
            CLex = 2*(-Gam)/(U*geo.c);
            info = {sprintf('t/c %.3f   camber/c %.3f',2*max(geo.yt),max(abs(geo.yc))), ...
                    sprintf('AoA (chord)= %.2f deg',rad2deg(geo.alphaC)), ...
                    sprintf('zero lift  = %.2f deg',rad2deg(geo.alphaC - aoa - betaK)), ...
                    sprintf('%s (cw)     = %.4f',char(915),-Gam), ...
                    sprintf('CL exact   = %.4f',CLex), ...
                    sprintf('CL thin    = %.4f (x%.3f)',th.CL,CLex/th.CL), ...
                    '', 'potential flow: no drag;', 'boundary layer in thin mode'};
        else
            set([hBL hTr hSep],'XData',NaN,'YData',NaN);
            solveDistributions();
            xb = [xe fliplr(xe)];  yb = [ycE + ytE, fliplr(ycE - ytE)];
            set(hBody,'XData',xb,'YData',yb);
            set(hCamber,'XData',xe,'YData',ycE);
            % viscous-inviscid solution first: it changes the flow field
            if blOn
                if dragPt(1) ~= 0 || quick, nIt = 8; else, nIt = VIT; end
                bl = viscousThin(nIt);
            else
                bl = [];  sigV(:) = 0;  gamV(:) = 0;  sigW(:) = 0;
            end
            if get(cStr,'Value')
                keep = @(z) boxOK(z) & ~insideThin(z);
                [x,y] = pf.traceStreamlines(@wThin,zLeft,ones(1,NLEFT),H,NSTEP,keep);
                set(hS,'XData',x,'YData',y,'Visible','on');
            else
                set(hS,'Visible','off');
            end
            set([hDiv hStag],'XData',NaN,'YData',NaN);
            alphaC = -atan(ycE(end));
            lead = (1 + cos(thm))./sin(thm);
            gA = 2*U*alphaC*lead;
            th = thinCp(gA,gam - gA);
            c = hypot(1,ycE(end));
            G = sum(gam.*dxp);
            set([hExU hExL],'XData',NaN,'YData',NaN);
            CLt = 2*G/(U*c);
            info = {sprintf('t/c %.3f   camber/c %.3f',2*max(ytE)/c, ...
                        max(abs(ycE - ycE(end)*xe))/c), ...
                    sprintf('nose r/c   = %.4f',kS^2*rN/c), ...
                    sprintf('AoA (chord)= %.2f deg',rad2deg(alphaC)), ...
                    sprintf('TE gap/c   = %+.4f',2*ytE(end)/c), ...
                    sprintf('%s (cw)     = %.4f',char(915),G), ...
                    sprintf('CL = 2%s/Uc = %.4f',char(915),CLt)};
            if abs(2*ytE(end)) > 0.002
                info{end+1} = 'body open: tick Keep body closed';
            end
            if min(ytE) < -1e-3
                info{end+1} = 'negative thickness: surfaces cross';
            end
            if blOn && bl.fail
                set([hExU hExL],'Visible','off');
                set([hBL hTr hSep],'XData',NaN,'YData',NaN);
                info{end+1} = 'boundary layer: no solution here';
            elseif blOn
                set(hStag,'XData',real(bl.zStag),'YData',imag(bl.zStag));
                info = [info, showVII()];
                % viscous surface pressure on the "exact" lines
                set(hExU,'XData',xm,'YData',1 - (bl.qu/U).^2,'Visible','on');
                set(hExL,'XData',xm,'YData',1 - (bl.ql/U).^2,'Visible','on');
                hLeg.String(1:2) = {'viscous (VII), upper','viscous (VII), lower'};
            else
                set([hExU hExL],'Visible','off');
            end
            drawEditor();
        end
        set(hInfo,'String',info);

        % thin-theory Cp curves
        set(hThU,'XData',th.x,'YData',th.u);  set(hThL,'XData',th.x,'YData',th.l);
        if get(cComp,'Value')
            set(hCt,'XData',th.x,'YData',th.t);
            set(hCcU,'XData',th.x,'YData',th.cu);  set(hCcL,'XData',th.x,'YData',th.cl);
            set(hCaU,'XData',th.x,'YData',th.au);  set(hCaL,'XData',th.x,'YData',th.al);
        else
            set([hCt hCcU hCcL hCaU hCaL],'XData',NaN,'YData',NaN);
        end

        % velocity vectors
        if get(cVec,'Value')
            [qu,qv] = pf.cappedArrows(physVel(gz),2*U,0.9*DG);
            set(hQ,'UData',qu,'VData',qv,'Visible','on');
        else
            set(hQ,'Visible','off');
        end
        set(hP,'Visible',onoff(get(cPart,'Value')));
    end

    function lines = showVII()
        % draw the coupled solution held in bl and return its info lines
        ds = [bl.dsurf{1}, NaN, bl.dsurf{2}, NaN, bl.wake, NaN, bl.wakeL];
        set(hBL,'XData',real(ds),'YData',imag(ds));
        set(hTr,'XData',real(bl.zTr),'YData',imag(bl.zTr));
        set(hSep,'XData',real(bl.zSep),'YData',imag(bl.zSep));
        CLv = bl.CL;  CLi = bl.CLi;
        if bl.conv
            st = sprintf('VII converged, %d it',bl.it);
        else
            st = sprintf('VII NOT converged (%.0e)',bl.res);
        end
        lines = {sprintf('Re %.3g  %s, %s* x%d',ReN,st,char(948),BLX), ...
                 sprintf('CL viscous %.4f (inv %.4f)',CLv,CLi), ...
                 ['transition up ' sideStr(bl.xTr(1),'laminar') ...
                  ', lo ' sideStr(bl.xTr(2),'laminar')]};
        if bl.sep
            sp = {'',''};
            for k = 1:2
                sp{k} = sideStr(bl.xSep(k),'att.');
                if isfinite(bl.xRe(k)), sp{k} = [sp{k} '-' sprintf('%.2f',bl.xRe(k))]; end
            end
            lines{end+1} = ['separated up ' sp{1} ', lo ' sp{2}];
        end
        lines{end+1} = sprintf('cd (wake) = %.5f   L/D = %.0f',bl.cd,CLv/bl.cd);
    end

    function s = sideStr(x,whenNone)
        % x/c of an event on one surface, or a word when it does not happen
        if isfinite(x), s = sprintf('%.2f',x); else, s = whenNone; end
    end

    function drawEditor()
        % the actual (scaled) densities; dragging edits the shape behind them
        set(hEs,'YData',kS*pchip(XC,sc,xf));  set(hEg,'YData',kG*pchip(XC,gc,xf));
        set(hEsp,'YData',kS*sc);  set(hEgp,'YData',kG*gc(1:end-1));
        ga = 2*aoaT*sqrt((1 - xf)./xf);  ga(1) = NaN;
        set(hEa,'YData',ga);
        sn = kS*sqrt(2*rN)*(1 - 3*xf)./sqrt(xf);  sn(1) = NaN;
        set(hEn,'YData',sn);
    end

    function s = onoff(tf)
        if tf, s = 'on'; else, s = 'off'; end
    end

    % =====================================================================
    %  Particles
    % =====================================================================
    function respawn(idx,atInflow)
        n = nnz(idx);
        x = XL(1) + diff(XL)*rand(1,n);
        if nargin > 1
            in = atInflow(idx);
            x(in) = XL(1) + 0.02*rand(1,nnz(in));
        end
        z = x + 1i*(YL(1) + diff(YL)*rand(1,n));
        if mode == 1, z = invmap(z); end
        pz(idx) = z;  pAge(idx) = 0;
        pLife(idx) = 400 + floor(400*rand(1,n));
    end

    function tick(~,~)
        if ~isgraphics(fig) || busy || ~get(cPart,'Value'), return; end
        try
            if mode == 1
                pz = pf.advectParticles(pz,@wParticleJ,DT,MAXHOP);
                zP = fmap(pz);  inside = abs(pz - mu) < 1.001*a;
            else
                pz = pf.advectParticles(pz,@wThin,DT,MAXHOP);
                zP = pz;  inside = insideThin(pz);
            end
            pAge = pAge + 1;
            out = real(zP) < XL(1) | real(zP) > XL(2) | imag(zP) < YL(1) | imag(zP) > YL(2);
            dead = out | inside | ~isfinite(zP) | pAge > pLife;
            if any(dead), respawn(dead,out); end
            if mode == 1
                zP = fmap(pz);  sp = abs(wZeta(pz)./dmap(pz));
            else
                zP = pz;  sp = abs(wThin(pz));
            end
            set(hP,'XData',real(zP),'YData',imag(zP),'CData',sp);
            drawnow limitrate nocallbacks
        catch err
            if isgraphics(fig), rethrow(err); end
        end
    end


    % =====================================================================
    %  Controls
    % =====================================================================
    function setMode(m)
        mode = m;
        set(hMode,'Value',m);
        isJ = mode == 1;
        set(pnlZ,'Visible',onoff(isJ));  set(pnlE,'Visible',onoff(~isJ));
        set([hExU hExL],'Visible',onoff(isJ));
        if isJ, hLeg.String(1:2) = {'exact, upper','exact, lower'}; end
        configureSliders();
        configureBottom();
        set([cBL hRe],'Enable',onoff(~isJ));     % the boundary layer is thin mode only
        if isJ
            set(hHint,'String',['Circle plane (bottom right): freestream + doublet + ' ...
                'vortex, mapped onto the airfoil. Solid Cp exact, dashed thin theory. ' ...
                'Potential flow only: the boundary layer and viscous polar are in thin mode.']);
        else
            set(hHint,'String',['Drag red source points (thickness) and green vortex ' ...
                'points (camber). Sliders scale them linearly; AoA adds the exact ' ...
                'flat-plate sheet, the nose slider an exact round nose. Boundary layer ' ...
                'uses the Lighthill-Van Dyke leading-edge correction.']);
        end
        readSliders();
        if ~isempty(pz), respawn(true(1,NPART)); end
    end

    function configureSliders()
        % labels, ranges and values of the three sliders for the current mode
        if mode == 1
            cfg = {['thickness  ' char(949) '/b'],[0 0.3],thk;
                   ['camber  ' char(948) '/b'],[-0.15 0.3],cmb;
                   'angle of attack [deg]',[-10 15],rad2deg(aoa)};
        else
            cfg = {['thickness scale (' char(963) ')'],[0 3],kS;
                   ['camber scale (' char(947) ')'],[-2 3],kG;
                   'angle of attack [deg]',[-10 15],rad2deg(aoaT)};
        end
        hs = [hThk hCmb hAoa];  hl = [hL1 hL2 hL3];
        for k = 1:3
            set(hl(k),'String',cfg{k,1});
            set(hs(k),'Min',-1e9,'Max',1e9);        % widen first: no range warning
            set(hs(k),'Value',cfg{k,3},'Min',cfg{k,2}(1),'Max',cfg{k,2}(2));
        end
    end

    function readSliders(fast)
        if nargin < 1, fast = false; end
        quick = fast;
        v = [get(hThk,'Value') get(hCmb,'Value') get(hAoa,'Value')];
        if mode == 1
            thk = round(v(1)*1000)/1000;
            cmb = round(v(2)*1000)/1000;
            aoa = deg2rad(round(v(3)*10)/10);
            s = {sprintf('%.3f',thk), sprintf('%.3f',cmb), sprintf('%.1f',rad2deg(aoa))};
        else
            kS = round(v(1)*100)/100;
            kG = round(v(2)*100)/100;
            aoaT = deg2rad(round(v(3)*10)/10);
            s = {sprintf('x %.2f',kS), sprintf('x %.2f',kG), sprintf('%.1f',rad2deg(aoaT))};
        end
        set(hThkV,'String',s{1});  set(hCmbV,'String',s{2});  set(hAoaV,'String',s{3});
        redraw();
        quick = false;
    end

    function apiSetJoukowski(t,c,adeg)
        if mode ~= 1, setMode(1); end
        set(hThk,'Value',t);  set(hCmb,'Value',c);  set(hAoa,'Value',adeg);
        readSliders();
    end

    function setBL()
        v = str2double(get(hRe,'String'));
        if isfinite(v) && v >= 1e4, ReN = v; end
        set(hRe,'String',sprintf('%.3g',ReN));
        blOn = logical(get(cBL,'Value'));
        redraw();
    end

    function configureBottom()
        % the small slider: trailing-edge angle (Joukowski) or nose radius (thin)
        set(hTE,'Min',-1e9,'Max',1e9);
        if mode == 1
            set(hTE,'Value',tauTE,'Min',0,'Max',20,'SliderStep',[0.025 0.1]);
        else
            set(hTE,'Value',rN,'Min',0,'Max',0.04,'SliderStep',[0.025 0.1]);
        end
        bottomLabel();
    end

    function bottomLabel()
        if mode == 1
            set(hTEl,'String',sprintf('TE angle %.1f%s',tauTE,char(176)));
        else
            set(hTEl,'String',sprintf('nose radius r/c %.4f',rN));
        end
    end

    function setBottom(fast)
        if nargin < 1, fast = false; end
        quick = fast;
        if mode == 1
            tauTE = round(get(hTE,'Value')*2)/2;
        else
            rN = round(get(hTE,'Value')*1e4)/1e4;
        end
        bottomLabel();
        redraw();
        quick = false;
    end

    function apiSetBL(on,re,bottom)
        % bottom: TE angle [deg] in Joukowski mode, nose radius in thin mode
        set(cBL,'Value',on);  set(hRe,'String',num2str(re));
        if nargin > 2
            if mode == 1, tauTE = bottom; else, rN = bottom; end
            configureBottom();
        end
        setBL();
    end

    function apiSetThin(ks,kg,adeg)
        if mode ~= 2, setMode(2); end
        set(hThk,'Value',ks);  set(hCmb,'Value',kg);  set(hAoa,'Value',adeg);
        readSliders();
    end

    function apiSetControls(s,g)
        sc = s;  gc = g;  gc(end) = 0;
        if get(cClose,'Value'), keepClosed(0); end
        redraw();
    end

    function preset(k)
        % NACA 0012 = exact round-nose term (r/c = 0.0159, which carries its
        % sqrt(x) part) + a smooth remainder for the editor:
        %   dy/dx(NACA) - dy/dx(nose) = 0.6(-0.126 - 0.7032x + 0.8529x^2
        %                               - 0.4144x^3) + 1.5*0.17814*sqrt(x)
        r0 = (0.6*0.2969)^2/2;
        dyt = 0.6*(-0.126 - 0.7032*XC + 0.8529*XC.^2 - 0.4144*XC.^3) ...
            + 1.5*sqrt(2*r0)*sqrt(XC);
        kS = 1;  kG = 1;  aoaT = 0;  rN = r0;
        switch k
            case 1
                sc = 2*dyt;  gc = zeros(size(XC));
            case 2
                sc = 2*dyt;  gc = 16*0.04*sqrt(XC.*(1 - XC));   % parabolic camber, 4 %
            otherwise
                sc = zeros(size(XC));  gc = zeros(size(XC));    % flat plate ...
                aoaT = deg2rad(5);  rN = 0;                      % ... at 5 deg, sharp
        end
        gc(end) = 0;
        if get(cClose,'Value'), keepClosed(0); end
        showThin();
    end

    function keepClosed(fixedIdx)
        % Make the net source strength zero by shifting every source point
        % by the same amount, except the dragged one (FIXEDIDX) and the
        % leading-edge point, which shapes the nose. PCHIP is not linear in
        % the point values, so the shift is refined a few times.
        free = true(size(XC));  free(1) = false;
        if fixedIdx > 0, free(fixedIdx) = false; end
        xq = linspace(0,1,2001);
        Ifree = trapz(xq,pchip(XC,double(free),xq));
        for it = 1:4
            I = trapz(xq,pchip(XC,sc,xq));
            sc(free) = sc(free) - I/Ifree;
        end
    end

    function toggleClosed()
        if get(cClose,'Value'), keepClosed(0); end
        redraw();
    end

    function showThin()
        % the distribution buttons work from either mode
        if mode == 2
            configureSliders();
            configureBottom();
            readSliders();
        else
            setMode(2);
        end
    end

    % =====================================================================
    %  Polar
    % =====================================================================
    function askPolar()
        % ask for the angle range, then sweep
        if mode == 2, def = {'-4','16','1'}; else, def = {'-4','12','1'}; end
        ans_ = inputdlg({'from [deg]','to [deg]','step [deg]'},'Polar',1,def);
        if isempty(ans_), return, end
        v = str2double(ans_);
        if any(~isfinite(v)) || v(3) <= 0 || v(2) <= v(1), return, end
        makePolar(v(1):v(3):v(2));
    end

    function T = makePolar(alist)
        % Sweep the angle of attack for the airfoil as it is set up now and
        % plot the polar. Thin mode: the coupled boundary layer at the
        % current Re (on even if its checkbox is off), stepping away from
        % 0 deg in both directions so that each solution starts from its
        % neighbour. Joukowski mode: exact potential flow (lift only).
        alist = sort(alist(:).');
        if mode == 2
            aSave = aoaT;  vsSave = vs;  blSave = blOn;
            [~,i0] = min(abs(alist));
            order = [i0:numel(alist), i0-1:-1:1];
            R = nan(numel(alist),9);
            vs = [];  vsStart = [];
            for q = 1:numel(order)
                k = order(q);
                if k == i0 - 1, vs = vsStart; end      % downward branch starts at the middle
                aoaT = deg2rad(alist(k));
                solveDistributions();
                set(hInfo,'String',{sprintf('polar: %g deg  (%d of %d)',alist(k),q,numel(alist))});
                drawnow limitrate
                o = viscousThin(VIT);
                if k == i0, vsStart = vs; end
                c = hypot(1,ycE(end));
                CLi = 2*sum(gam.*dxp)/(U*c);
                if o.fail
                    R(k,:) = [alist(k), rad2deg(-atan(ycE(end))), CLi, NaN(1,5), 0];
                else
                    R(k,:) = [alist(k), rad2deg(-atan(ycE(end))), o.CLi, o.CL, o.cd, ...
                        o.xTr(1), o.xTr(2), o.xSep(1), o.conv || o.res < 0.02];
                end
            end
            aoaT = aSave;  vs = vsSave;  blOn = blSave;
            T = array2table(R,'VariableNames',{'alpha','alphaChord','CLinv','CL','cd', ...
                'xtrUpper','xtrLower','xsepUpper','converged'});
        else
            aSave = aoaT;  aJ = aoa;
            R = nan(numel(alist),4);
            for k = 1:numel(alist)
                aoa = deg2rad(alist(k));
                setJoukowski();
                [~,~,~,~,geo] = exactCp();
                th = thinFromGeometry(geo.yt,geo.yc,geo.alphaC);
                R(k,:) = [alist(k), rad2deg(geo.alphaC), 2*(-Gam)/(U*geo.c), th.CL];
            end
            aoa = aJ;  aoaT = aSave;
            T = array2table(R,'VariableNames',{'alpha','alphaChord','CL','CLthin'});
        end
        assignin('base','profilPolar',T);
        showPolar(T);
        redraw();
    end

    function showPolar(T)
        pf_ = findobj(0,'Type','figure','Tag','profilPolar');
        if isempty(pf_)
            pf_ = figure('Tag','profilPolar','Name','Profil - polar','NumberTitle','off', ...
                'Color','w','Position',[120 80 1100 760]);
            try theme(pf_,'light'); catch, end   %#ok<*NOSEMI>
        end
        clf(pf_);
        cV = [0.13 0.33 0.62];  cI = [0.55 0.55 0.55];
        if ismember('cd',T.Properties.VariableNames)
            ok = T.converged == 1;
            ttl = sprintf('Thin airfoil, Re = %.3g  (open markers: not converged)',ReN);
            mk = @(ax,x,y) [plot(ax,x,y,'-','Color',cV,'LineWidth',1.4), ...
                plot(ax,x(ok),y(ok),'o','Color',cV,'MarkerFaceColor',cV,'MarkerSize',5), ...
                plot(ax,x(~ok),y(~ok),'o','Color',cV,'MarkerFaceColor','w','MarkerSize',5)];
            ax1 = subplot(2,2,1,'Parent',pf_);  hold(ax1,'on');  grid(ax1,'on');
            plot(ax1,T.alphaChord,T.CLinv,'--','Color',cI,'LineWidth',1.2);
            mk(ax1,T.alphaChord,T.CL);
            xlabel(ax1,'\alpha (chord) [deg]');  ylabel(ax1,'C_L');
            legend(ax1,{'inviscid','viscous'},'Location','northwest');
            ax2 = subplot(2,2,2,'Parent',pf_);  hold(ax2,'on');  grid(ax2,'on');
            mk(ax2,T.cd,T.CL);  xlabel(ax2,'c_d');  ylabel(ax2,'C_L');  title(ax2,'drag polar');
            ax3 = subplot(2,2,3,'Parent',pf_);  hold(ax3,'on');  grid(ax3,'on');
            mk(ax3,T.alphaChord,T.cd);  xlabel(ax3,'\alpha (chord) [deg]');  ylabel(ax3,'c_d');
            ax4 = subplot(2,2,4,'Parent',pf_);  hold(ax4,'on');  grid(ax4,'on');
            mk(ax4,T.alphaChord,T.CL./T.cd);  xlabel(ax4,'\alpha (chord) [deg]');  ylabel(ax4,'L/D');
            sgtitle(pf_,ttl);
        else
            ax1 = axes('Parent',pf_);  hold(ax1,'on');  grid(ax1,'on');
            plot(ax1,T.alphaChord,T.CL,'-o','Color',cV,'MarkerFaceColor',cV,'LineWidth',1.4,'MarkerSize',5);
            plot(ax1,T.alphaChord,T.CLthin,'--','Color',cI,'LineWidth',1.2);
            xlabel(ax1,'\alpha (chord) [deg]');  ylabel(ax1,'C_L');
            legend(ax1,{'exact (Joukowski, Kutta)','thin airfoil theory'},'Location','northwest');
            title(ax1,'Joukowski airfoil: lift (potential flow, no drag)');
        end
        uicontrol(pf_,'Style','pushbutton','String','Save CSV...','Units','normalized', ...
            'Position',[0.86 0.005 0.12 0.04],'Callback',@(~,~) savePolar(T));
    end

    function savePolar(T)
        [fn,fp] = uiputfile('*.csv','Save polar','polar.csv');
        if isequal(fn,0), return, end
        writetable(T,fullfile(fp,fn));
    end

    function s = getState()
        s = struct('mode',mode,'mu',mu,'a',a,'beta',betaK,'Gamma',Gam, ...
            'sigma',sig,'gamma',gam,'yt',ytE,'yc',ycE,'x',xe,'info',{info}, ...
            'bl',bl,'nKT',nKT,'rN',rN);
    end

    function out = getCp()
        out = struct('exU',[get(hExU,'XData'); get(hExU,'YData')], ...
            'exL',[get(hExL,'XData'); get(hExL,'YData')], ...
            'thU',[get(hThU,'XData'); get(hThU,'YData')], ...
            'thL',[get(hThL,'XData'); get(hThL,'YData')]);
    end

    % =====================================================================
    %  Mouse: dragging editor points
    % =====================================================================
    function onDown(~,~)
        dragPt = [0 0];
        if mode ~= 2, return; end
        cp = get(axE,'CurrentPoint');  p = cp(1,1:2);
        yl = get(axE,'YLim');
        if p(1) < 0 || p(1) > 1 || p(2) < yl(1) || p(2) > yl(2), return; end
        ys = diff(yl);
        ds = hypot(XC - p(1),(kS*sc - p(2))/ys);
        dg = hypot(XC(1:end-1) - p(1),(kG*gc(1:end-1) - p(2))/ys);
        [m1,i1] = min(ds);  [m2,i2] = min(dg);
        if min(m1,m2) > PICK, return; end
        if m1 <= m2, dragPt = [1 i1]; else, dragPt = [2 i2]; end
    end

    function onMove(~,~)
        if dragPt(1) == 0 || busy, return; end
        busy = true;
        try
            cp = get(axE,'CurrentPoint');
            yl = get(axE,'YLim');
            v = min(max(cp(1,2),yl(1)),yl(2));
            % a point at scale 0 cannot be moved: there is nothing to scale
            if dragPt(1) == 1 && abs(kS) > 1e-3
                sc(dragPt(2)) = v/kS;
                if get(cClose,'Value'), keepClosed(dragPt(2)); end
            elseif dragPt(1) == 2 && abs(kG) > 1e-3
                gc(dragPt(2)) = v/kG;
            end
            redraw();
            drawnow limitrate
        catch err
            busy = false;
            rethrow(err);
        end
        busy = false;
    end

    function onUp(~,~)
        if dragPt(1) == 0, return; end
        dragPt = [0 0];
        redraw();
    end
end
