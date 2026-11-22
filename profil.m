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
%     dragged, so that int sigma = 0.
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

    % ---- state ----------------------------------------------------------
    mode = 1;                               % 1 Joukowski, 2 thin airfoil
    thk = 0.12;  cmb = 0.06;  aoa = deg2rad(5);
    mu = 0;  a = b;  betaK = 0;  Gam = 0;  Uc = U;  Ap = 0;  Bp = 0;
    sc = zeros(size(XC));  gc = zeros(size(XC));   % editor shapes (before scaling)
    kS = 1;  kG = 1;  aoaT = 0;             % thin mode: thickness/camber scale, AoA
    sig = zeros(1,NPAN);  gam = sig;  ytE = zeros(1,NPAN+1);  ycE = ytE;
    info = {};
    pz = zeros(1,0);  pAge = pz;  pLife = pz;
    dragPt = [0 0];  busy = false;

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
    legend(axC,[hExU hExL hThU hThL hCt hCcU hCaU], ...
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
    hEsp = line('Parent',axE,'XData',XC,'YData',sc,'LineStyle','none','Marker','o','MarkerSize',8, ...
        'MarkerFaceColor',[0.84 0.22 0.18],'MarkerEdgeColor','k',np{:});
    hEgp = line('Parent',axE,'XData',XC(1:end-1),'YData',gc(1:end-1),'LineStyle','none', ...
        'Marker','s','MarkerSize',8,'MarkerFaceColor',[0.2 0.6 0.3],'MarkerEdgeColor','k',np{:});
    line('Parent',axE,'XData',1,'YData',0,'LineStyle','none','Marker','s','MarkerSize',8, ...
        'MarkerEdgeColor',[0.2 0.6 0.3],'LineWidth',1.5,np{:});
    text(axE,0.70,2.75,[char(947) '(1) = 0: Kutta condition'],'FontSize',9,np{:});

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
    for h = [hThk hCmb hAoa]
        set(h,'Callback',@(~,~) readSliders());
        addlistener(h,'ContinuousValueChange',@(~,~) readSliders());
    end

    btn = @(s,x,w,cb) uicontrol(fig,'Style','pushbutton','String',s,'Units','normalized', ...
        'Position',[x 0.105 w 0.04],'FontSize',10,'Callback',cb);
    btn('Symmetric',0.04,0.075,@(~,~) preset(1));
    btn('Cambered',0.12,0.075,@(~,~) preset(2));
    btn('Plate, 5 deg',0.20,0.075,@(~,~) preset(3));

    hInfo = uicontrol(fig,'Style','text','Units','normalized','Position',[0.41 0.10 0.19 0.19], ...
        'BackgroundColor',[0.97 0.97 0.97],'FontSize',9,'HorizontalAlignment','left', ...
        'FontName','Consolas');
    hHint = uicontrol(fig,'Style','text','Units','normalized','Position',[0.17 0.015 0.43 0.08], ...
        'BackgroundColor','w','ForegroundColor',[0.35 0.35 0.35],'FontSize',9, ...
        'HorizontalAlignment','left');

    set(fig,'WindowButtonDownFcn',@onDown,'WindowButtonMotionFcn',@onMove, ...
        'WindowButtonUpFcn',@onUp);

    setappdata(fig,'api',struct('setMode',@setMode,'setJoukowski',@apiSetJoukowski, ...
        'setThin',@apiSetThin,'setControls',@apiSetControls,'preset',@preset, ...
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
        % circle through the trailing-edge point zeta = b, Kutta circulation
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
        z = zeta + b^2./zeta + X0;
    end

    function d = dmap(zeta)
        d = 1 - b^2./zeta.^2;
    end

    function [zeta,inside] = invmap(z)
        % the preimage outside the circle; points inside the airfoil have none
        zz = z - X0;
        s = sqrt(zz.^2 - 4*b^2);
        r1 = (zz + s)/2;  r2 = (zz - s)/2;
        use2 = abs(r2 - mu) > abs(r1 - mu);
        zeta = r1;  zeta(use2) = r2(use2);
        inside = abs(zeta - mu) < a;
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
            'alphaC',aoa - angle(eC),'zLE',zLE,'zTE',zTE);
    end

    % =====================================================================
    %  Thin airfoil theory
    % =====================================================================
    function th = thinFromGeometry(yt,yc,alphaC)
        % sources from the thickness slope, vortices from Glauert's solution
        sig = 2*U*diff(yt)./dxp;
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
        lead = (1 + cos(thm))./sin(thm);
        sig = kS*U*pchip(XC,sc,xm);
        gUser = kG*U*pchip(XC,gc,xm);
        gam = gUser + 2*U*aoaT*lead;
        ytE = [0 cumsum(sig.*dxp)]/(2*U);
        % downwash on the chord line: panels for the editor part, and the
        % exact uniform -U*alpha of the flat-plate sheet (panels would smear
        % its leading-edge singularity)
        v = -(Lmat*gUser.').'/(2*pi) - U*aoaT;
        ycE = [0 cumsum(v.*dxp)]/U;
    end

    function th = thinCp(gA,gC)
        ut = (Lmat*sig.').'/(2*pi);
        th.x = xm;
        th.t = -2*ut/U;
        th.au = -gA/U;  th.al = gA/U;
        th.cu = -gC/U;  th.cl = gC/U;
        th.u = th.t + th.au + th.cu;
        th.l = th.t + th.al + th.cl;
    end

    function w = wThin(z)
        % source and vortex panels on the chord line, constant density each
        sz = size(z);  z = z(:);
        s = ((sig + 1i*gam)/(2*pi)).';
        w = U + log((z - xe(1:end-1))./(z - xe(2:end)))*s;
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
            % body and circle
            zc = mu + a*exp(1i*linspace(0,2*pi,400));
            zb = fmap(zc);
            set(hBody,'XData',real(zb),'YData',imag(zb));
            set(hZc,'XData',real(zc),'YData',imag(zc));
            set(hZmu,'XData',real(mu),'YData',imag(mu));
            % streamlines are traced in the circle plane and mapped
            if get(cStr,'Value')
                keep = @(zeta) abs(zeta - mu) > 1.002*a & boxOK(fmap(zeta));
                [xz,yz] = pf.traceStreamlines(@wZeta,invmap(zLeft),ones(1,NLEFT),H,NSTEP,keep);
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
            CLex = 2*(-Gam)/(U*geo.c);
            info = {sprintf('t/c        = %.3f',2*max(geo.yt)), ...
                    sprintf('camber/c   = %.3f',max(abs(geo.yc))), ...
                    sprintf('AoA (chord)= %.2f deg',rad2deg(geo.alphaC)), ...
                    sprintf('zero lift  = %.2f deg',rad2deg(geo.alphaC - aoa - betaK)), ...
                    sprintf('%s (cw)     = %.4f',char(915),-Gam), ...
                    sprintf('CL exact   = %.4f',CLex), ...
                    sprintf('CL thin    = %.4f',th.CL), ...
                    sprintf('exact/thin = %.3f',CLex/th.CL)};
        else
            solveDistributions();
            xb = [xe fliplr(xe)];  yb = [ycE + ytE, fliplr(ycE - ytE)];
            set(hBody,'XData',xb,'YData',yb);
            set(hCamber,'XData',xe,'YData',ycE);
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
            info = {sprintf('t/c        = %.3f',2*max(ytE)/c), ...
                    sprintf('camber/c   = %.3f',max(abs(ycE - ycE(end)*xe))/c), ...
                    sprintf('AoA (chord)= %.2f deg',rad2deg(alphaC)), ...
                    sprintf('TE gap/c   = %+.4f',2*ytE(end)/c), ...
                    sprintf('%s (cw)     = %.4f',char(915),G), ...
                    sprintf('CL = 2%s/Uc = %.4f',char(915),2*G/(U*c)), ...
                    '', ...
                    'thin theory only'};
            if abs(2*ytE(end)) > 0.002
                info{7} = 'body open: tick Keep body closed';
            end
            if min(ytE) < -1e-3
                info{8} = 'negative thickness: surfaces cross';
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

    function drawEditor()
        % the actual (scaled) densities; dragging edits the shape behind them
        set(hEs,'YData',kS*pchip(XC,sc,xf));  set(hEg,'YData',kG*pchip(XC,gc,xf));
        set(hEsp,'YData',kS*sc);  set(hEgp,'YData',kG*gc(1:end-1));
        ga = 2*aoaT*sqrt((1 - xf)./xf);  ga(1) = NaN;
        set(hEa,'YData',ga);
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
        configureSliders();
        if isJ
            set(hHint,'String',['Joukowski: the circle plane (bottom right) holds a ' ...
                'freestream, a doublet and a vortex, exactly like przeplywer. ' ...
                'The map z = ' char(950) ' + b' char(178) '/' char(950) ...
                ' turns the circle into the airfoil; red + marks the critical points ' ...
                char(177) 'b. Solid Cp: exact; dashed: thin theory for the same shape. ' ...
                'To shape an airfoil with sources and vortices, pick a preset or switch Mode.']);
        else
            set(hHint,'String',['Thin airfoil: drag the red source points (thickness) and ' ...
                'green vortex points (camber). The sliders scale the two distributions ' ...
                'linearly, and the angle of attack adds the exact flat-plate vortex sheet, ' ...
                'so the airfoil turns nose-up in the horizontal freestream. Cp is linear ' ...
                'theory: the three parts add up exactly to the dashed total.']);
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

    function readSliders()
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
    end

    function apiSetJoukowski(t,c,adeg)
        if mode ~= 1, setMode(1); end
        set(hThk,'Value',t);  set(hCmb,'Value',c);  set(hAoa,'Value',adeg);
        readSliders();
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
        xx = max(XC,0.004);
        dyt = 0.6*(0.2969*0.5./sqrt(xx) - 0.126 - 2*0.3516*xx + 3*0.2843*xx.^2 ...
            - 4*0.1036*xx.^3);                   % NACA 0012 thickness slope
        kS = 1;  kG = 1;  aoaT = 0;
        switch k
            case 1
                sc = 2*dyt;  gc = zeros(size(XC));
            case 2
                sc = 2*dyt;  gc = 16*0.04*sqrt(XC.*(1 - XC));   % parabolic camber, 4 %
            otherwise
                sc = zeros(size(XC));  gc = zeros(size(XC));    % flat plate ...
                aoaT = deg2rad(5);                               % ... at 5 deg
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
            readSliders();
        else
            setMode(2);
        end
    end

    function s = getState()
        s = struct('mode',mode,'mu',mu,'a',a,'beta',betaK,'Gamma',Gam, ...
            'sigma',sig,'gamma',gam,'yt',ytE,'yc',ycE,'x',xe,'info',{info});
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
