function [J,Un,xn,beta0x] = betaSolve_2009(beta0,H,x,U,hb,n,E,m,dx0,rho_i,g,rho_sw,sigma_b,dUdx,c0,x0,hb0,W0,A0,U0,h0)
    
    % interpolate beta0 to grid spacing
    beta0x = 0:round(x0(end)/(length(beta0))):round(x0(end)/(length(beta0)))*(length(beta0)-1); 
    beta = interp1(beta0x,beta0,x0,'pchip');
    
    % use initial calving front position
    xcf = x0(c0); 
    
    % calculate the thickness required to remain grounded at each grid cell
    Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
    % find the location of the grounding line and use a floating
    % geometry from the grounding line to the calving front
    if ~isempty(find(Hf-H>0,1,'first'))
        if length(Hf)>=find(Hf-H>0,1,'first')+1
            xgl = interp1(Hf(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1)...
                -H(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),...
                x(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
        else
            xgl = x(find(Hf-H>0,1,'first')-1);
        end
    else
        xgl=xcf;
    end
    if xgl>xcf % grounding line can't be past calving front
        xgl=xcf;
    end

    % create coordinate system that hits cf and gl exactly
    % has resolution dxmax near the ice divide
    % has resolution dxmin from gl to c
    % and has smooth variation between
    xl = round(xgl/dx0); %number of ideal grid spaces needed to reach the grounding line
    dx = xgl/xl; %new grid spacing (should be ~dx0)
    xn = 0:dx:xgl; %new distance vector
    if xcf-xgl > 0
        xl = round((xcf-xgl)/dx0);
        dx = (xcf-xgl)/xl;
        xn = [xn xn(end)+dx:dx:xcf];
    end
    clear dx; dxn = [xn(2:end)-xn(1:end-1) xn(end)-xn(end-1)];

    % get geometry on new coordinates
    c = length(xn); gl = dsearchn(xn',xgl); % indeces for cf and gl
    %if the crevasses never intersect sea level
    if isempty(c) == 1 %set the calving front to a default minimum ice thickness value
        c = find(H<Hc,1,'first');
    end
    if isempty(c)==1 % set c to length of x if criteria still not met
        c=length(x);
        disp('calving criteria not met');
    end

    hb = interp1(x0,hb0,xn,'linear','extrap');
    W = interp1(x0,W0,xn,'linear','extrap');
    H = interp1(x,H,xn,'linear','extrap');
    U = interp1(x,U,xn,'linear','extrap');
    dUdx = interp1(x,dUdx,xn,'linear','extrap');
    A = interp1(x0,A0,xn,'linear','extrap');
    x = xn; dx = dxn;
    beta(gl+1:end) = 0; % zero bed roughness where ice is floating

    % calculate surface elevation
    h = hb+H; % surface elevation (m a.s.l.)
    h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); %adjust the surface elevation of ungrounded ice to account for buoyancy
    H(h<0)=0-hb(h<0); h(h<0)=0; % surface cannot go below sea level
    h(h-H<hb) = hb(h-H<hb)+H(h-H<hb); % thickness cannot go beneath bed elevation

    % calculate the effective pressure (ice overburden pressure minus water
    % pressure) assuming an easy & open connection between the ocean and
    % ice-bed interface
    sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
    N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
    N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
    N = [N_ground N_marine];
    N(N<0)=0; % cannot have negative values
    
    %calculate the linearization terms & effective viscosity required for 
    %inversion of the stress coefficient matrix
    dUdx(1) = (U(2)-U(1))./(x(2)-x(1)); % forward difference
    dUdx(2:c-1) = (U(3:c)-U(1:c-2))./(x(3:c)-x(1:c-2)); % central difference
    dUdx(c) = (U(c-1)-U(c))/(x(c-1)-x(c)); % backward difference at c
       
    % solve for new velocity
    [Un,~,~,~,~] = U_convergence(x,U,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,1);
    xn=x;
    
    % calculate cost of parameter solutions
    % modified from Morlighem et al., 2010; Larour et al., 2012; Kyrke-Smith et al., 2018
    U_err = 33/3.1536e7; % m/s
    h_err = 22; % m
    K = log(beta); K(~isfinite(K))=0;
    J = nanmean(sqrt((Un-interp1(x0,U0,x)).^2)-U_err)/(nanmean(U0))+... % speed misfit term
        nanmean(abs(gradient(gradient(K)))); % regularization term to penalize changes in beta gradient
        %nanmean(sqrt((h-interp1(x0,h0,x)).^2)-h_err)./nanmean(h0)+... % surface elevation misfit term

end