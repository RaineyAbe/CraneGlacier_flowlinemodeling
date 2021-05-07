function [J,beta,x,U,H,h] = betaSolve(H,x,U,hb,n,E,fwd,m,dx0,rho_i,g,rho_sw,rho_fw,sigma_b,dUdx,c0,x0,hb0,W0,A0,U_2017,h_2017,xcf_2017,smr0,smb0,Q0)
      
try

    % used 2009 observed calving front position
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
        xl = round((xcf-xgl)/dx);
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
    hb = interp1(x0,hb0,xn,'pchip');
    W = interp1(x0,W0,xn,'pchip');
    H = interp1(x,H,xn,'pchip');
    U = interp1(x,U,xn,'pchip');
    dUdx = interp1(x,dUdx,xn,'pchip');
    A = interp1(x0,A0,xn,'pchip');
    x = xn; dx = dxn;

    % calculate surface elevation
    h = hb+H; % surface elevation (m a.s.l.)
    h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); %adjust the surface elevation of ungrounded ice to account for buoyancy
    H(h<0)=0-hb(h<0); h(h<0)=0; % surface cannot go below sea level

    % calculate the effective pressure (ice overburden pressure minus water
    % pressure) assuming an easy & open connection between the ocean and
    % ice-bed interface
    sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
    N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
    N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
    N = [N_ground N_marine];
    N(N<0)=0; % cannot have negative values

    % Set up H & dHdx on staggered grid
    Hm(1:c-1) = (H(2:c) + H(1:c-1))./2; % forward difference
    Hm(c) = (H(c)+H(c-1))./2; % backward difference at c
    dHmdx(1:c-1) = (H(2:c)-H(1:c-1))./(x(2:c)-x(1:c-1)); % forward difference
    dHmdx(c) = (H(c)-H(c-1))./(x(c)-U(c-1)); % backward difference at c    

    %calculate the linearization terms & effective viscosity required for 
    %inversion of the stress coefficient matrix
    if n == 3
        gamma=zeros(1,c); % pre-allocate gamma
        for k=1:c
            gamma(k) = U(k).^((1-n)/n); % linearization term for lateral resistance
        end
        gamma(1) = gamma(2); % set linearization term at the divide (U(1) = 0)
        gamma(gamma>1e+06) = 1e+06; % set the limit so gamma does not approach infinity (minimum U = 1e-09 m s^-1)

        % get A, U, & the effective viscosity on the staggered grid for the 
        % longitudinal stress calculation
        Am(1:c-1) = (A(2:c)+A(1:c-1))./2; % forward difference
        Am(c) = (A(c)+A(c-1))./2; % backward difference at c

        Um(1:c-1) = (U(2:c)+U(1:c-1))./2; % forward difference
        Um(c) = (U(c)+U(c-1))./2; % backward difference at c

        dUmdx(1:c-1) = (U(1:c-1)-U(2:c))/(x(1:c-1)-x(2:c)); % forward difference
        dUmdx(c) = (U(c-1)-U(c))/(x(c-1)-x(c)); % backward difference at c

        dUdx(1) = (U(2)-U(1))./(x(2)-x(1)); % forward difference
        dUdx(2:c-1) = (U(3:c)-U(1:c-2))./(x(3:c)-x(1:c-2)); % central difference
        dUdx(c) = (U(c-1)-U(c))/(x(c-1)-x(c)); % backward difference at c

        vm = ((E.*Am).^(-1/n)).*(abs(dUmdx)).^((1-n)/n);
        vm(vm>8e+16) = 8e+16; %set a maximum value for very low strain rates

        if m > 1
            eta=zeros(1,c); % pre-allocate eta
            for k=1:c
                eta(k) = U(k).^((1-m)/m); %linearization term for basal resistance
            end
            eta(1) = eta(2); %set linearization term at the divide (U(1) = 0)

            %set the limit so eta does not approach infinity (minimum U = 1e-09 m s^-1)
            if m == 2
                eta(eta>3.16e+04) = 3.16e+04;
            end
            if m == 3
                eta(eta>1e+06) = 1e+06;
            end
        else
            eta = ones(1,c); %if m=1, the basal resistance term does not need to be linearized
        end
    else
        disp(['Adjust maximum value for the lateral resistance linearization term (gamma)']);
    end

    %set-up coefficient vectors for the linearized stress terms over the calving front
    %[C(k)*U(k-1)+E(k)*U(k)+G(k)*U(k+1)=Td]  
    % coefficients up to calving front
    G_minus(2:c-1) = (2./(dx(2:c-1).^2)).*Hm(1:c-2).*vm(1:c-2); %for U(k-1)
    G_plus(2:c-1) = (2./(dx(2:c-1).^2)).*Hm(2:c-1).*vm(2:c-1); %for U(k+1)
    T(2:c-1) = (rho_i.*g.*H(2:c-1).*(h(1:c-2)-h(3:c))./(x(1:c-2)-x(3:c))); %gravitational driving stress
    % upper boundary condition
    T(1) = (rho_i.*g.*H(1).*(h(1)-h(2))./(x(1)-x(2)));      
    % calving front condition
    G_minus(c) = -1;
    G_plus(c) = 0;
    T(c) = (E*A(c).*(((rho_i.*g./4).*((H(c).*(1-(rho_i./rho_sw))-sigma_b./(rho_i.*g)))).^n)).*dx(c); 
    %remove any NaNs from the coefficient vectors
    G_minus(isnan(G_minus)) = 0;
    G_plus(isnan(G_plus)) = 0;
    T(isnan(T)) = 0;

    % Solve for beta using the G term (solved from the equation below)
    % [G_minus(k)*U(k-1)+G(k)*U(k)+G_plus(k)*U(k+1)=T(k)]
    % G(k) = (T(k) - G_minus(k)*U(k-1) - G_plus(k)*U(k+1)))/U(k)
    %   where G(k) = -2./(dx(k).^2)*(Hm(k)*vm(k)+Hm(k-1)*vm(k-1)) 
    %             - (beta(k)*N(k)*eta(k))
    %             - (gamma(k)*H(k)/W(k))*(5/(2*A(k)*W(k))^(1/3)); 
    %
    %        beta(k) = [-G(k) - 2./(dx(k).^2)*(Hm(k)*vm(k)+Hm(k-1)*vm(k-1)) 
    %         - (gamma(k)*H(k)/W(k))*(5/(2*A(k)*W(k))^(1/3))]/(N(k)*eta(k));
    % Solve for G (for U(k))
    G(2:c) = (T(2:c)-G_minus(2:c).*Um(1:c-1)-G_plus(2:c).*Um(2:c))./U(2:c);
    G(1) = G(2);

    % Solve the basal roughness factor, beta
    beta(2:c) = (-G(2:c)-(2./(dx(2:c).^2)).*(Hm(2:c).*vm(2:c)+Hm(1:c-1).*vm(1:c-1))...
        -(2.*gamma(2:c).*H(2:c)./W(2:c)).*((5./(A(2:c).*W(2:c))).^(1/n)))./(N(2:c).*eta(2:c));  
    beta(1) = beta(2); 
    % average beta over the approximate stress-coupling length (SCL)
%     scl = 1350; % m 
%     edges = (0:scl:x(end));
%     xmid = 0.5*(edges(1:end-1)+edges(2:end));
%     betan(1) = mean(beta(1):beta(dsearchn(x',xmid(1))));
%     for k=2:length(edges)-1
%         Ix1 = dsearchn(x',xmid(k-1)); Ix2 = dsearchn(x',xmid(k));
%         betan(k) = nanmean(beta(Ix1:Ix2));
%     end
%     betan(length(edges)) = nanmean(beta(dsearchn(x',xmid(end-1)):end));
%     beta = interp1(edges,betan,x);
    
    beta(beta<0)=0; % beta cannot be less than 0
    beta(~isfinite(beta))=0; % cannot be infinite
    beta(gl+1:end) = 0; % zero roughness where ice is ungrounded

    % run flowline model for one year with resulting beta
    % define time stepping (s)
    dt = 0.01*3.1536e7;
    t_start = 0*3.1536e7;
    t_end = 8*3.1536e7;
    t = (t_start:dt:t_end);

    for i=1:length(t)

        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
        crev_s = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % surface crevasse penetration depth (m)
        Hab = H+rho_sw/rho_i*(hb); % height above buoyancy (m)
        crev_b = rho_i/(rho_sw-rho_i).*(Rxx./(rho_i*g)-Hab); % basal crevasse depth (m)
        % calving front located where the inland-most crevasse intersects sea level
        if i==1 % use observed calving front position for first iteration
            xcf = x0(c0);
        else 
            xcf_s = interp1(h-crev_s,x,0,'linear','extrap'); % (m along centerline)
            xcf_b = interp1(h-crev_b,x,0,'linear','extrap'); % (m along centerline)
                if xcf_s<0; xcf_s=NaN; end
                if xcf_b<0; xcf_b=NaN; end
            % calving front = whichever calving criteria occurs the
            % furthest inland
%             if xcf_s<xcf_b
%                 xcf = xcf_s;
%             else
%                 xcf = xcf_b;
%             end
            xcf = xcf_s;
        end

        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding line to the calving front
        if ~isempty(find(Hf-H>0,1,'first'))
            %xgl = x(find(Hf-H>0,1,'first')-1);
            xgl = interp1(Hf(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1)...
               -H(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),...
               x(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
        else
            xgl=xcf;
        end
        if xgl>xcf % grounding line can't be past calving front
            xgl=xcf;
        elseif xgl<0 % grounding line can't be negative
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
        A = interp1(x0,A0,xn,'linear','extrap');
        beta = interp1(x,beta,xn,'linear','extrap'); beta(beta<0)=0;
        x = xn; dx = dxn;

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

        % Solve for new velocity
        [U,dUdx,~,~,~,~,~] = U_convergence(x,U,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);

        % calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;

        % implement SMB & SMR
        smr = zeros(1,length(x));
        if gl<c
            smr(gl+1) = smr0;
            for k=gl+1:length(x)
                smr(k) = smr0-0.01*(smr0)*k;
            end
        end
        smb = interp1(x0,smb0+Q0,x);

        % calculate the  change in ice thickness from continuity
        clearvars dHdt
        dHdt(1) = (-1/W(1))*(F(1)-F(2))/(x(1)-x(2)); % forward difference
        dHdt(2:c-1) = (-1./W(2:c-1)).*(F(1:c-2)-F(3:c))./(x(1:c-2)-x(3:c)); % central difference
        dHdt(c:length(x)) = (-1./W(c:length(x))).*(F(c-1:length(x)-1)-F(c:length(x)))./(x(c-1:length(x)-1)-x(c:length(x))); % backward difference
        dH = dHdt.*dt;

        % new thickness (change from dynamics, SMB, & SMR)
        Hn = H+dH+(smb.*dt)+(smr.*dt);
        Hn(Hn < 0) = 0; % remove negative values
        H = Hn; % set as the new thickness value

        % stop the model if it behaves unstably (monitored by ice thickness and speed)
        if any(~isfinite(H(1:c))) || any(~isfinite(U(1:c))) || any(~isfinite(h(1:c)))
            disp('non finite values');
            break;
        end

    end

    % calculate cost of parameter solutions
    % modified from Morlighem et al., 2010; Larour et al., 2012; Kyrke-Smith et al., 2018
    U_err = 33/3.1536e7; % m/s
    h_err = 22; % m
    xcf_err = 15; % m
    K = log(beta); K(~isfinite(K))=NaN;
    J = 1.2*nanmean(sqrt((U-interp1(x0,U_2017,x)).^2)-U_err)/nanmean(interp1(x0,U_2017,x))+... % speed misfit term
       nanmean(sqrt((h-interp1(x0,h_2017,x)).^2)-h_err)./nanmean(interp1(x0,h_2017,x))+... % surface elevation misfit term
       nanmean(sqrt((xcf-xcf_2017).^2)-xcf_err)./xcf_2017+... % calving front location misfit term
       0.5*nanmean(abs(gradient(gradient(K)))); % regularization term to penalize changes in beta gradient
    
catch 

   J=NaN; beta = NaN; x=NaN; U=NaN;

end


end

