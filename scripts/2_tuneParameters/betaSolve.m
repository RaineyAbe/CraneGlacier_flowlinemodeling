function [J,U,x,xcf,beta0x] = betaSolve(A0,A,beta0,H,x,U,hb,n,E,m,dx0,rho_i,g,rho_sw,rho_fw,d_fw,sigma_b,dUdx,c0,x0,hb0,W0,U_2018,xcf_2018,SMB0,Q0,SMR0,H_max,U_min,F0)

% interpolate beta0 to grid spacing
beta0x = 0:round(x0(end)/(length(beta0))):round(x0(end)/(length(beta0)))*(length(beta0)-1); 
beta = interp1(beta0x,beta0,x0,'pchip');
    
% define time stepping (s)
dt = 0.001*3.1536e7;
t_start = 0*3.1536e7;
t_end = 9*3.1536e7;
t = (t_start:dt:t_end);

%try
    for i=1:length(t)
        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
        crev_s = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*d_fw); % surface crevasse penetration depth (m)
        Hab = H+rho_sw/rho_i*(hb); % height above buoyancy (m)
        crev_b = rho_i/(rho_sw-rho_i).*(Rxx./(rho_i*g)-Hab); % basal crevasse depth (m)
        % calving front located where the inland-most crevasse intersects sea level
        if i==1 % use observed calving front position for first iteration
            xcf = x0(c0);
        else
            if length(h)>=find(h-crev_s<0,1,'first')+1
                xcf_s = interp1(h(find(h-crev_s<0,1,'first')-1:find(h-crev_s<0,1,'first')+1)...
                    -crev_s(find(h-crev_s<0,1,'first')-1:find(h-crev_s<0,1,'first')+1),...
                    x(find(h-crev_s<0,1,'first')-1:find(h-crev_s<0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
            else
                xcf_s = interp1(h-crev_s,x,0,'linear','extrap');
            end
            if length(h)>=find(h-crev_b<0,1,'first')+1
                xcf_b = interp1(h(find(h-crev_b<0,1,'first')-1:find(h-crev_b<0,1,'first')+1)...
                    -crev_b(find(h-crev_b<0,1,'first')-1:find(h-crev_b<0,1,'first')+1),...
                    x(find(h-crev_b<0,1,'first')-1:find(h-crev_b<0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
            else
                xcf_b = interp1(h-crev_b,x,0,'linear','extrap');
            end
            if xcf_s<0; xcf_s=NaN; end
            if xcf_b<0; xcf_b=NaN; end
            % calving front = whichever calving criteria occurs the
            % furthest inland
            if xcf_s<xcf_b
                xcf = xcf_s;
            else
                xcf = xcf_b;
            end
            xcf=xcf_s; % use only surface crevasses if basal crevasses
            %lead to unstable model behavior
            if xcf<20e3 || xcf > 70e3 || isnan(xcf)
                xcf = interp1(feval(fit(x',(h-crev_s)','poly1'),x),x,0,'linear','extrap');
                %xcf = x0(c0);
            end
        end
        
        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding linU_widthavge to the calving front
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
        XGL(i) = xgl; % save grounding line position over time
        
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
        if isempty(c) == 1 % set the calving front to a default minimum ice thickness value
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
        beta = interp1(x,beta,xn,'linear','extrap');
        x = xn; dx = dxn; clear xn dxn; %EE: added the clear statement 14/09/21
        
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
        [U,dUdx,~,~,~,~,~,~,~,~,~,~,~] = U_convergence(x,U,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);
        %EE: removed U0 and dUdx from inputs to U_convergence 15/09/21
        
        % calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        %F(1)=F(2)+F0;
        
        % implement SMB, SMR, delta_SMB, & delta_SMR
        smr = zeros(1,c);
        for k=gl+1:c
            smr(k) = SMR0-0.001*(SMR0)*(k-gl+1);
        end
        smb = interp1(x0,SMB0+Q0,x);
            
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
        if max(H) > H_max
            disp(['Adjust dt']);
            break;
        end
        if mean(U) < U_min
            disp('Too slow!');
            break;
        end
        if any(~isfinite(H(1:c))) || any(~isfinite(U(1:c))) || any(~isfinite(h(1:c)))
            disp('non finite values');
            break;
        end
        
    end
    
    % calculate cost of parameter solutions
    % modified from Morlighem et al., 2010; Larour et al., 2012; Kyrke-Smith et al., 2018
    U_err = 33/3.1536e7; % m/s
    h_err = 22; % m
    K = log(beta); K(~isfinite(K))=0;
    J = nanmean(sqrt((U-interp1(x0(1:dsearchn(x0',xcf_2018)),U_2018(1:dsearchn(x0',xcf_2018)),x)).^2)-U_err)/nanmean(U_2018(1:dsearchn(x0',xcf_2018)))+... % speed misfit term
        nanmean(abs(gradient(gradient(K)))); % regularization term to penalize changes in beta gradient        
        %nanmean(sqrt((h-interp1(x0,h_2018,x)).^2)-h_err)./nanmean(h_2018); % surface elevation misfit term

%catch
%    J=NaN;
%end

end