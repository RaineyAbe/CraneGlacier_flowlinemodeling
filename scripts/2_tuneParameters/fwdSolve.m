function [J,x,U,xcf] = fwdSolve(dUdx,U,E,A,m,n,rho_i,g,rho_fw,rho_sw,fwd,x0,c0,h,hb,x,H,dx0,hb0,W0,A0,beta0,plot_timeSteps,smr0,smb0,Q0,sigma_b,H_max,U_min,xcf_2018,U_2018)

% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 9*3.1536e7;
t = (t_start:dt:t_end);

try
    
    % run flowline model
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
            xcf=xcf_s;
            if xcf<20e3 || xcf > 70e3 || isnan(xcf)
                xcf = interp1(feval(fit(x',(h-crev_s)','poly1'),x),x,0,'linear','extrap');
                %xcf = x0(c0);
            end
        end

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
        beta = interp1(x0,beta0,xn,'linear','extrap'); beta(gl+1:end)=0;
        x = xn; dx = dxn;

        % calculate surface elevation
        h = hb+H; % surface elevation (m a.s.l.)
        h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); %adjust the surface elevation of ungrounded ice to account for buoyancy
        H(h<0)=0-hb(h<0); h(h<0)=0; % surface cannot go below sea level
        h(h-H<hb) = hb(h-H<hb)+H(h-H<hb); % thickness cannot go beneath bed elevation

        % plot geometry, speed, & grounding line and calving front positions
        col = parula(length(t)+20); % color scheme for plots
        if plot_timeSteps
            if t(i)==t_start
                figure(1); clf
                set(gcf,'Position',[0 100 1300 400]);
                ax1 = axes('Position',[0.06 0.12 0.27 0.78]); % glacier geometry
                hold on; grid on;
                set(gca,'FontSize',14,'linewidth',2);
                %title('Glacier Geometry');
                legend('Location','northeast');
                xlim([0 75]); ylim([min(hb)-100 max(h)+200]);
                xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
                % ice surface
                plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','2009');
                % calving front
                plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % floating bed
                plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % bed elevation
                plot(x0./10^3,hb0,'k','linewidth',2,'HandleVisibility','off');
                % mean sea level
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
                ax2 = axes('Position',[0.38 0.1 0.28 0.8]); % ice speed
                hold on; grid on;
                set(gca,'FontSize',14,'linewidth',2);
                title('Ice Speed'); legend('Location','northeast');
                xlim([0 65]); ylim([0 1500]);
                xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
                % ice speed
                plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','2009');
                ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
                hold on; grid on;
                set(gca,'FontSize',14,'linewidth',2);
                title('Terminus & Grounding Line Positions'); legend('Location','best');
                xlim([30 65]); ylim([0 t_end./3.1536e7]);
                xlabel('Distance Along Centerline (km)'); ylabel('Year');
                % terminus & grounding line positions
                plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','2009');
                plot(ax3,x(gl)./10^3,t(i)./3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off');
            elseif mod(i-1,round(length(t)/10))==0 % display every length(t)/20
                figure(1);
                % ice surface
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(round(t(i)./3.1536e7)+2009));
                plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                plot(ax1,[x(c);x(c)]/10^3,[h(c);h(c)-H(c)],'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % calving front
                plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % floating bed (gl:c)
                plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(round(t(i)./3.1536e7)+2009));
                % calving front position
                plot(ax3,x(c)/10^3,t(i)/3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(round(t(i)./3.1536e7)+2009)); hold on;
                % grounding line position
                plot(ax3,x(gl)/10^3,t(i)/3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off'); hold on;
            end
        end
        
        % calculate the effective pressure (ice overburden pressure minus water
        % pressure) assuming an easy & open connection between the ocean and
        % ice-bed interface
        sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
        N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
        N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
        N = [N_ground N_marine];
        N(N<0)=0; % cannot have negative values

        % Solve for new velocity
        [U,dUdx,~,~,~] = U_convergence(x,U,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);

        % calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        F(1)=F(2);

        % implement SMB, SMR, delta_SMB, & delta_SMR
        smr = zeros(1,c);
        for k=gl+1:c
            smr(k) = smr0-0.001*(smr0)*(k-gl+1);
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
        if max(H) > H_max
            disp([num2str(fwd),' m failed']);
            J = NaN;           
            break;
        end
        if mean(U) < U_min/3.1536e7
            disp([num2str(fwd),' m failed']);
            J = NaN;
            break;
        end
        if any(~isfinite(H(1:c))) || any(~isfinite(U(1:c))) || any(~isfinite(h(1:c)))
            disp([num2str(fwd),' m failed']);
            J = NaN;  
            break;
        end

        if t(i)==t_end
            % calculate cost of resulting terminus position
            % modified from Morlighem et al., 2010; Larour et al., 2012; Kyrke-Smith et al., 2018
            U_err = 33/3.1536e7; % m/s
            h_err = 22; % m
            xcf_err = 15; % m
            J = abs((xcf-xcf_2018-xcf_err)/xcf_2018);%+... % surface elevation misfit term
            %nanmean(sqrt((U-interp1(x0,U_2018,x)).^2)-U_err)/nanmean(U_2018);%+... % speed misfit term
        end
        
    end

    
catch
    J=NaN; x=NaN; U=NaN; xcf=NaN;
    disp([num2str(fwd),' m failed']);
end

end