%% Script to tune the basal roughness factor beta, 
% using observations of glacier speed and geometry.
% Adapted from Ellyn Enderlin's flowline model (Enderlin et al., 2013)
% Rainey Aberle
% 2020-2021

%% 0. define time and space independent variables
  
clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';
cd([homepath,'inputs-outputs/']);

save_betabest = 1;  % = 1 to save the beta with the best U RMSE
save_results = 1;   % = 1 to save all solution results 

% grid spacing
dx0 = 100:50:3000;

% Load Crane Glacier initialization variables
load('Crane_flowlineModelInitialization.mat');
        
% densities and g
rho_i = 917; % ice density (kg m^-3)
rho_sw = 1028; % ocean water density (kg m^-3)
rho_fw = 1000; % fresh water density (kg m^-3)
g = 9.81; % acceleration (m s^-2)

% stress parameters (unitless)
m = 3; % basal sliding exponent
n = 3; % flow law exponent

% calving parameters
Hc = 400; % m -> set the calving front to a default minimum ice thickness value if calving criteria is not met
fwd = 30; % fresh water depth in crevasses (m)
sigma_b = 0; % back pressure (Pa) due to sea ice or sikkusak

% maximum & minimum thickness & speed cut-off to check for instability
H_max = 2000; % maximum thickness (m)
H_min = 100;  % minimum thickness (m)
U_min = 200;  % minimum mean speed (m a^-1)

% initial conditions
H0 = h0-hb0; % ice thickness (m)
dUdx0 = [(U0(1:end-1)-U0(2:end))./(x0(1:end-1)-x0(2:end)) 0]; % strain rate (1/s)
% find the location of the grounding line and the end of the ice-covered domain
Hf = -(rho_sw./rho_i).*hb0; % flotation thickness (m)
gl0 = find(Hf-H0>0,1,'first')-1; % grounding line location 
H0(gl0:length(x0))=h0(gl0:end)*rho_sw/(rho_sw-rho_i); % buoyant thickness using surface
H0(H0>=(h0-hb0))=h0(H0>=(h0-hb0))-hb0(H0>=(h0-hb0)); % thickness can't go beneath bed elevation
H0(c0+1:end) = 0; % zero ice thickness past calving front

% extend variables past where there are NaNs (set to last value)
A0(find(isnan(A0),1,'first'):end) = A0(find(isnan(A0),1,'first')-1);
U0(find(isnan(U0),1,'first'):end) = U0(find(isnan(U0),1,'first')-1);
hb0(find(isnan(hb0),1,'first'):end) = hb0(find(isnan(hb0),1,'first')-1);
W0(find(isnan(W0),1,'first'):end) = W0(find(isnan(W0),1,'first')-1);    
E0 = 1*ones(1,length(x0)); % enhancement factor

%% 1. run the flowline model

clear beta

% loop through grid spatial resolutions
for j=1:length(dx0)
    
    % set dx
    beta(j).dx = dx0(j);
    
    % initialize variables
    x=x0; H=H0; hb=hb0; U=U0; A=A0; E=E0; dUdx=dUdx0;
    
    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
    % calving front located where the inland-most crevasse intersects sea level
    xcf = x(c0);
    
    % calculate the thickness required to remain grounded at each grid cell
    Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
    % find the location of the grounding line and use a floating
    % geometry from the grounding line to the calving front
    if ~isempty(find(Hf-H>0,1,'first'))
        xgl = interp1(Hf(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1)...
            -H(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),...
            x(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
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
    xl = round(xgl/dx0(j)); %number of ideal grid spaces needed to reach the grounding line
    dx = xgl/xl; %new grid spacing (should be ~dx0)
    xn = 0:dx:xgl; %new distance vector
    if xcf-xgl > 0
        xl = round((xcf-xgl)/dx0(j));
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
    A = interp1(x0,A0,xn,'pchip');
    E = interp1(x,E,xn,'pchip');
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
    
    % solve for new beta, save results
    [beta(j).beta,dUdx,Un_rev] = betaSolve(H,c,x,U,n,A,E,m,dx,rho_i,g,h,rho_sw,sigma_b,W,N);
    beta(j).x = x; beta(j).U = Un_rev;
    
    % plot geometry, speed, & calving front position
    if j==1
        col = parula(length(dx0)+2); % color scheme for plots
        figure(1); clf
        set(gcf,'Position',[0 100 1300 400]);
        ax1 = axes('Position',[0.05 0.1 0.25 0.8]); % glacier geometry
        hold on; grid on; set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
        title('Glacier Geometry'); xlim([0 50]); ylim([min(hb)-100 max(h)+200]);
        xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
        % ice surface
        plot(x(1:c)./10^3,h(1:c),'color',col(j,:),'linewidth',2,'displayname',num2str(dx0(j)));
        % calving front
        plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(j,:),'linewidth',2,'HandleVisibility','off');
        % floating bed
        plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(j,:),'linewidth',2,'HandleVisibility','off');
        % bed elevation
        plot(x0./10^3,hb0,'k','linewidth',2,'HandleVisibility','off');
        % mean sea level
        plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
        ax2 = axes('Position',[0.35 0.1 0.25 0.8]); % ice speed
        hold on; grid on; set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
        title('Ice Speed'); xlim([0 50]); ylim([0 2500]);
        xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
        plot(x(1:c)./10^3,Un_rev(1:c).*3.1536e7,'color',col(j,:),'linewidth',2);
        ax3 = axes('Position',[0.7 0.1 0.25 0.8]); % beta
        hold on; grid on; set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
        title('Basal Roughness Factor \beta'); xlim([0 50]);
        xlabel('Distance Along Centerline (km)'); ylabel('\beta (s^{1/m} m^{-1/m})');
        plot(x(1:c-1)./10^3,beta(j).beta(1:c-1),'color',col(j,:),'linewidth',2);
    else
        % ice surface
        plot(ax1,x(1:c)./10^3,h(1:c),'color',col(j,:),'linewidth',2);
        % calving front
        plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(j,:),'linewidth',2);
        % floating bed
        plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(j,:),'linewidth',2);
        % ice speed
        plot(ax2,x(1:c)./10^3,Un_rev(1:c).*3.1536e7,'color',col(j,:),'linewidth',2);
        % beta
        plot(ax3,x(1:c-1)./10^3,beta(j).beta(1:c-1),'color',col(j,:),'linewidth',2);
    end
    if j==length(dx0)
        figure(1);
        % colorbar
        colormap(parula(length(col))); cb=colorbar('Ticks',[0 0.5 1],...
            'YTickLabel',[{num2str(dx0(1)/10^3)} {num2str(dx0(round(length(dx0)/2))/10^3)} {num2str(dx0(end)/10^3)}],...
            'Position',[.96 .35 .025 .3410],'Fontname','arial','fontweight','bold');
        cb.Title.String = 'dx (km)';
    end
    
    % Calculate RMSE of resulting speed using beta solution
    beta(j).RMSE = sqrt((sum((Un_rev-U).^2)./length(U)));
    RMSE(j) = sqrt((sum((Un_rev-U).^2)./length(U)));
    
    % display results
    disp(['dx = ',num2str(dx0(j)),' m: ']);
    disp(['RMSE = ',num2str(RMSE(j)),' m/s']);
    
end

% define best spatial resolution as the minimum spatial resolution where
% the gradient is less than the velocity spatial resolution
for i=1:length(beta)
    if all(abs(gradient(RMSE(i:end)*3.1536e7))<=240)
        Ibest=i;
        break;
    end
end
if length(Ibest)>1 % if lowest RMSE exists for multiple dx, use lowest dx
    Ibest(2:end)=[];
end
figure(4); clf
set(gcf,'Position',[441   145   894   652]);
subplot(2,2,1:2);
    set(gca,'fontsize',14,'fontname','arial','linewidth',2);
    xlabel('Distance Alonc Centerline (km)'); ylabel('\beta (s^{1/m} m^{-1/m})'); 
    title('\beta_{best}'); hold on; grid on; 
    plot(beta(Ibest).x./10^3,beta(Ibest).beta,'linewidth',2);  
subplot(2,2,3);
    set(gca,'fontsize',14,'fontname','arial','linewidth',2);
    xlabel('dx (km)'); ylabel('RMSE (m a^{-1})'); title('RMSE of \beta Solutions');
    hold on; plot(dx0./10^3,RMSE.*3.1536e7,'.','markersize',20); grid on;
    plot(dx0(Ibest)/10^3,RMSE(Ibest)*3.1536e7,'*','markersize',20,'linewidth',2);  
subplot(2,2,4);
    set(gca,'fontsize',14,'fontname','arial','linewidth',2);
    xlabel('x (m along centerline)'); ylabel('U (m a^{-1})'); title('U Solution');
    hold on; grid on; legend('Location','best');
    plot(x0(1:c0)./10^3,U0(1:c0).*3.1536e7,'k','displayname','observed','linewidth',2);
    plot(beta(Ibest).x./10^3,beta(Ibest).U.*3.1536e7,'--k','displayname','solution','linewidth',2);    

% save beta with lowest U RMSE 
if save_betabest
    cd([homepath,'inputs-outputs/']);
    optBeta.beta = beta(Ibest).beta; 
    optBeta.x = beta(Ibest).x;
    optBeta.dx = beta(Ibest).dx;
    save('optimalBeta.mat','optBeta');
    % save to initialization file
    beta0 = interp1(beta(Ibest).x,beta(Ibest).beta,x0);
    save('Crane_flowlineModelInitialization.mat','beta0','-append');
    disp('optimalBeta saved.');
end

% save solution results
if save_results
    cd([homepath,'inputs-outputs/']);
    betaSolutions = beta; 
    save('betaSolutions.mat','betaSolutions');
    disp('betaSolutions saved.');
end

%% solve the stress balance equations to obtain the basal roughness factor (beta)
function [beta,dUdx,Un_rev] = betaSolve(H,c,x,U,n,A,E,m,dx,rho_i,g,h,rho_sw,sigma_b,W,N)
        
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
    T(c) = (E(c)*A(c).*(((rho_i.*g./4).*((H(c).*(1-(rho_i./rho_sw))-sigma_b./(rho_i.*g)))).^n)).*dx(c); %SOMETHING WRONG HERE
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
    beta(1) = beta(2); beta(c)=0;
    beta(beta<0)=0; % beta cannot be less than 0
    
    % Run the forward U_convergence with the resulting beta to check success
        %set-up coefficient vectors for the linearized stress terms over the calving front
        %[C(k)*U(k-1)+E(k)*U(k)+G(k)*U(k+1)=Td]  
        % upper boundary condition
        G_minus_rev(1) = 0;
        G_rev(1) = -(beta(1).*N(1).*eta(1))-...
            (((2*gamma(1).*H(1))./W(1)).*((5/(E(1)*A(1).*W(1))).^(1/n)));
        G_plus_rev(1) = 0;
        T_rev(1) = (rho_i.*g.*H(1).*(h(1)-h(2))./(x(1)-x(2)));    
        % coefficients up to calving front
        G_minus_rev(2:c-1) = (2./(dx(2:c-1).^2)).*Hm(1:c-2).*vm(1:c-2); %for U(k-1)
        G_rev(2:c-1) = (-2./(dx(2:c-1).^2)).*(Hm(1:c-2).*vm(1:c-2)+Hm(1:c-2).*vm(1:c-2))-...
            (beta(2:c-1).*(N(2:c-1)).*eta(2:c-1))-...
            (((2*gamma(2:c-1).*H(2:c-1))./W(2:c-1)).*((5./(E(2:c-1).*A(2:c-1).*W(2:c-1))).^(1/n))); %for U(k)        
        G_plus_rev(2:c-1) = (2./(dx(2:c-1).^2)).*Hm(2:c-1).*vm(2:c-1); %for U(k+1)
        T_rev(2:c-1) = (rho_i.*g.*H(2:c-1).*(h(3:c)-h(1:c-2))./(x(3:c)-x(1:c-2))); %gravitational driving stress  
        % calving front condition
        G_minus_rev(c) = -1;
        G_rev(c) = 1;
        G_plus_rev(c) = 0;
        T_rev(c) = (E(c)*A(c).*(((rho_i.*g./4).*((H(c).*(1-(rho_i./rho_sw))-sigma_b./(rho_i.*g)))).^n)).*dx(c);
        %remove any NaNs from the coefficient vectors
        G_minus_rev(isnan(G_minus_rev)) = 0;
        G_rev(isnan(G_rev)) = 0;
        G_plus_rev(isnan(G_plus_rev)) = 0;
        T_rev(isnan(T_rev)) = 0;

        %create a sparse tri-diagonal matrix for the velocity coefficient vectors
        M = diag(G_minus_rev(2:c),-1) + diag(G_rev(1:c)) + diag(G_plus_rev(1:c-1),1);
        M = sparse(M);
            
        %make sure Td is a column vector for the inversion 
        if size(T_rev)==[1,c]
            T_rev=T_rev';
        end

        %use the backslash operator to perform the matrix inversion to solve for ice velocities
        Un_rev = M\T_rev; %velocity (m s^-1)
        %remove NaNs
        Un_rev(isnan(Un_rev)) = 0;
        Un_rev(Un_rev<0)=0; % U cannot be less than zero
        
        %make sure Un is a row vector so it can be compared with U
        if size(Un_rev) == [c,1]
            Un_rev=Un_rev';
        end
        
        %calculate new strain rates (equal to the velocity gradient)
        dUndx_rev(1) = (Un_rev(2)-Un_rev(1))/(x(2)-x(1)); % backward difference
        dUndx_rev(2:c-1) = (Un_rev(3:c)-Un_rev(1:c-2))./(x(3:c)-x(1:c-2)); % central difference
        dUndx_rev(c) = (Un_rev(c)-Un_rev(c-1))./(x(c)-x(c-1)); % forward differe ce

        %make sure dUdx is a row vector
        if size(dUndx_rev) == [c,1]
            dUndx_rev=dUndx_rev';
        end

        %check if the difference in speed between iteratons (U vs. Un_rev) meets a set tolerance   
        if abs(sum(U)-sum(Un_rev))<0.1*abs(sum(U))==0 %determine if U has converged sufficiently
            %use sufficiently converged values for speeds & strain rates
           disp('U did not converge sufficiently.');
        end
               
            
end 

