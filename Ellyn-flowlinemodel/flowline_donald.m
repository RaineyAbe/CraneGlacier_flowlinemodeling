function [t,XC,SMBsave,sol,x0,hb0] = flowline_donald

clear; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
plotflag = 1;

%% define time and space independent variables

% densities and g
p.rho_i = 917; % ice density (kg m^-3)
p.rho_sw = 1028; % ocean water density (kg m^-3)
p.rho_fw = 1000; % fresh water density (kg m^-3)
p.g = 9.81; % acceleration (m s^-2)

% stress parameters
p.m = 1/3; % basal sliding exponent
p.C = 5.624e6; % basal roughness factor ((m/s)^(-1/m)) 
p.n = 3; % flow law exponent
p.A = 4.22e-25; % rate factor (Pa^-n s^-1)

% numerical parameters
p.velreg = 1e-25; % velocity regulariser

% climate parameters
smb1 = 4/(86400*365); t1 = 200;
smb2 = 0.8/(86400*365); t2 = 1500;
smb3 = 4/(86400*365);

% grid spacing
dxmin = 200; % calving front grid spacing (m)
dxpad = 2000; % distance around these points with min grid spacing (m)
nramp = 5; % number of points over which to ramp up/down resolution
dxmax = 1000; % target max grid spacing (m) 

% time stepping (s)
dt = 0.1*86400*365;
t_start = 0*86400*365;
t_end = 2000*86400*365;
t = [t_start:dt:t_end];

% glacier geometry
x0 = [0:100:100e3]; % initial space vector
hb0 = 300-1500*x0/100e3; % bed topo (m)
W0 = 0*x0+5e3; % width (m)

% add gaussian bump
bumptop = 50e3;
bumpwidth = 5e3;
bumpheight = 100;
hb0 = hb0 + bumpheight*exp(-((x0-bumptop)/bumpwidth).^2);

% initialisation
H = 1000-1100*x0/100e3; H(H<0)=0; % ice thickness (m)
U = (100/(86400*365))*x0/100e3; % ice velocity (m/s)
x = x0; % initial space vector
hb = hb0; % bed topo (m)
W = W0; % width (m)

% calving
runtype = 'faf'; % fraction above flotation
faf = 1.05; % fraction above flotation

% saving
ts = [199,700,1499,1600];
si = 1;

%% Run the flowline model

if plotflag == 1,
    figure();
end

for i=1:length(t),
    
    % find the calving front location - subgridscale scheme
    if strmatch(runtype,'faf');
        xcf = interp1(-faf*p.rho_sw*hb/p.rho_i-H,x,0,'pchip');
    end
    XC(i) = xcf;
    
    % create coordinate system that hits cf exactly
    % has resolution dxmin around these points
    % has resolution dxmax away from these points
    % and has smooth variation between
    dxramp = linspace(dxmax,dxmin,nramp+2);
    ramplength = sum(dxramp(2:end-1));
    intlength = xcf-dxpad-ramplength;
    dxint = intlength/round(intlength/dxmax);
    dxn = [dxint*ones(1,round(intlength/dxmax)),dxramp(2:end-1),dxmin*ones(1,dxpad/dxmin)];
    xn = [0,cumsum(dxn)];
    
    % get geometry on new coordinates
    hb = interp1(x0,hb0,xn,'pchip');
    W = interp1(x0,W0,xn,'pchip');
    H = interp1(x,H,xn,'pchip');
    U = interp1(x,U,xn,'pchip');
    x = xn;
    
    % new surface mass balance
    if t(i)/(365*86400)<t1,
        SMB = smb1*ones(1,length(x));
        SMBsave(i) = smb1;
    elseif t(i)/(365*86400)<t2,
        SMB = smb2*ones(1,length(x));
        SMBsave(i) = smb2;
    else
        SMB = smb3*ones(1,length(x));
        SMBsave(i) = smb3;
    end

    % calculate surface elevation and slope
    h = hb+H; % surface elevation (m)
    
    if plotflag == 1,
        if mod(i-1,40)==0,
            subplot(1,2,1); cla;
            plot(x/10^3,h,'r','linewidth',1); hold on;
            plot(x(end)*[1,1]/10^3,[hb(end),h(end)],'r','linewidth',1);
            plot(x0/10^3,hb0,'k','linewidth',1);
            plot([x(end),x0(end)]/10^3,[0,0],'k--');
            title(['$t$ = ',num2str(t(i)/(365*86400)),' yrs']);
            xlim([0 70]); ylim([-600 1400]);
            xlabel('$x$ (km)'); ylabel('$h$ (m)');
            subplot(1,2,2); cla;
            plot(t(1:i)/(365*86400),XC(1:i)/10^3);
            xlim([0,t(end)]/(365*86400));
            xlabel('$t$ (yrs)'); ylabel('calving front position (km)');
            drawnow;
        end
    end
        
    % solve for velocity
    U = U_convergence(x,U,H,h,W,p);
    
    % save variables if wanted
    if ~isempty(find(t(i)/(365*86400)-ts==0)),
        sol(si).x = x;
        sol(si).U = U;
        sol(si).h = h;
        sol(si).t = t(i)/(365*86400);
        sol(si).hb = hb;
        sol(si).smb = SMBsave(i);
        sol(si).xc = x(end);
        si = si+1;
    end
       
    % calculate ice flux
    F = U.*H.*W; % ice flux (m^3 s^-1)
    
    % change in thickness term
    c = length(x);
    clearvars dHdt
    dHdt(1) = (-1/W(1))*(F(2)-F(1))/(x(2)-x(1)) + SMB(1);
    dHdt(2:c-1) = (-1./W(2:c-1)).*(F(3:c)-F(1:c-2))./(x(3:c)-x(1:c-2)) + SMB(2:c-1);
    dHdt(c) = (-1/W(c))*(F(c)-F(c-1))/(x(c)-x(c-1)) + SMB(c);
    
    % new thickness
    H = H + dHdt*dt; 
        
end

end

%% solve the stress balance to obtain speed
function U = U_convergence(x,U,H,h,W,p)

b = 1;
    
while b<=50,

    % get the H*visc term on the staggered grid
    % note Hnu(i) is the same as Hnu(i+1/2) in my notes
    c = length(x);
    Hnu(1:c-1) = 0.5*p.A^(-1/p.n)*(H(2:c)+H(1:c-1))...
                       .*(abs((U(2:c)-U(1:c-1))./(x(2:c)-x(1:c-1)))+p.velreg).^(1/p.n-1);

    % linearisation for lateral resistance
    gamma = U.^(1/p.n-1);
    % set maximum otherwise approaches infty as U->0
    Umin = 1e-9;
    gammamax = Umin^(1/p.n-1);
    gamma(gamma>gammamax) = gammamax;
    % set linearisation at ice divide
    % why? just because U=0 there? this would already be covered by gammamax    
    gamma(1) = gamma(2);

    % linearisation for basal resistance
    eta = U.^(p.m-1);
    % set maximum otherwise approaches infty as U->0
    etamax = Umin^(p.m-1);
    eta(eta>etamax) = etamax;
    % set linearisation at ice divide
    % why? just because U=0 there? this would already be covered by etamax
    eta(1) = eta(2);

    % set-up coefficient vectors for the linearized stress terms
    % as in Enderlin 2013: C(k)*U(k-1)+E(k)*U(k)+G(k)*U(k+1) = T(k)
    % set divide coefficients (imposes boundary condition U(1)=0)
    C(1) = 0;
    G(1) = 0;
    T(1) = 0;
    E(1) = 1;
    % coefficients up to calving front
    C(2:c-1) = 2*(2./(x(3:c)-x(1:c-2))).*Hnu(1:c-2)./(x(2:c-1)-x(1:c-2));
    G(2:c-1) = 2*(2./(x(3:c)-x(1:c-2))).*Hnu(2:c-1)./(x(3:c)-x(2:c-1));
    E(2:c-1) = 2*(-2./(x(3:c)-x(1:c-2))).*(Hnu(2:c-1)./(x(3:c)-x(2:c-1))+Hnu(1:c-2)./(x(2:c-1)-x(1:c-2)))...
               -p.C*eta(2:c-1)-((gamma(2:c-1).*H(2:c-1))./W(2:c-1)).*((5./(2*p.A*W(2:c-1))).^(1/p.n));
    T(2:c-1) = p.rho_i*p.g*H(2:c-1).*(h(3:c)-h(1:c-2))./(x(3:c)-x(1:c-2));        
    % at calving front
    C(c) = -1;
    E(c) = 1;
    G(c) = 0;
    T(c) = (p.A*(((p.rho_i*p.g/4)*H(c)*(1-(p.rho_i/p.rho_sw))))^p.n)*(x(c)-x(c-1));

    % create a sparse tri-diagonal matrix for the velocity coefficient vectors
    M = diag(C(2:c),-1) + diag(E(1:c)) + diag(G(1:c-1),1);
    M = sparse(M);

    % perform the matrix inversion to solve for new ice velocities
    Un = (M\T')';

    % check for convergence
    if abs(sum(U-Un))<0.1*abs(sum(U)),
        % if converged, output new velocities
        U = Un; 
        return % break the U iterations
    else
        % if not sufficiently converged, set Un to U and solve the stress balance again
        U = Un;
        b = b+1;      
    end

end

end

