%% Script to tune the basal roughness factor beta, 
% using observations of glacier speed and geometry.
% Adapted from Ellyn Enderlin's flowline model (Enderlin et al., 2013)
% Rainey Aberle
% 2020-2021

%% define constants and spatial variables
  
clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean/']);
cd([homepath,'inputs-outputs/']);

save_params = 1;  % = 1 to save resulting optimal parameters

% Load Crane Glacier initialization variables
load('Crane_flowlineModelInitialization.mat');
cl.X = load('Crane_centerline.mat').x; cl.Y = load('Crane_centerline.mat').y;
%SCL = load('Crane_SCL_results.mat').

% Load 2017 ice surface speed observation
U_2017 = double(load('Crane_centerlineSpeeds_2007-2017.mat').U(19).speed');
h_2017 = double(load('Crane_surfaceElevationObs.mat').h(35).surface);
h_2017(isnan(h_2017))=0;
xcf_2017 = x0(dsearchn([cl.X cl.Y],[load('LarsenB_centerline.mat').centerline.x(32),load('LarsenB_centerline.mat').centerline.y(32)]));
        
% grid spacing
dx0 = 200; % m

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

%% solve for optimal model parameters 

% define method
method = 'grad'; % = 'brute' to solve by brute force
                 % = 'grad' to solve by gradient descent

% define parameter bounds
Ebounds = [1 5];
fwdbounds = [0 20];

if strcmp(method,'grad')
    % solve by gradient descent
    % start timer
    tic
    % initialize variables
    x=x0; H=H0; hb=hb0; U=U0; dUdx=dUdx0;
    % solve for parameters which minimize the cost function for H, U, and xcf
    cd([homepath,'scripts/']);
    % Create function handle for the beta solve function
    fh = @(b)betaSolve(H,x,U,hb,n,b(1),b(2),m,dx0,rho_i,g,rho_sw,rho_fw,sigma_b,dUdx,c0,x0,hb0,W0,A0,U_2017,h_2017,xcf_2017,smr0,smb0,Q0);
    b0 = [3 15]; % initial guesses for E and fwd
    % Solve for the parameter values which minimize the cost function J
    [b,fval] = fmincon(fh,b0,[],[],[],[],[Ebounds(1) fwdbounds(1)],[Ebounds(2) fwdbounds(2)]);
    % initialize variables
    x=x0; H=H0; hb=hb0; U=U0; A=A0; dx=dx0; dUdx=dUdx0;
    % implement and grab results
    E = b(1); fwd = b(2);
    [J,beta,x,U,H,h] = betaSolve(H,x,U,hb,n,E,fwd,m,dx0,rho_i,g,rho_sw,rho_fw,sigma_b,dUdx,c0,x0,hb0,W0,A0,U_2017,h_2017,xcf_2017,smr0,smb0,Q0);
    % stop timer
    toc
elseif strcmp(method,'brute')
    % solve by brute force and plot cost
    % start timer
    tic
    % initialize variables
    x=x0; H=H0; hb=hb0; U=U0; dUdx=dUdx0;
    E = Ebounds(1):Ebounds(2); fwd = fwdbounds(1):0.5:fwdbounds(2);
    J = zeros(length(E),length(fwd));
    cd([homepath,'scripts/']);
    for i=1:length(E)
        for j=1:length(fwd)
            [J(i,j),~,~,~,~,~] = betaSolve(H,x,U,hb,n,E(i),fwd(j),m,dx0,rho_i,g,rho_sw,rho_fw,sigma_b,dUdx,c0,x0,hb0,W0,A0,U_2017,h_2017,xcf_2017,smr0,smb0,Q0); 
        end
    end
    % plot results for cost
    figure(3); clf; hold on; grid on;
    set(gca,'fontsize',16,'fontname','arial','linewidth',2);
    colormap(cmocean('thermal'));
    imagesc(E,fwd,J');
    xlabel('E'); ylabel('fwd (m)'); xlim([min(E) max(E)]); ylim([min(fwd) max(fwd)]);
    c = colorbar; caxis([0 1]); c.Label.String = 'Cost';
    [IEbest,Ifwdbest] = find(J==min(J(:)),1);
    E = E(IEbest); fwd = fwd(Ifwdbest);
    plot(E,fwd,'xw','markersize',10,'linewidth',2);
    % grab values using best parameters
    [Jbest,beta,x,U,H,h] = betaSolve(H,x,U,hb,n,E,fwd,m,dx0,rho_i,g,rho_sw,rho_fw,sigma_b,dUdx,c0,x0,hb0,W0,A0,U_2017,h_2017,xcf_2017,smr0,smb0,Q0);
    % stop timer
    toc
end

% plot results
figure(4); clf
set(gcf,'Position',[441   145   900   800]);
subplot(3,1,1:2);
    hold on; grid on; set(gca,'fontsize',16,'fontname','arial','linewidth',2);
    ylabel('\beta (s^{1/m} m^{-1/m})'); legend('Location','east');
    ylim([-2 max(beta)+2]); title(['E=',num2str(E),', fwd=',num2str(fwd),' m']);
    plot(x./10^3,beta,'linewidth',2,'displayname','\beta');
    yyaxis right; ylabel('U m/yr'); ylim([-200 1400]);
    plot(x./10^3,U.*3.1536e7,'k','linewidth',2,'displayname','U_{sol}'); 
    plot(x0./10^3,U_2017.*3.1536e7,'--k','linewidth',2,'displayname','U_{obs}');
subplot(3,1,3);
    set(gca,'fontsize',16,'fontname','arial','linewidth',2);
    xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
    hold on; grid on; legend('Location','east');
    plot(x/10^3,h,'k','linewidth',2,'displayname','h_{sol}');
    plot(x0/10^3,h_2017,'--k','linewidth',2,'displayname','h_{obs}');
    
% save resulting optimal parameters
if save_params
    cd([homepath,'inputs-outputs/']);
    save('optimalParameters.mat','E','fwd','beta','x');
    disp('optimal parameters saved.');
end