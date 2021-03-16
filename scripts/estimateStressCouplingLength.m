%% Estimate Stress-Coupling Length
% Rainey Aberle 2021
% Script to estimate the glacier stress-coupling length using methods 
% adapted from Enderlin et al. (2016) at Crane Glacier, Antarctic Peninsula. 
% 1. Load centerline observations of ice speed, ice thickness, glacier width
% 2. 

clear all; close all;

save_results = 1; % = 1 to save results

homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/';
cd([homepath,'inputs-outputs/']);

% 1. Load centerline observations of ice speed, ice thickness, glacier width
% and define constants

    x = load('Crane_flowline_initialization.mat').x0; % spatial grid (m along centerline)
    hb = load('Crane_flowline_initialization.mat').hb0; % glacier bed elevation (m)
    h = load('Crane_flowline_initialization.mat').h0; % ice surface elevation (m)
    W = load('Crane_flowline_initialization.mat').W0; % glacier width (m)
    U = load('Crane_flowline_initialization.mat').U0; % ice surface speed (m/s)
    A = load('Crane_flowline_initialization.mat').A0; % rate factor (Pa^-n s^-1)
    c = load('Crane_flowline_initialization.mat').c0; % calving front position
    x(c+1:end)=[]; hb(c+1:end)=[]; h(c+1:end)=[]; W(c+1:end)=[];
    U(c+1:end)=[]; A(c+1:end)=[];
    H = h-hb; % ice thickness (m)
    
    % Define constants
    n = 3; % flow law exponent (unitless)
    m = 1; % basal sliding exponent (unitless)
    rho_i = 917; % density of ice (kg/m^3)
    rho_sw = 1028; % ocean water density (kg m^-3)
    g = 9.81; % gravitational acceleration (m/s^2)
    E = 1; % enhancement factor (unitless)
    beta = 2; % basal roughness factor (s^(1/m) m^(-1/m))
    
    % Define averaging window
    w = [30 105]; % indices of spatial grid to use in averaging

% 2. Calculate gradients in speed and ice surface elevation 

    % fit a linear regression to speed observations
    U_lin = feval(fit(x(w(1):w(2))',U(w(1):w(2))','poly1'),x)'; % (m/s)
    
    % calculate strain rates dUdx (1/s)
    dUdx = zeros(1,length(x));
    dUdx(1) = (U(2)-U(1))./(x(2)-x(1)); % forward difference
    dUdx(2:end-1) = (U(3:end)-U(1:end-2))./(x(3:end)-x(1:end-2)); % central difference
    dUdx(end) = (U(end)-U(end-1))./(x(end)-x(end-1)); % backward difference
    % use a linear regression for dUdx
    dUdx_lin = feval(fit(x(w(1):w(2))',dUdx(w(1):w(2))','poly1'),x)'; 
    
    % fit a linear regression to ice surface elevation observations
    h_lin = feval(fit(x(w(1):w(2))',h(w(1):w(2))','poly1'),x)'; % (m/m)
    
    % calculate surface slope (unitless)
    dhdx = zeros(1,length(x));
    dhdx(1) = (h(2)-h(1))./(x(2)-x(1)); % forward difference
    dhdx(2:end-1) = (h(3:end)-h(1:end-2))./(x(3:end)-x(1:end-2)); % central difference
    dhdx(end) = (h(end)-h(end-1))./(x(end)-x(end-1)); % backward difference
    % use a linear regression for dhdx
    dhdx_lin = feval(fit(x(w(1):w(2))',dhdx(w(1):w(2))','poly1'),x)';   
    
    % plot speed, strain rates, surface, & surface slope
    figure(1); clf
    set(gcf,'position',[10 245 749 405]);    
    subplot(1,2,1);
        set(gca,'fontsize',12,'fontname','arial','linewidth',2);
        hold on; grid on; legend('Location','best');
        plot(x,U,'.k','markersize',10,'displayname','U');
        plot(x,U_lin,'--k','linewidth',2,'displayname','U_{lin}');
        ylabel('m a^{-1}'); title('Ice Surface Speed');
        yyaxis right;
        plot(x,dUdx,'.r','markersize',10,'displayname','\partialU/\partialx');
        plot(x,dUdx_lin,'--r','linewidth',2,'displayname','\partialU/\partialx_{lin}');
        xlabel('distance along centerline (km)'); ylabel('s^{-1}');
    subplot(1,2,2);
        set(gca,'fontsize',12,'fontname','arial','linewidth',2);
        hold on; grid on; legend('Location','best');
        plot(x,h,'.k','markersize',10,'displayname','h');
        plot(x,h_lin,'--k','linewidth',2,'displayname','h_{lin}');
        ylabel('m.a.s.l.'); title('Ice Surface Elevation');
        yyaxis right;
        plot(x,dhdx,'.r','markersize',10,'displayname','\partialh/\partialx');
        plot(x,dhdx_lin,'--r','linewidth',2,'displayname','\partialh/\partialx_{lin}');
        xlabel('distance along centerline (km)'); ylabel('m m^{-1}');    
    
% 3. Calculate resistive stress

    % calculate spatial grid spacing dx (m)
    dx = zeros(1,length(x));
    dx(1) = x(2)-x(1); % forward difference
    dx(2:end-1) = x(3:end)-x(1:end-2); % central difference
    dx(end) = x(end)-x(end-1); % backward difference   
    
    % estimate effective viscosity v (Pa s)
    v = ((E.*A).^(-1/n)).*(abs(dUdx_lin)).^((1-n)/n);
    
    % calculate the effective pressure N (ice overburden pressure minus water
    % pressure) assuming an easy & open connection between the ocean and
    % ice-bed interface
    sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
    N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
    N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
    N = [N_ground N_marine];
    N(N<0)=0; % cannot have negative values
    
    % calculate along-flow resistive stress as a sum of the local 
    % longitudinal stresses and basal drag
    % Td = Tlat + Tb + Tlon
    % Rxx = Tlon + Tb
    Td = rho_i*g.*H.*dhdx; % driving stress (Pa)
    Tlat = H./W.*(5.*U./(E.*A.*W)).^(1/n); % lateral resistance (Pa)
    Tlon = (2./dx.*(H.*v.*dUdx_lin)); % longitudinal resistance (Pa)
    Tb = -Td-Tlat+Tlon; % basal drag (Pa)
    Rxx = Tlon+Tb;  % (Pa) 
        Rxx(end-2:end)=Rxx(end-3);
    
    % fit a fourier transform to the data to estimate period
    fitmethod = 'poly7';
    f_Rxx = feval(fit(x',Rxx',fitmethod),x);
    f_H = feval(fit(x',H',fitmethod),x);
    [~,RxxTlocs]=findpeaks(f_Rxx);
    T_Rxx = mean(diff(x(RxxTlocs))*0.1); 
    if any(isnan(T_Rxx))
        T_Rxx = x(RxxTlocs);
    end
    [~,Hlocs]=findpeaks(f_H);
    T_H = mean(diff(x(Hlocs))*0.1);    
    if isnan(T_H)
        T_H = x(Hlocs);
    end
    % display
    disp(['T_Rxx = ',num2str(round(T_Rxx)),' m']);
    disp(['T_H = ',num2str(round(T_H)),'m']);   
    disp(['mean(SCL:H) = ',num2str(T_Rxx/nanmean(H))]);
    
    % plot results
    figure(2); clf
    set(gcf,'position',[750 80 615 715]);        
    subplot(2,1,1);
        set(gca,'fontsize',14,'fontname','arial','linewidth',2);
        hold on; grid on; 
        plot(x/10^3,Rxx./10^3,'linewidth',2);
        plot(x/10^3,f_Rxx/10^3,'linewidth',2);
        for i=1:length(RxxTlocs)
            plot([x(RxxTlocs(i)) x(RxxTlocs(i))]/10^3,[-2000 2000],'--','color',[136 86 167]/255,'linewidth',2);
        end
        ylim([min(Rxx)/10^3 max(Rxx)/10^3]);
        ylabel('R_{xx} (kPa)');
    subplot(2,1,2);
        set(gca,'fontsize',14,'fontname','arial','linewidth',2);
        hold on; grid on; 
        plot(x/10^3,H,'linewidth',2); 
        plot(x/10^3,f_H,'linewidth',2);
        for i=1:length(RxxTlocs)
            plot([x(RxxTlocs(i)) x(RxxTlocs(i))]/10^3,[-2000 2000],'--','color',[136 86 167]/255,'linewidth',2);
        end
        ylim([min(H)-100 max(H)+100]);
        xlabel('distance along centerline (km)'); ylabel('Thickness (m)'); 
        
% save resulting resistive stress, peaks, fit        
if save_results
    cd([homepath,'inputs-outputs/']);
    save('Crane_SCL_results.mat','x','Rxx','RxxTlocs','f_Rxx');
    disp('results saved.');
end




    