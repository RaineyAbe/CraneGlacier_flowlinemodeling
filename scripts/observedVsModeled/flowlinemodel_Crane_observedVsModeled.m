%% Glacier Flowline Model
% Script to tune the upper boundary and SMB to nearly match observations of
% ice surface elevation (and dHDt)
%
% Rainey Aberle
% Fall 2020
% Adapted from Ellyn Enderlin's flowline model demo code

clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)

%% define time and space independent variables
    
dx0 = 150; % desired grid spacing (m)
dx=dx0;
            
use_binavg = 1;     % = 1 to use average within bins proportional to dx0
figure_save = 1;    % = 1 to save figure

% define home path in directory
    homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/';
    cd([homepath,'inputs-outputs']);
    addpath([homepath,'scripts/tuningCalving']); % add path to U_convergence
    
% Load Crane Glacier initialization variables
    load('Crane_flowline_initialization.mat');
    n = find(isnan(W0),1,'first');
    W0(n:end) = W0(n-1).*ones(1,length(W0(n:end))); clear n
    A0(end+1:length(x0)) = A0(end).*ones(1,length(x0)-length(A0)); 
    
% Load observations
    % Surface
    h_obs0 = load('Crane_SurfaceObservations_2009-2018.mat').h;
    k = [1 2 3 5 10 13 15 18 29 36]; % indices of surfaces to use for annual tracking
    years=2009:2018; 
    h_obs = zeros(length(years),length(x0));
    for i=1:length(2009:2018)
        surfs = [];
        for j=1:length(h_obs0)
            if contains(h_obs0(j).date,num2str(years(i)))
                if j==16
                    continue;
                else
                    surfs = [surfs; h_obs0(j).surface];
                end
            end
        end
        h_obs(i,:) = mean(surfs,1,'omitnan');
    end
    h0 = h_obs(1,:); 
    clear k
    
    % Speed
    u_obs0 = load('Crane_CenterlineSpeeds_2007-2017.mat').U;
    k = [5 9 12 16 19 22 24 26 28];
    u_obs = zeros(length(years)-1,length(x0));
    for i=1:length(2009:2018)-1
        speeds = [];
        for j=1:length(u_obs0)
            if contains(num2str(fix(u_obs0(j).year)),num2str(years(i)))
                speeds = [speeds; u_obs0(j).speed'];
            end
        end
        if ~isempty(speeds)
            u_obs(i,:) = mean(speeds,1,'omitnan');
        end
    end
    U0 = feval(fit(x0(~isnan(u_obs(1,:)))',u_obs(1,~isnan(u_obs(1,:)))','linearinterp'),x0'); 

    % Terminus
    term = load('Crane_TerminusPosition_2002-2019.mat').term;
    for i=1:length(term)
        termx_obs(i) = term(i).x;
        termDate_obs(i) = term(i).decidate;
    end
    % start in 2009
    termx_obs(1:4)=[]; termDate_obs(1:4)=[];
    years = 2009:2018; 
    % fit a quadratic to the terminus positions
    termx_obs = feval(fit(termDate_obs',termx_obs','poly2'),years');
    clear term    
    
% densities and g
    rho_i = 917; % ice density (kg m^-3)
    rho_sw = 1028; % ocean water density (kg m^-3)
    rho_fw = 1000; % fresh water density (kg m^-3)
    g = 9.81; % acceleration (m s^-2)

% time stepping (s)
    dt = 0.1*3.1536e7; 
    t_start = 0*3.1536e7; 
    t_end = 10*3.1536e7;    
    t = (t_start:dt:t_end);

% stress parameters (unitless)
    m = 1; % basal sliding exponent
    n = 3; % flow law exponent
    E = 1; % enhancement factor

% calving parameters
    Hc = 400; % m -> set the calving front to a default minimum ice thickness value
    fwd = 30; % fresh water depth in crevasses (m)
    
% maximum thickness cut-off to check for instability
    H_max = 2000; %maximum thickness (m)

% regrid the initialization data to match the desired grid spacing, 
    xi = x0; % rename initialization distance vector
    L = 70e3; % length of glacier (m)
    xi = 0:dx0:L; % desired distance vector (m from ice divide)  
    
    % calving front location
    c = dsearchn(transpose(xi),x0(term_2009.x)); % 2009 terminus location (index)
         
    % If the desired grid spacing is smaller than the original, use the
    % interp1 function to determine each spatial vector.
    % Otherwise, take the average within each bin at every point for each
    % spatial variable. 
    if length(xi)>length(x0)
        hi = interp1(x0,h0,xi);
        hbi = interp1(x0,hb0,xi);
        Wi = interp1(x0,W0,xi);
        Ui = interp1(x0,U0,xi);
        Ai = interp1(x0,A0,xi);
        betai = interp1(x0,beta0,xi);            
    else
        % Use a staggered grid for bin averages
        xm = (xi(1:end-1)+xi(2:end))./2;
        for k=1:length(xi)
            if k==1
              hi(k) = mean(h0(1:dsearchn(x0',xm(k))));
              hbi(k) = mean(hb0(1:dsearchn(x0',xm(k))));
              Wi(k) = mean(W0(1:dsearchn(x0',xm(k))));
              Ui(k) = mean(U0(1:dsearchn(x0',xm(k))));
              Ai(k) = mean(A0(1:dsearchn(x0',xm(k))));
              betai(k) = mean(beta0(1:dsearchn(x0',xm(k))));               
            elseif k==length(xi)
              hi(k) = mean(h0(dsearchn(x0',xm(k-1)):c));
              hbi(k) = mean(hb0(dsearchn(x0',xm(k-1)):c));
              Wi(k) = mean(W0(dsearchn(x0',xm(k-1)):c));
              Ui(k) = mean(U0(dsearchn(x0',xm(k-1)):c));
              Ai(k) = mean(A0(dsearchn(x0',xm(k-1)):c));
              betai(k) = mean(beta0(dsearchn(x0',xm(k-1)):c));              
            else
              hi(k) = mean(h0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
              hbi(k) = mean(hb0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
              Wi(k) = mean(W0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
              Ui(k) = mean(U0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
              Ai(k) = mean(A0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
              betai(k) = mean(beta0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));                            
            end
        end
    end 
    
    dUdxi = [(Ui(2:end)-Ui(1:end-1))./(xi(2:end)-xi(1:end-1)) 0]; % strain rate
    Hi = hi-hbi; % thickness (m)    

    % add the calving front conditions for each spatial variable
    hi(c+1:length(xi)) = zeros(1,length(hi(c+1:length(xi))));
    Ui(c+1:length(xi)) = zeros(1,length(Ui(c+1:length(xi))));
    Hi(c+1:length(xi)) = zeros(1,length(Hi(c+1:length(xi))));
    dUdxi(c+1:length(xi)) = zeros(1,length(dUdxi(c+1:length(xi))));
    betai(c+1:length(xi)) = zeros(1,length(betai(c+1:length(xi))));
    
    % find the location of the grounding line and the end of the ice-covered domain
    Hf = -(rho_sw./rho_i).*hbi; % flotation thickness (m)
    gl = find(Hf-Hi<0,1,'last')-1; %grounding line location 
    Hi(gl:end)=hi(gl:end)*rho_sw/(rho_sw-rho_i); % buoyant thickness using surface
    for i=1:length(xi) % thickness can't go beneath bed elevation
        if Hi(i) >= hi(i)-hbi(i)
            Hi(i) = hi(i)-hbi(i);
        end
    end    
    
    % add a dummy ice end (hi & Hi)
    for i=c+1:length(xi)
        %hi(i) = hi(i-1)-5; % decrease by 5m until at 0m  
        Hi(i) = Hi(i-1)-30; % decrease by 20m until 0 
        if Hi(i)>=hi(i)-hbi(i)
            Hi(i)=hi(i)-hbi(i); % can't go beneath bed elevation
        end
    end  
    hi(hi<0)=0; % surface can't go below sea level
    Hi(Hi<0)=0; % no negative thicknesses
    
    % find the end of the ice-covered domain (m along centerline)
    ice_end = (find(Hi<=0,1,'first')); 
    
    % extend other variables from c+1:ice_end (Ui,Ai)
    Ui(c+1:ice_end) = Ui(c).*ones(1,length(xi(c+1:ice_end)));
    Ai(c+1:ice_end) = Ai(c).*ones(1,length(xi(c+1:ice_end)));
    betai(gl:end) = zeros(1,length(betai(gl:end)));
    
    % use 90% the observed velocity at the upper bounds
    Ui(1:round(0.2*length(Ui))) = 0.9.*Ui(1:round(0.2*length(Ui)));
    
    % rename initial vectors
    x=xi; h=hi; hb=hbi; W=Wi; H=Hi; A=Ai; beta=betai; U=Ui; dUdx=dUdxi;
    
%% Run the flowline model   

cd /Users/raineyaberle/Desktop/Research/write-ups/AGU2020/ObsVid/

for i=1:length(t)

    % set up figures, plot geometries at t==0, then every t/10 iterations
    if t(i)==t_start
        ii=1; % index for full years
        years=2009:2019;
        col = parula(length(t)+15); %Color scheme for plots
        figure; set(gcf,'Position',[0 0 1100 1000]);
        % Modeled
        subplot(3,3,1); % modeled glacier geometry
            hold on; grid on;
            set(gca,'FontSize',11,'linewidth',1); 
            xlim([0 60]); ylim([min(hb)-100 max(h)+100]); title('Glacier Geometry'); 
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)'); 
            % Add text label
            ta = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
                ' a ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);
            % ice surface
            plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','2009');
            % calving front
            plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');          
            % floating bed
            plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
            if i==1
                % bed elevation
                plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off'); 
                % mean sea level
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
            end
        subplot(3,3,2); % modeled ice surface speed
            hold on; grid on; 
            set(gca,'FontSize',11,'linewidth',1); 
            title('Ice Surface Speed');  
            xlim([0 60]); ylim([0 2200]);
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})'); 
            plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','2009');
            % Add text label            
            tb = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
                ' b ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);           
        subplot(3,3,3); % modeled terminus position
            hold on; grid on; 
            set(gca,'FontSize',11,'linewidth',1);
            title('Terminus Position'); ylabel('Year');
            xlabel('Distance Along Centerline (km)'); 
            xlim([42 52]); ylim([2008 2020]); legend('Location','eastoutside');
            % 2009 terminus position
            plot(x(c)./10^3,2009,'.','markersize',20,'color',col(i,:),...
                'linewidth',1.5,'displayname','2009'); 
            % Add text label            
            tc = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.86+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
                ' c ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);   
       % Observed
       term = dsearchn(x0',termx_obs(i)); % terminus index on original grid
       subplot(3,3,4); % observed glacier geometry
            hold on; grid on; set(gca,'FontSize',11,'linewidth',1); 
            xlim([0 60]); ylim([min(hb)-100 max(h)+100]);
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)'); 
            % ice surface
            plot(x0(1:term)./10^3,h_obs(i,1:term),'color',col(i,:),'linewidth',2);
            % calving front
            plot(x0(term)*[1,1]/10^3,[hb0(term),h0(term)],'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');          
            if i==1
                % bed elevation
                plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off');
                % mean sea level
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off'); 
            end
            % Add text label            
            td = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
                ' d ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);            
       subplot(3,3,5); % observed ice surface speed 
            hold on; grid on; set(gca,'FontSize',11,'linewidth',1); 
            xlim([0 60]); ylim([0 2200]);
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})'); 
            plot(x0(1:term)./10^3,u_obs(i,1:term).*3.1536e7,'color',col(i,:),'linewidth',2);
            % Add text label            
            te = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
                ' e ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);            
       subplot(3,3,6); % observed terminus position
            hold on; grid on; set(gca,'FontSize',11,'linewidth',1);
            xlabel('Distance Along Centerline (km)'); ylabel('Year');
            xlim([42 52]); ylim([2008 2020]); legend('Location','eastoutside');
            plot(x0(term)./10^3,years(i),'.','markersize',20,'color',col(i,:),...
                'linewidth',1.5,'displayname',num2str(years(i)));      
            % Add text label            
            tf = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.88+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
                ' f ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);           
        % Misfit
        h_misfit_mean = zeros(1,length(years)-1);
        u_misfit_mean = zeros(1,length(years)-2);
        term_misfit = zeros(1,length(years));
        subplot(3,3,7); % surface misfit
            hold on; grid on; set(gca,'fontsize',11,'linewidth',1); 
            xlim([-100 100]); ylim([2008 2020]);
            ylabel('Year'); xlabel('Mean Misfit (m)'); 
            hobs = ~isnan(h_obs(i,:));
            h_misfit_mean(i) = mean(h(hobs)-interp1(x0(hobs)',h_obs(i,hobs)',x(hobs)),'omitnan');
            plot(h_misfit_mean(i),years(i),'x','markersize',15,'color',...
                col(i,:),'linewidth',2);  
            % Add text label            
            tg = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
                ' g ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);           
        subplot(3,3,8); %speed misfit
            hold on; grid on; set(gca,'fontsize',11,'linewidth',1); 
            ylim([2008 2020]); xlim([-600 100]);
            ylabel('Year'); xlabel('Mean Misfit (m yr^{-1})');
            uobs = ~isnan(u_obs(ii,:));
            u_misfit_mean(ii) = mean(U'-interp1(x0(uobs)',u_obs(ii,uobs)',x','linear'),'omitnan');
            plot(u_misfit_mean(i).*3.1536e7,years(i),'x','markersize',15,'color',...
                col(i,:),'linewidth',2);  
            % Add text label            
            th = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
                ' h ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);
        subplot(3,3,9); % terminus misfit
            hold on; grid on; set(gca,'fontsize',11,'linewidth',1); 
            ylabel('Year'); xlabel('Mean Misfit (km)'); legend('Location','eastoutside');
            ylim([2008 2020]); xlim([-4 3]);
            term_misfit(i) = x(c)-x0(term); 
            plot(term_misfit(i)./10^3,years(i),'x','markersize',15,'color',...
                col(i,:),'linewidth',2,'displayname',num2str(years(i))); 
            % Add text label            
            ti = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.88+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
                ' i ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);           
       % Insert text annotations
       tmod = annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Modeled', ...
            'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.05 .87 0 0],...
            'FontSize',18,'FontName','Arial','fontweight','bold','TextColor',[0.75 0 0]);
       tobs = annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Observed', ...
            'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.05 .58 0 0],...
            'FontSize',18,'FontName','Arial','fontweight','bold','TextColor',[0.75 0 0]);
       tmis = annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Misfit', ...
            'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.05 .26 0 0],...
            'FontSize',18,'FontName','Arial','fontweight','bold','TextColor',[0.75 0 0]);
        if figure_save==1
             % save image on each iteration
                saveas(gcf,[num2str(t(i)./3.1536e7+2009),'.png'],'png');
                disp([num2str(t(i)./3.1536e7+2009),'.png saved']);
        end
        
    elseif mod(i-1,round(length(t)/10))==0
        ii = ii+1;
        figure(1); 
        % Modeled
        subplot(3,3,1); % modeled glacier geometry
            hold on; grid on;
            % ice surface
            plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','2009');
            % calving front
            plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');          
            % floating bed
            plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
        subplot(3,3,2); % modeled ice surface speed
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})'); 
            plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','2009');     
        subplot(3,3,3); % modeled terminus position
            plot(x(c)/10^3,years(ii),'.','markersize',20,'color',col(i,:),'linewidth',1.5,'displayname',num2str(t(i)./3.1536e7+2009)); 
        % Observed
            if ii<11
                term = dsearchn(x0',termx_obs(ii)); % terminus index on original grid
                subplot(3,3,4); % observed glacier geometry
                % ice surface
                plot(x0(1:term)./10^3,h_obs(ii,1:term),'color',col(i,:),'linewidth',2);
                % calving front
                plot(x0(term)*[1,1]/10^3,[hb0(term),h0(term)],'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');          
            end
            if ii<=9
                subplot(3,3,5); % observed ice surface speed 
                xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})'); 
                plot(x0(1:term)./10^3,u_obs(ii,1:term).*3.1536e7,'color',col(i,:),'linewidth',2);
            end
            subplot(3,3,6); % observed terminus position
                plot(x0(term)./10^3,years(ii),'.','markersize',20,'color',col(i,:),'linewidth',1.5,'displayname',num2str(years(ii)));  
        % Misfit
            if ii<11
                subplot(3,3,7); % surface misfit
                hobs = ~isnan(h_obs(ii,:));
                h_misfit_mean(ii) = mean(h(hobs)-interp1(x0(hobs)',h_obs(ii,hobs)',x(hobs),'linear'),'omitnan');
                plot(h_misfit_mean(ii),years(ii),'x','markersize',15,'color',...
                col(i,:),'linewidth',2);  
            end
            if ii<=9
                subplot(3,3,8); % speed misfit
                uobs = ~isnan(u_obs(ii,:));
                u_misfit_mean(ii) = mean(U'-interp1(x0(uobs)',u_obs(ii,uobs)',x','linear'),'omitnan');
                plot(u_misfit_mean(ii).*3.1536e7,years(ii),'x','markersize',15,'color',...
                col(i,:),'linewidth',2);         
            end
            subplot(3,3,9); % terminus misfit
            term_misfit(ii) = x(c)-x0(term);
            plot(term_misfit(ii)/10^3,years(ii),'x','markersize',15,'color',...
                col(i,:),'linewidth',2,'displayname',num2str(years(ii)));  
       if figure_save==1
           % save image on each iteration
             if years(ii)==2019
                saveas(gcf,[num2str(t(i)./3.1536e7+2009),'.png'],'png');
                saveas(gcf,[num2str(t(i)./3.1536e7+2009),'_2.png'],'png');
                saveas(gcf,[num2str(t(i)./3.1536e7+2009),'_3.png'],'png');
                disp([num2str(t(i)./3.1536e7+2009),'.png saved']);
             else
                saveas(gcf,[num2str(t(i)./3.1536e7+2009),'.png'],'png');
                disp([num2str(t(i)./3.1536e7+2009),'.png saved']);
             end
       end
    end 

    % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; %flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding line to the calving front
        gl = find(Hf-H<0,1,'last')-1; %grounding line location 
        ice_end = (find(H<=0,1,'first')); %end of ice-covered domain
        if isempty(ice_end)
            ice_end = length(x);
            disp('ice end criteria not met.');
        end

    % calculate the glacier's surface elevation and slope
        h = hb+H; %h = surface elevation (m a.s.l.)
        h(gl:length(x)) = (1-rho_i/rho_sw).*H(gl:length(x)); %adjust the surface elevation of ungrounded ice to account for buoyancy
        dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; %surface slope (unitless)

    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E*A(1:length(dUdx)))),n); %resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); %crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); %calving front located where the inland-most crevasse intersects sea level
        %if the crevasses never intersect sea level
        if isempty(c) == 1
            c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
        end
        if isempty(c)==1
            c=length(x);
            disp('calving criteria not met');
        end
        % if the crevasses first intersect sea level inland of the grounding line 
        if c <= gl
            c = gl; %set the grounding line as the calving front
        end
        % use observed terminus position for first time increment
        if i==1
            c=dsearchn(transpose(x),x0(term_2009.x));            
        end 

    %calculate the effective pressure (ice overburden pressure minus water
    %pressure) assuming an easy & open connection between the ocean and
    %ice-bed interface
    sl = find(hb<=0,1,'first'); %find where the glacier base first drops below sea level
    N_ground = rho_i*g*H(1:sl); %effective pressure where the bed is above sea level (Pa)
    N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); %effective pressure where the bed is below sea level (Pa)
    N = [N_ground N_marine];
    N(N<0)=1; %cannot have negative values

    % Solve for new velocity
    [U,dUdx,T] = U_convergence(x,U,u_obs,dUdx,dhdx,H,A,E,N,W,dx,c,ice_end,n,m,beta,rho_i,rho_sw,g); 

    % calculate ice flux
    F = U.*H.*W; % ice flux (m^3 s^-1)
    F(isnan(F))=0;

    % calculate the  change in ice thickness from continuity
    dHdt = -(1./W).*gradient(F,x);
    dH = dHdt.*dt;

    % surface mass balance
    yr = round(t(i)/3.1536e7)+1;
    clear smb sigma_smb % clear to avoid changing size with changing x

    % interpolate smb0 to centerline, add tributary flux Q0 to smb
    smb = [interp1(x0,smb0(yr).smb_linear'+Q0,x)]./3.1536e7; % m/s
        smb(c+1:length(x)) = zeros(1,length(smb(c+1:length(x)))); % zeros past c
    sigma_smb = [interp1(x0,smb0(yr).sigma_smb+Q0_err,x)]./3.1536e7; % m/s
        sigma_smb(c+1:length(x)) = zeros(1,length(smb(c+1:length(x)))); % zeros past c

    % adjust smb to minimize misfit of surface observations 
    smb = smb-0.1e-5; % m/s
    smb(1:30) = smb(1:30)-0.12e-5; 
    smb(50:70) = smb(50:70)+0.05e-5;
    smb(50:100) = smb(50:100)+0.08e-5; 
    smb(115:290) = smb(115:290)-0.06e-5; 

    % new thickness (change from dynamics, SMB, & submarine melting)
    Hn = H+dH+(smb.*dt); 
    Hn(Hn < 0) = 0; % remove negative values 
    %Hn(c+1:end) = 0; % remove thickness values past calving front
    H = Hn; %set as the new thickness value

    % stop the model if it behaves unstably (monitored by ice thickness)
    if max(H) > H_max
        disp(['Adjust dt']);
        break
    end

    % find the precise location of the grounding line (where H=Hf)
    %xf = interp1((H(1:ice_end)-Hf(1:ice_end)),x(1:ice_end),0,'spline','extrap'); 
    xf = find(Hf-H>0,1,'first')-1;

    %adjust the grid spacing so the grounding line is continuously tracked
    xl = round(xf/dx0); %number of ideal grid spaces needed to reach the grounding line
    dx = xf/xl; %new grid spacing (should be ~dx0)
    xn = 0:dx:L; %new distance vector    

    %adjust the space-dependent variables to the new distance vector
    hb = interp1(x0,hb0,xn);
    W = interp1(x0,W0,xn);
    H = interp1(x,H,xn,'linear','extrap'); %ice thickness (m)
    Hf = interp1(x,Hf,xn,'linear','extrap');
    U = interp1(x(1:c),U(1:c),xn,'linear','extrap'); %speed (m s^-1)
    A = interp1(x,A,xn,'linear','extrap'); %rate factor (Pa^-n s^-1)
    beta = interp1(x,beta,xn,'linear','extrap'); % basal roughness factor

    %find the location of the grounding line and end of the ice-covered domain for the adjusted data
    gl = find(Hf-H<0,1,'last');
    ice_end = (find(H<=0,1,'first'));
    if isempty(ice_end)
        ice_end = length(x);
        disp('ice end criteria not met.')
    end

    %rename the distance vector
    x = xn; %distance from the divide (m)

    %calculate the new surface elevation and slope
    h = hb+H; %grounded ice surface elevation (m a.s.l.)
    h(gl:length(x)) = (1-rho_i/rho_sw).*H(gl:length(x)); %floating ice surface elevation (m a.s.l.)

    % calculate new strain rate
    dUdx = [(U(2:end)-U(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % strain rate
    %dUdx(c+1:end) = zeros(1,length(dUdx(ice_end+1:end))); 

    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E*A(1:length(dUdx)))),n); %resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); %crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); %calving front located where the inland-most crevasse intersects sea level
        %if the crevasses never intersect sea level
        if isempty(c) == 1
            c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
        end
        if isempty(c)==1
            c=length(x);
            disp('calving criteria not met');
        end
        % if the crevasses first intersect sea level inland of the grounding line 
        if c <= gl
            c = gl; %set the grounding line as the calving front
        end
        % use observed terminus position for first time increment
        if i==1
            c=dsearchn(transpose(x),x0(term_2009.x));            
        end 


end 
