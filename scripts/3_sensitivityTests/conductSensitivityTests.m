%% Glacier Flowline Model: Sensitivity Tests for Crane Glacier, Antarctic Peninsula
% Rainey Aberle 
% Fall 2021
% 
%   0. Define time and space independent variables by loading the flowline 
%       initialization file and regridding to the desired grid spacing. 
%   1. Run sensitivity tests for surface mass balance (SMB), submarine
%       melting rate (SMR), and fresh water depth in crevasses (FWD) independently. 

%% 0. define necessary paths in directory

clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory and add necessary paths
homepath = '/Users/raineyaberle/Desktop/Research/CraneModeling/CraneGlacier_flowlinemodeling/';
addpath([homepath,'scripts/']); % add path to U_convergence
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean']);
addpath([homepath,'inputs-outputs/']);

%% 1. conduct sensitivity tests 
% (A) Run through the the designated number of model years with no change 
% (B) Run through simulated 2009-2019, then continue evolving the model  
%   under new scenarios: 
%   (1) SMB, DFW, and TF independent change
%   -- Switch between scenarios on lines 63-65 below.
%   (2) SMB_enh: increased air temperature --> increased SMR due to 
%   increased surface runoff
%   -- Run the SMB independent scenario with SMB_enhance = 1.
%   (3) SMB_enh + TF: increased air and ocean temperature --> increased 
%   surface runoff --> increased SMR AND increased TF
%   -- Run the SMB_enh scenario and TF scenario simulataneously

close all; 

saveFinal = 1;         % = 1 to save final conditions in homepath/3_sensitivityTests/results/
plotTimeSteps = 0;     % = 1 to plot geometry, speed, cf/gl positions every decade
plotMisfits = 0;       % = 1 to plot misfit with 2018 conditions
plotClimateParams = 1; % = 1 to plot climate parameters

% load no change conditions
load('2100_noChange.mat'); % load no change variables
    
% set up changes in SMB, DFW, & TF
% - decrease maximum SMB by increments of 0.5 m a-1 down to -10 m a-1
% - increase DFW by increments of 1 m up to 10 m
% - increase TF by increments of 0.1 ^oC up to 1 ^oC
delta_SMB0 = (0:-1:-10)./3.1536e7; % m/s change in SMB at the calving front (used to increase gradient)
delta_DFW0 = 0:10; % m change in DFW 
delta_TF0 = 0:0.1:1; % ^oC change in TF

% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 91*3.1536e7;

% loop through scenarios
for j=1:length(delta_SMB0)

    load('flowlineModelInitialization.mat');
    
    % Switch scenarios on and off
    delta_SMB = delta_SMB0(j);
    delta_DFW = 0;%delta_DFW0(j);
    delta_TF = delta_TF0(j);
    SMB_enhance = 1; % = 1 to increase SMR due to decreased SMB    

    %try
        % run flowline model
        [x,U,h,hb,H,gl,c,xcf,dUdx,Fgl,XCF,XGL] = flowlineModel(homepath,plotTimeSteps,plotMisfits,plotClimateParams,dt,t_start,t_end,beta0,DFW0,delta_SMB,delta_DFW,delta_TF,SMB_enhance);
        
        % save geometry & speed
        if saveFinal
            % store final geometry & speed
            h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; DFW2 = DFW0+delta_DFW;
            % no change
            if delta_DFW==0 && delta_TF==0 && delta_SMB==0 
                cd([homepath,'scripts/3_sensitivityTests/results/1_SMB_DFW_TF/']);
                save('SMB0_DFW0m_TF0_geom.mat','h2','H2','c2','U2','gl2','x2','DFW2','Fgl2','XGL2','XCF2');
                cd([homepath,'inputs-outputs/']);
                h1=h; H1=H; hb1=hb; c1=c; U1=U; x1=x; gl1=dsearchn(x1',XGL(end)); DFW1=DFW0; Fgl1=Fgl; XGL1=XGL; XCF1=XCF; 
                save('2100_noChange.mat','h1','H1','hb1','c1','U1','gl1','x1','DFW1','Fgl1','XGL1','XCF1');
                disp('2100_noChange saved');
            % 2) SMB_enhanced
            elseif SMB_enhance==1 && delta_TF==0
                cd([homepath,'scripts/3_sensitivityTests/results/2_SMB_enh/']);
                fileName = ['SMB',num2str(delta_SMB*3.1536e7),'_enh_geom.mat'];
                save(fileName,'h2','H2','c2','U2','gl2','x2','DFW2','Fgl2','XGL2','XCF2');
                disp('geometry saved (2)');
            % 3) SMB_enhanced + TF
            elseif SMB_enhance==1 && delta_TF~=0 
                cd([homepath,'scripts/3_sensitivityTests/results/3_SMB_enh+TF/']);
                fileName = ['SMB',num2str(delta_SMB*3.1536e7),'_enh_dTF',num2str(delta_TF),'geom.mat'];
                save(fileName,'h2','H2','c2','U2','gl2','x2','DFW2','Fgl2','XGL2','XCF2');
                disp('geometry saved (3)');           
            % 1) SMB, DFW, and TF independent
            else
                cd([homepath,'scripts/3_sensitivityTests/results/1_SMB_DFW_TF/']);
                fileName = ['SMB',num2str(delta_SMB*3.1536e7),'_DFW',num2str(DFW2),'m_TF',num2str(delta_TF),'_geom.mat'];
                save(fileName,'h2','H2','c2','U2','gl2','x2','DFW2','Fgl2','XGL2','XCF2');
                disp('geometry saved (1)');
            end
        else
            disp('geometry not saved.');
        end

        % plot results
        b=1;
        while b==1
            cd([homepath,'inputs-outputs/']);
            load('2100_noChange.mat');
            h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; DFW2=DFW0+delta_DFW; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; % store final geometry & speed
            gl1=dsearchn(x2',XGL(end)); 
            figure(10); clf % sensitivity test changes
            hold on; grid on;
            set(gcf,'Position',[491 80 886 686]);
            set(gca,'FontSize',14,'linewidth',2,'fontweight','bold');
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            legend('Location','east'); xlim([0 100]); ylim([min(hb0)-100 max(h)+100]);
            title(['TF = + ',num2str(round(delta_TF,1)),'^oC, SMB = + ',...
                num2str(round(delta_SMB*3.1536e7,1)),'m/a, DFW = ',num2str(DFW2),'m']);
            ax1=get(gca);
                % ice surface
                plot(x1(1:c1)/10^3,h1(1:c1),'-k','linewidth',2,'displayname','no change');
                plot(x2(1:c2)/10^3,h2(1:c2),'color',[0.8 0 0],'linewidth',2,'displayname','change');
                % calving front
                plot(x1(c1)*[1,1]/10^3,[h1(c1)-H1(c1),h1(c1)],'-k','linewidth',2,'HandleVisibility','off');
                plot(x2(c2)*[1,1]/10^3,[h2(c2)-H2(c2),h2(c2)],'color',[0.8 0 0],'linewidth',2,'HandleVisibility','off');
                % floating bed
                plot(x1(gl1:c1)/10^3,h1(gl1:c1)-H1(gl1:c1),'-k','linewidth',2,'HandleVisibility','off');
                plot(x2(gl2:c2)/10^3,h2(gl2:c2)-H2(gl2:c2),'color',[0.8 0 0],'linewidth',2,'HandleVisibility','off');
                % bed elevation
                plot(x1/10^3,hb1,'k','linewidth',2,'HandleVisibility','off');
            % inset plot of terminus
            ax2 = axes('Position',[0.62 0.62 0.28 0.28]);
                hold on; grid on; set(gca,'linewidth',2,'fontweight','bold','fontsize',9);
                delta_L = x2(c2)-x1(c1); % change in length (m)
                delta_H = mean(H2(1:c2))-mean(H1(1:c1)); % change in ice thickness (m)
                delta_U = mean(U2(1:c2))-mean(U1(1:c1)); % change in ice speed (m/s) 
                title(['\Delta L=',num2str(round(delta_L)),'m, ','\DeltaH_{\mu}=',...
                    num2str(round(delta_H)),' m, ','\DeltaU_{\mu}=',...
                    num2str(round(delta_U.*3.1536e7)),' m/a']);
                xlim([x(gl)/10^3-5 x(c)/10^3+5]); ylim([min(hb0)-100 300]);
                % ice surface
                plot(x1(1:c1)/10^3,h1(1:c1),'-k','linewidth',2,'displayname','no change');
                plot(x2(1:c2)/10^3,h2(1:c2),'color',[0.8 0 0],'linewidth',2,'displayname','no change');
                % calving front
                plot(x1(c1)*[1,1]/10^3,[h1(c1)-H1(c1),h1(c1)],'-k','linewidth',2,'HandleVisibility','off');
                plot(x2(c2)*[1,1]/10^3,[h2(c2)-H2(c2),h2(c2)],'color',[0.8 0 0],'linewidth',2,'HandleVisibility','off');
                % floating bed
                plot(x1(gl1:c1)/10^3,h1(gl1:c1)-H1(gl1:c1),'-k','linewidth',2,'HandleVisibility','off');
                plot(x2(gl2:c2)/10^3,h2(gl2:c2)-H2(gl2:c2),'color',[0.8 0 0],'linewidth',2,'HandleVisibility','off');
                % bed elevation
                plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off');
                % mean sea level
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off'); drawnow
            b=b+1;
        end

    %catch
    %    disp(['iteration ',num2str(j),' failed']);
    %end
    
end
