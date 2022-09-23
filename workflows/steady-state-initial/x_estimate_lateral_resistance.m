%% Script to estimate lateral resistance:
%   1. Pre-collapse: between the 2002 terminus position and the end of the fjord
%   2. Pre-collapse: grounding line --> 2002 terminus position
%   3. Model simulations: furthest retreated (grounding line -> calving front)
%   4. Model simulations: furthest advanced (grounding line -> calving front)
% Rainey Aberle
% Fall 2022

% lateral resistance = 2*H/W * (5*U/E*A*W)^(1/n)

clear all; close all;

% Define path in directory to CraneGlacier_flowlinemodeling
basepath = '/Users/raineyaberle/Research/MS/CraneGlacier_flowlinemodeling/';
cd([basepath,'inputs-outputs/']);

% Load model initialization variables
load('model_initialization_pre-collapse.mat');

% Load centerline
cl.x = load('Crane_centerline.mat').x; cl.y = load('Crane_centerline.mat').y; 
if size(cl.x)==[186 1]
    cl.x=cl.x';
end
if size(cl.y)==[186 1]
    cl.y=cl.y';
end

% Define necessary parameters 
rho_i = 917; % density of ice [kg/m^3]
rho_sw = 1000; % density of seawater [kg/m^3]
n = 3; % flow law exponent [unitless]

%% 1. Pre-collapse: between the 2002 terminus position and the end of the fjord

% end of fjord location index
f0 = 295; 

% dx
dx = 200;

% floating thickness from c0:f0
% rho_i*H_i = rho_sw*H_sw --> H = -rho_sw/(rho_i - rho_sw) * h
H = -rho_sw/(rho_i-rho_sw) * h0(c0:f0);

% lateral resistance
Rxy = 2*H./W0(c0:f0) .* nthroot(5.*U0(c0:f0) ./ (A0(c0:f0).*W0(c0:f0)), n);
Rxy_cf_fm_sum = sum(Rxy); % Pa
disp(['Sum Rxy from 2002 cf position to fjord mouth = ',num2str(Rxy_cf_fm_sum/10^3),' kPa']);

% plot
figure(1); clf;
hold on; set(gca,'fontsize',12,'linewidth',1);
plot(x0(c0:f0)/10^3, Rxy/10^3, '-b', 'linewidth',2, 'displayname', 'R_{xy}');
grid on; legend;
xlabel('distance along centerline [km]');
ylabel('lateral resistance [kPa]');

figure(2); clf
subplot(1,2,1); 
hold on; set(gca,'fontsize',12,'linewidth',1);
plot(x0/10^3, U0.*3.1536e7, '-k','linewidth',2);
xlabel('distance along centerline [km]');
ylabel('speed [m/y]');
grid on
subplot(1,2,2);
hold on; set(gca,'fontsize',12,'linewidth',1);
plot(x0/10^3, h0, '-b','linewidth',2);
plot(x0(c0:f0)/10^3, h0(c0:f0)-H,'-c','linewidth',2);
plot(x0/10^3, b0,'-k', 'linewidth',2);
xlabel('distance along centerline [km]');
ylabel('elevation [m]');
grid on

%% 2. Pre-collapse: between the most retreated and the 2002 calving front positions

% initial thickness (grounded everywhere)
H0 = h0-b0; % ice thickness [m]

% load calving front positions
term = load('observed_terminus_positions.mat').term;
c_min = dsearchn(x0',term.x(4));

% lateral resistance
Rxy = 2*H0(c_min:c0)./W0(c_min:c0) .* nthroot(5.*U0(c_min:c0) ./ (A0(c_min:c0).*W0(c_min:c0)), n);
Rxy_2007_sum = sum(Rxy); % Pa
disp(['Sum Rxy from most retreated to 2002 cf position = ',num2str(Rxy_2007_sum/10^3),' kPa']);

% plot
figure(1); clf;
hold on; set(gca,'fontsize',12,'linewidth',1);
plot(x0(c_min:c0)/10^3, Rxy/10^3, '-b', 'linewidth',2, 'displayname', 'R_{xy}');
grid on; legend;
xlabel('distance along centerline [km]');
ylabel('lateral resistance [kPa]');

%% 3. Model simulations: 2018 (grounding line -> calving front)

% load modeled 2018 conditions
mod_cond = load('modeled_conditions_2007-2018.mat').mod_cond;
H_2018 = mod_cond(end).H;
x_2018 = mod_cond(end).x;
U_2018 = mod_cond(end).U;
gl_2018 = mod_cond(end).gl;
c_2018 = mod_cond(end).c;
W_2018 = interp1(x0, W0, x_2018);
A_2018 = interp1(x0, A0, x_2018);

% lateral resistance
Rxy_2018 = (2*H_2018(gl_2018:c_2018)./W_2018(gl_2018:c_2018) ...
    .* nthroot(5.*U_2018(gl_2018:c_2018) ./ (A_2018(gl_2018:c_2018).*W_2018(gl_2018:c_2018)), n));
Rxy_2018_sum = sum(Rxy_2018); % Pa
disp(['Sum Rxy_2018 = ',num2str(Rxy_2018_sum/10^3),' kPa']);

%% 4. Model simulations: furthest advanced for unperturbed scenarios(grounding line -> calving front)

% load modeled 2100 conditions under the highest SMB_enh and F_T
% perturbation scenario
SMBenh_TF_2100 = load([basepath,'workflows/steady-state-initial/results/1_SMB_DFW_TF/SMB0_DFW5m_TF0_geom.mat']);
H_2100 = SMBenh_TF_2100.H2; 
x_2100 = SMBenh_TF_2100.x2; 
U_2100 = SMBenh_TF_2100.U2; 
gl_2100 = SMBenh_TF_2100.gl2; 
c_2100 = SMBenh_TF_2100.c2; 
A_2100 = interp1(x0, A0, x_2100);
W_2100 = interp1(x0, W0, x_2100);

% lateral resistance
Rxy_2100 = (2*H_2100(gl_2100:c_2100)./W_2100(gl_2100:c_2100) ...
    .* nthroot(5.*U_2100(gl_2100:c_2100) ./ (A_2100(gl_2100:c_2100).*W_2100(gl_2100:c_2100)), n));
Rxy_2100_sum = sum(Rxy_2100); % Pa
disp(['Sum Rxy_2100 = ',num2str(Rxy_2100_sum/10^3),' kPa']);

% print results
disp('  ');
disp('Lateral resistance lost after ice shelf collapse:')
disp([num2str((Rxy_cf_fm_sum+3000+Rxy_2007_sum)/10^3), ' kPa']);
disp('  ');
disp('Lateral resistance between grounding line and most advanced cf position:')
disp([num2str(Rxy_2100_sum/10^3),' kPa']);



