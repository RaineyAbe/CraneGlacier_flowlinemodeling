%% Script to save files associated with data publishing
% Rainey Aberle
% 2022

clear all; close all;
base_path = '/Users/raineyaberle/Research/MS/CraneGlacier_flowlinemodeling/';
out_path = [base_path, 'data_package/'];

%% 1. Observations along the centerline

% -----Centerline
X0 = load([base_path,'inputs-outputs/Crane_centerline.mat']).x; 
Y0 = load([base_path,'inputs-outputs/Crane_centerline.mat']).y;
% define x as distance along centerline
x0 = zeros(1,length(X0));
for i=2:(length(X0))
    x0(i) = sqrt((X0(i)-X0(i-1))^2 + (Y0(i)-Y0(i-1))^2) + x0(i-1);
end

% -----Bed elevations
b0 = load([base_path,'inputs-outputs/observed_bed_elevations.mat']).HB.hb0; % OIB picks
b0_width_ave = load([base_path,'inputs-outputs/observed_bed_elevations_width_averaged.mat']).b_adj; % width-averaged
b0_width_ave_smooth = movmean(b0_width_ave, 20); % width-averaged, smoothed
b0_bathym = load([base_path,'inputs-outputs/observed_bed_elevations.mat']).bathym;

% -----Surface elevations
h0 = load([base_path, 'inputs-outputs/observed_surface_elevations.mat']).h;
h0(1) = [];

% -----Width
W0 = load([base_path, 'inputs-outputs/observed_glacier_width.mat']).width.W; 

% -----Surface speeds
U0 = load([base_path, 'inputs-outputs/observed_surface_speeds.mat']).U; 

% -----Calving front positions
xcf0 = load([base_path, 'inputs-outputs/observed_terminus_positions.mat']).term;

% -----Interpolate to 100 m grid spacing
L = x0(end); % length
dx = 100; % grid spacing
x = 0:dx:L; % new x grid
cl.x = x; 
cl.Easting = interp1(x0, X0, x);
cl.Northing = interp1(x0, Y0, x);
% bed elevation
b(1).b_centerline = interp1(x0, b0, x);
b(1).date = '20181016';
b(1).source = 'NASA OIB manual picks';
b(1).units = 'm';
b(2).b_centerline = interp1(x0, b0_width_ave, x);
b(2).date = '20181016';
b(2).source = 'NASA OIB manual picks, width-averaged';
b(2).units = 'm';
b(3).b_centerline = interp1(x0, b0_width_ave_smooth, x);
b(3).date = '20181016';
b(3).source = 'NASA OIB manual picks, width-averaged, smoothed';
b(3).units = 'm';
b(4).b_centerline = interp1(x0, b0_bathym, x);
b(4).date = '2005-2006';
b(4).source = 'Fjord bathymetry (Rebesco et al., 2014)';
b(4).units = 'm';
% surface elevations
for i=1:length(h0)
    h(i).date = h0(i).date;
    h(i).units = h0(i).units;
    h(i).source = h0(i).source;
    if length(h0(i).h_centerline) > 1
        I = find(~isnan(h0(i).h_centerline));
        h(i).h_centerline = interp1(x0(I), h0(i).h_centerline(I), x);
    else
        h(i).h_centerline = NaN;
    end
    if length(h0(i).h_width_ave) > 1
        I = find(~isnan(h0(i).h_width_ave));
        h(i).h_width_ave = interp1(x0(I), h0(i).h_width_ave(I), x);
    else
        h(i).h_width_ave = NaN;
    end

end
% width
W.width = interp1(x0, W0, x);
W.date = '20191013';
W.source = 'Landsat 8 manual picks';
W.units = 'm';
% surface speeds
for i=1:length(U0)
    U(i).date = U0(i).date;
    U(i).units = "m/yr";
    U(i).source = U0(i).source;
    if length(U0(i).U_centerline) > 1
        I = find(~isnan(U0(i).U_centerline));
        if length(I) > 1
            U(i).U_centerline = interp1(x0(I), U0(i).U_centerline(I), x);
        else
            U(i).U_centerline = NaN;
        end
    else
        U(i).U_centerline = NaN;
    end
    if length(U0(i).U_width_ave) > 1
        I = find(~isnan(U0(i).U_width_ave));
        if length(I) > 1
            U(i).U_width_ave = interp1(x0(I), U0(i).U_width_ave(I), x);
        else 
            U(i).U_width_ave = NaN;
        end
    else
        U(i).U_width_ave = NaN;
    end
    if length(U0(i).U_err) > 1
        I = find(~isnan(U0(i).U_err));
        U(i).U_err = double(interp1(x0(I), U0(i).U_err(I), x));
    else
        U(i).U_err = NaN;
    end
    % convert to m/yr
    U(i).U_centerline = double(U(i).U_centerline .* 3.1536e7);
    U(i).U_width_ave = U(i).U_width_ave .* 3.1536e7; 
    U(i).U_err = U(i).U_err .* 3.1536e7;
end
% calving front position
for i=1:length(xcf0.X)
    xcf(i).date = xcf0.date(i);
    xcf(i).Easting = xcf0.X(i);
    xcf(i).Northing = xcf0.Y(i);
    xcf(i).x = interp1(Y0, x0', xcf(i).Northing);
    xcf(i).source = 'Landsat 7-8 picks (Dryak and Enderlin, 2020)';
end

% ----Plot results
figure(1); clf;
subplot(2,1,1);
hold on; grid on;
plot(x/10^3, b(1).b_centerline, '-k', 'linewidth', 2);
col = parula(length(h));
for i=1:length(h)
    plot(x/10^3, h(i).h_centerline, 'color', col(i,:), 'linewidth',1);
end
col = winter(length(xcf));
for i=1:length(xcf)
    plot(xcf(i).x/10^3, 0, 'x', 'color', col(i,:), 'linewidth', 1, 'markersize', 10);
end
subplot(2,1,2);
hold on; grid on;
col = parula(length(U));
for i=1:length(U)
    plot(x/10^3, U(i).U_centerline, 'color', col(i,:), 'linewidth', 1);
end

% -----Save to file
save([out_path, 'observed_centerline.mat'], 'cl');
save([out_path, 'observed_bed_elevations.mat'], 'b');
save([out_path, 'observed_width.mat'], 'W');
save([out_path, 'observed_surface_elevations.mat'], 'h');
save([out_path, 'observed_calving_front_positions.mat'], 'xcf');
save([out_path, 'observed_surface_speeds.mat'], 'U');
cd(out_path);
disp('observed conditions saved to file.');

%% 2. Modeled conditions and parameters along the centerline

% -----Downscaled RACMO variables, interpolated to regridded centerline coordinates
load([base_path, 'inputs-outputs/downscaled_RACMO_variables_1995-2019.mat']);
for i=1:length(years)
    RACMO_ds(i).date = string(years(i));
    RACMO_ds(i).units = "m/yr";
    RACMO_ds(i).runoff = interp1(x0, ro.linear(i,:), x);
    RACMO_ds(i).snowfall = interp1(x0, sf.linear(i,:), x);
    RACMO_ds(i).snowmelt = interp1(x0, sm.linear(i,:), x);
    RACMO_ds(i).SMB = interp1(x0, smb.linear(i,:), x);
end

% -----Centerline averaged, interpolated to regridded centerline cooordinates
RACMO_ds(i+1).date = "time-averaged";
RACMO_ds(i+1).units = "m/yr";
RACMO_ds(i+1).runoff = interp1(x0, ro.downscaled_average_linear, x);
RACMO_ds(i+1).snowfall = interp1(x0, sf.downscaled_average_linear, x);
RACMO_ds(i+1).snowmelt = interp1(x0, sm.downscaled_average_linear, x);
RACMO_ds(i+1).SMB = interp1(x0, smb.downscaled_average_linear, x);

% -----Plot
figure(2); clf;
set(gcf, 'position', [-1000, 0, 960, 725]);
subplot(2, 2, 1); hold on; grid on; title('Runoff'); ylabel('m/yr'); legend('position', [0.04 0.25 0.02 0.6]);
subplot(2, 2, 2); hold on; grid on; title('Snowfall'); 
subplot(2, 2, 3); hold on; grid on; title('Snowmelt'); ylabel('m/yr'); xlabel('Distance along centerline [km]');
subplot(2, 2, 4); hold on; grid on; title('SMB'); xlabel('Distance along centerline [km]');
col = parula(length(RACMO_ds));
for i=1:length(RACMO_ds)
    if i<length(RACMO_ds)
        % runoff
        subplot(2, 2, 1);
        plot(x/1e3, RACMO_ds(i).runoff, 'color', col(i,:), 'linewidth',1, 'displayname', num2str(years(i)));
        % snowfall
        subplot(2, 2, 2);
        plot(x/1e3, RACMO_ds(i).snowfall, 'color', col(i,:), 'linewidth',1, 'displayname', num2str(years(i)));
        % snowmelt
        subplot(2, 2, 3);
        plot(x/1e3, RACMO_ds(i).snowmelt, 'color', col(i,:), 'linewidth',1, 'displayname', num2str(years(i)));
        % SMB
        subplot(2, 2, 4);
        plot(x/1e3, RACMO_ds(i).SMB, 'color', col(i,:), 'linewidth',1, 'displayname', num2str(years(i)));
    else
        % runoff
        subplot(2, 2, 1);
        plot(x/1e3, RACMO_ds(i).runoff, '-k', 'linewidth', 2, 'displayname', 'time-ave');
        % snowfall
        subplot(2, 2, 2);
        plot(x/1e3, RACMO_ds(i).snowfall, '-k', 'linewidth', 2, 'displayname', 'time-ave');
        % snowmelt
        subplot(2, 2, 3);
        plot(x/1e3, RACMO_ds(i).snowmelt, '-k', 'linewidth', 2, 'displayname', 'time-ave');
        % SMB
        subplot(2, 2, 4);
        plot(x/1e3, RACMO_ds(i).SMB, '-k', 'linewidth', 2, 'displayname', 'time-ave');
    end 
end

% -----Save to file
save([out_path, 'downscaled_RACMO_variables.mat'], 'RACMO_ds');
disp('downscaled_RACMO_variables.mat saved to file');


%% 3. Modeling results

results_path = [base_path, 'workflows/steady-state-initial/results/'];

% -----Define time variable
t_start = 0*3.1536e7; % 2002
t_mid = 18*3.1536e7; % 2020
t_end = 98*3.1536e7; % 2100
dt1 = 0.0005*3.1536e7; 
dt2 = 0.001*3.1536e7; 
t = [t_start:dt1:t_mid t_mid+dt2:dt2:t_end];
% convert to deciyears
t = t./3.1536e7 + 2002;

% -----Unperturbed
% loop through results, split into each scenario
cd([results_path, '1_SMB_DFW_TF/']);
files = dir('*geom.mat');
for i=1:length(files)
    if contains(files(i).name,'SMB0') && contains(files(i).name,'DFW0') && contains(files(i).name,'TF0_')
        files(i).change = 0;
        files(i).changeIn = 'SMB';
        files(i).name = files(i).name; 
        files(length(files)+1).change = 0;
        files(length(files)).changeIn = 'DFW';
        files(length(files)).name = files(i).name;
        files(length(files)+1).change = 0; 
        files(length(files)).changeIn = 'TF';
        files(length(files)).name = files(i).name;        
    % (a) SMB
    elseif contains(files(i).name,'DFW0m') && contains(files(i).name,'TF0_')
        files(i).change = str2double(files(i).name(regexp(files(i).name,'B')+1:...
            regexp(files(i).name,'_D')-1)); 
        files(i).changeIn = 'SMB';  
        files(i).name = files(i).name;
    % (b) DFW
    elseif contains(files(i).name,'SMB0') && contains (files(i).name,'TF0_')
        files(i).change = str2double(files(i).name(regexp(files(i).name,'W')+1:...
            regexp(files(i).name,'m_')-1)); % m 
        files(i).changeIn = 'DFW';  
        files(i).name = files(i).name;
    % (c) F_T
    elseif contains(files(i).name,'SMB0') && contains (files(i).name,'DFW0m')        
        files(i).change = str2double(files(i).name(regexp(files(i).name,'TF')+2:...
            regexp(files(i).name,'_geom')-1)); % m 
        files(i).changeIn = 'TF';  
        files(i).name = files(i).name;        
    end
end
% sort by changeIn and change
files = struct2table(files);
files = sortrows(files,[8,7]);
% save indices for SMB & SMR
ISMB = flipud(find(strcmp('SMB',table2array(files(:,8)))));
IDFW = find(strcmp('DFW',table2array(files(:,8))));
IFT = find(strcmp('TF',table2array(files(:,8))));
files = table2struct(files);
% DFW
for i=1:length(IDFW)
    file = load(files(IDFW(i)).name);
    DFW(i).change_DFW = file.DFW2;
    DFW(i).change_DFW_units = "m";
    DFW(i).H = interp1(file.x2, file.H2, x);
    DFW(i).h = interp1(file.x2, file.h2, x);
    DFW(i).Hh_units = "m";
    DFW(i).U = interp1(file.x2, file.U2, x) .* 3.1536e7;
    DFW(i).U_units = "m/yr";
    DFW(i).x_cf_final = x(dsearchn(x', file.x2(file.c2)));
    DFW(i).x_gl_final = x(dsearchn(x', file.x2(file.gl2)));
    DFW(i).t = t; 
    DFW(i).t_units = "deciyear";
    DFW(i).Q_gl = file.Fgl2; 
    DFW(i).Q_gl_units = "Gt/yr";
    DFW(i).x_cf = file.XCF2;
    DFW(i).x_gl = file.XGL2;
end

% -----SMB
for i=1:length(ISMB)
    file = load(files(ISMB(i)).name);
    SMB(i).change_SMB = files(ISMB(i)).change;
    SMB(i).change_SMB_units = "m/yr";
    SMB(i).H = interp1(file.x2, file.H2, x);
    SMB(i).h = interp1(file.x2, file.h2, x);
    SMB(i).Hh_units = "m";
    SMB(i).U = interp1(file.x2, file.U2, x) .* 3.1536e7;
    SMB(i).U_units = "m/yr";
    SMB(i).x_cf_final = x(dsearchn(x', file.x2(file.c2)));
    SMB(i).x_gl_final = x(dsearchn(x', file.x2(file.gl2)));
    SMB(i).t = t; 
    SMB(i).t_units = "deciyear";
    SMB(i).Q_gl = file.Fgl2; 
    SMB(i).Q_gl_units = "Gt/yr";
    SMB(i).x_cf = file.XCF2;
    SMB(i).x_gl = file.XGL2;
end

% -----F_T
for i=1:length(IFT)
    file = load(files(IFT(i)).name);
    FT(i).change_FT = files(IFT(i)).change;
    FT(i).H = interp1(file.x2, file.H2, x);
    FT(i).h = interp1(file.x2, file.h2, x);
    FT(i).Hh_units = "degrees_C";
    FT(i).U = interp1(file.x2, file.U2, x) .* 3.1536e7;
    FT(i).U_units = "m/yr";
    FT(i).x_cf_final = x(dsearchn(x', file.x2(file.c2)));
    FT(i).x_gl_final = x(dsearchn(x', file.x2(file.gl2)));
    FT(i).t = t; 
    FT(i).t_units = "deciyear";
    FT(i).Q_gl = file.Fgl2; 
    FT(i).Q_gl_units = "Gt/yr";
    FT(i).x_cf = file.XCF2;
    FT(i).x_gl = file.XGL2;
end

% -----SMB_enh
cd([results_path, '2_SMB_enh/']);
files = dir('*.mat');
for i=1:length(files)
    files(i).change = str2double(files(i).name(regexp(files(i).name,'B')+1:...
            regexp(files(i).name,'_enh')-1)); 
end
% sort by changeIn and change
files = struct2table(files);
files = sortrows(files,7, 'descend');
files = table2struct(files);
for i=1:length(files)
    file = load(files(i).name);
    SMB_enh(i).change_SMB_enh = files(i).change;
    SMB_enh(i).change_SMB_enh_units = "m/yr";
    SMB_enh(i).H = interp1(file.x2, file.H2, x);
    SMB_enh(i).h = interp1(file.x2, file.h2, x);
    SMB_enh(i).Hh_units = "m";
    SMB_enh(i).U = interp1(file.x2, file.U2, x) .* 3.1536e7;
    SMB_enh(i).U_units = "m/yr";
    SMB_enh(i).x_cf_final = x(dsearchn(x', file.x2(file.c2)));
    SMB_enh(i).x_gl_final = x(dsearchn(x', file.x2(file.gl2)));
    SMB_enh(i).t = t; 
    SMB_enh(i).t_units = "deciyear";
    SMB_enh(i).Q_gl = file.Fgl2; 
    SMB_enh(i).Q_gl_units = "Gt/yr";
    SMB_enh(i).x_cf = file.XCF2;
    SMB_enh(i).x_gl = file.XGL2;
end

% -----SMB_enh & FT
cd([results_path, '3_SMB_enh+TF/']);
files = dir('*.mat');
for i=1:length(files)
    files(i).change_SMB_enh = str2double(files(i).name(regexp(files(i).name,'B')+1:...
            regexp(files(i).name,'_enh')-1)); 
    files(i).change_FT = str2double(files(i).name(regexp(files(i).name,'F')+1:...
            regexp(files(i).name,'_geom')-1));
end
% sort by changeIn and change
files = struct2table(files);
files = sortrows(files,7, 'descend');
files = table2struct(files);
for i=1:length(files)
    file = load(files(i).name);
    SMB_enh_FT(i).change_SMB_enh = files(i).change_SMB_enh;
    SMB_enh_FT(i).change_SMB_enh_units = "m/yr";
    SMB_enh_FT(i).change_FT = files(i).change_FT;
    SMB_enh_FT(i).change_FT_units = "degrees_C";
    SMB_enh_FT(i).H = interp1(file.x2, file.H2, x);
    SMB_enh_FT(i).h = interp1(file.x2, file.h2, x);
    SMB_enh_FT(i).Hh_units = "m";
    SMB_enh_FT(i).U = interp1(file.x2, file.U2, x) .* 3.1536e7;
    SMB_enh_FT(i).U_units = "m/yr";
    SMB_enh_FT(i).x_cf_final = x(dsearchn(x', file.x2(file.c2)));
    SMB_enh_FT(i).x_gl_final = x(dsearchn(x', file.x2(file.gl2)));
    SMB_enh_FT(i).t = t; 
    SMB_enh_FT(i).t_units = "deciyear";
    SMB_enh_FT(i).Q_gl = file.Fgl2; 
    SMB_enh_FT(i).Q_gl_units = "Gt/yr";
    SMB_enh_FT(i).x_cf = file.XCF2;
    SMB_enh_FT(i).x_gl = file.XGL2;
end

% -----Save results
save([out_path, 'modeling_results_unperturbed.mat'], 'DFW');
save([out_path, 'modeling_results_SMB.mat'], 'SMB');
save([out_path, 'modeling_results_FT.mat'], 'FT');
save([out_path, 'modeling_results_SMB_enh.mat'], 'SMB_enh');
save([out_path, 'modeling_results_SMB_enh_FT.mat'], 'SMB_enh_FT');
cd(out_path);
disp('modeling results saved to file');
