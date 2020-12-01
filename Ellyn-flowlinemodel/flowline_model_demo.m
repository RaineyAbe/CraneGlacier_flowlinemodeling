function flowline_model_demo
%flowline_model_demo Basic flowline model set-up to execute the
%   demonstration (*demo*) files included with the script. The code is
%   explained in detail in the flowline_model_demo_userguide.pdf file.
%
%   The input files and parameterizations can be modified to reproduce
%   observational data.
%
%   Copyright 2020 Ellyn M. Enderlin, ellynenderlin@boisestate.edu
%   $Version: 1.1 $  $Date: 30-March-2020 11:33:00 $

%% Define time- and space-independent variables
warning off; %turn off warnings (velocity coefficient matrix is close to singular)

%clear all existing data and close all figures
close all; clear all; 

%Domain length and desired grid spacing
L = 60000; %domain length (m)
dx0 = 100; %desired grid spacing (m)

%Ice, ocean water, and fresh water densities and gravitational acceleration
rho_i = 917; %ice density (kg m^-3)
rho_sw = 1028; %ocean water density (kg m^-3)
rho_fw = 1000; %fresh water density (kg m^-3)
g = 9.8; %acceleration (m s^-2)

%Stress parameters
m = 3; %basal sliding exponent
beta = 2e+04; %basal roughness factor ((m/s)^(-1/m)) 
n = 3; %flow law exponent *NOTE: cannot use nthroot for Rxx calculation if n is even
A_min = 3.5e-25; %minimum rate factor value (Pa^-n s^-1)
A_max = 9.3e-25; %maximum rate factor value (Pa^-n s^-1)
E = 1.0; %enhancement factor (accounts for fabric development)

%Climate-related parameters
ELA = 800; %equilibrium line altitude (m)
fwd = 4; %fresh water depth in crevasses (m)

%Maximum thickness cut-off to check for instability
H_max = 1500; %maximum thickness (m)

%Time step (s)
dt = 315360; %make sure this satisfies the CFL condition (dt <=dx/U)

%% Adjust variables to a set starting grid & set-up the time counters
%load the starting file
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/Ellyn-flowlinemodel/';
cd(homepath);
load flowline_model_demo_initialization.mat; 
%x = distance from the ice divide (m), H = thickness (m), U = speed (m s^-1),
%dUdx = strain rate (s^-1), dx = grid spacing (m) 

%regrid the initialization data using the interp1 function to match the desired grid spacing
x_initial = x; %rename initialization distance vector
x = 0:dx0:L; %desired distance vector (m from ice divide)
H = interp1(x_initial,H,x,'linear','extrap');
U = interp1(x_initial,U,x,'linear','extrap');
dUdx = gradient(U,x);

%set the initial grid spacing to the desired grid spacing
dx = dx0;

%get the bed (hb) and width (w) profiles using the desired grid spacing
[hb,W] = flowline_model_demo_geometry(L,dx,homepath); %zb = bed elevation (m a.s.l.), w = width (m)

%set the timers
year = 1; %keeps track of the number of model years 
year_previous = 0; %used to plot profiles once per model year
year_start = 0; %tracks the time elapsed
year_end = year_start+5; %stops time for the model

%set-up empty plots for the data
gl = find((-(rho_sw./rho_i).*hb)-H>0,1,'first'); %find the initial grounding line
c = find(H<20,1,'first'); %set the initial calving front to a default minimum ice thickness value
%set-up empty figures for the plots
fig2 = figure;
xlabel('Distance from ice divide (m)','fontsize',16); ylabel('Speed (m d^{-1})','fontsize',16); set(fig2,'position',[500 100 500 400]);%velocity profiles
fig1 = figure; xlabel('Distance from ice divide (m)','fontsize',16); ylabel('Elevation (m)','fontsize',16); set(fig1,'position',[0 100 500 400]); %elevation profiles
fig4 = figure; set(fig4,'position',[1000 100 500 400]); %grounding line time series
subplot(2,1,1); xlabel('Time (yrs)','fontsize',16); ylabel('x_{grounding line} (m from divide)','fontsize',16);
subplot(2,1,2); xlabel('Time (yrs)','fontsize',16); ylabel('U_{grounding line} (m d^{-1})','fontsize',16);
[x_gl,H_gl,U_gl] = flowline_model_demo_plots(year_start,year_end,year,x,hb,U,gl,c,rho_i,rho_sw,H,fig1,fig2,fig4);

%% Run the flowline model 
i=1;
while i
    %calculate the thickness required to remain grounded at each grid cell
    Hf = -(rho_sw./rho_i).*hb; %flotation thickness (m)
    
    %find the location of the grounding line and the end of the ice-covered domain
    gl = find(Hf-H>0,1,'first')-1; %grounding line location
    ice_end = (find(H<=0,1,'first')); %end of ice-covered domain
 
    %calculate the glacier's surface elevation and slope
    h = hb+H; %h = surface elevation (m a.s.l.)
    h(gl:length(x)) = (1-rho_i/rho_sw).*H(gl:length(x)); %adjust the surface elevation of ungrounded ice to account for buoyancy
    dhdx = gradient(h,x); %surface slope (unitless)
    
    %increase the rate factor along flow (default is a linear increase)
    j=1;
    for j = 1:length(x);
        %assume that the ice temperature (& rate factor) are a function of
        %distance from the divide and are temporally-fixed
        A(j) = A_min + (A_max-A_min).*(x(j)/x(end)); 
    end
    
    %find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    Rxx = 2*nthroot((dUdx./(E*A)),n); %resistive stress (Pa)
    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); %crevasse penetration depth (m)
    c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); %calving front located where the inland-most crevasse intersects sea level
    %if the crevasses never intersect sea level
    if isempty(c) == 1;
        c = find(H<20,1,'first'); %set the calving front to a default minimum ice thickness value
    end
    %if the crevasses first interset sea level inland of the grounding line 
    if c <= gl;
        c = gl; %set the grounding line as the calving front
    end
    
    %calculate the effective pressure (ice overburden pressure minus water
    %pressure) assuming an easy & open connection between the ocean and
    %ice-bed interface
    sl = find(hb<=0,1,'first'); %find where the glacier base first drops below sea level
    N_ground = rho_i*g*H(1:sl); %effective pressure where the bed is above sea level (Pa)
    N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); %effective pressure where the bed is below sea level (Pa)
    N = [N_ground N_marine];
    N(N<0)=0; %cannot have negative values
        
    %calculate the thickness on the staggered grid for use in stress balance equations
    Hm = (H(1:end-1) + H(2:end))./2; %use forward differences
    Hm = [Hm,0];
    Hm(Hm<0) = 0; %cannot have negative values
    
    %solve the stress balance equations to obtain speed values (U)
    [U,dUdx,vm] = U_convergence(x,U,dUdx,Hm,H,A,E,N,W,dhdx,dx,c,ice_end,n,m,beta,rho_i,rho_sw,g,year,year_end,Rxx);
    
    %calculate the ice flux
    F = U.*H.*W; %ice flux (m^3 s^-1)
    F(isnan(F))=0;
    
    %calculate the  change in ice thickness from continuity
    dHdt = -(1./W).*gradient(F,x);
    dH = dHdt.*dt;
    
    %obtain values for the surface mass balance & submarine melting terms
    [SMB,submelt] = flowline_model_demo_SMB(x,h,ELA,gl,homepath); 
    %SMB = surface mass balance (m s^-1), submelt = submarine melt rate (m s^-1)
    
    %new thickness (change from dynamics, SMB, & submarine melting)
    Hn = H + dH + submelt.*dt + SMB.*dt; 
    Hn(Hn < 0) = 0; %remove any negative thickness values
    H = Hn; %set as the new thickness value
    
    %stop the model if it behaves unstably (monitored by ice thickness)
    if max(H) > H_max;
        disp(['Adjust dt']);
        cd(homepath); 
        save flowline_model_demo_max_thick.mat x hb W dx H h U dUdx A N fwd gl c dt year;

        return
    end
    
    %advance the model time
    year = year + dt./31536000;

    %find the new end of the ice-covered domain
    ice_end = (find(H<=0,1,'first'));
    
    %find the precise location of the grounding line (where H = Hf)
    xf = interp1((H(1:ice_end)-Hf(1:ice_end)),x(1:ice_end),0,'spline','extrap'); 
    
    %adjust the grid spacing so the grounding line is continuously tracked
    xl = round(xf/dx0); %number of ideal grid spaces needed to reach the grounding line
    dx = xf/xl; %new grid spacing (should be ~dx0)
    xn = 0:dx:L; %new distance vector
    
    %adjust geometry to the new distance vector
    [hb,W] = flowline_model_demo_geometry(L,dx,homepath); 
    
    %adjust the space-dependent variables to the new distance vector
    H = interp1(x,H,xn,'linear','extrap'); %ice thickness (m)
    Hf = interp1(x,Hf,xn,'linear','extrap');
    U = interp1(x,U,xn,'linear','extrap'); %speed (m s^-1)
    dUdx = interp1(x,dUdx,xn,'linear','extrap'); %strain rate (s^-1)
    A = interp1(x,A,xn,'linear','extrap'); %rate factor (Pa^-n s^-1)
    
    %find the location of thegrounding line and end of the ice-covered domain for the adjusted data
    gl = find(Hf-H>0,1,'first');
    ice_end = (find(H<=0,1,'first'));
    
    %rename the distance vector
    x = xn; %distance from the divide (m)
    
    %calculate the new surface elevation and slope
    h = hb+H; %grounded ice surface elevation (m a.s.l.)
    h(gl:length(x)) = (1-rho_i/rho_sw).*H(gl:length(x)); %floating ice surface elevation (m a.s.l.)
    
    %find the new calving front location for plotting purposes
    Rxx = 2*nthroot((dUdx./(E*A)),n); %resistive stress (Pa)
    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); %crevasse penetration depth (m)
    c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); %calving front located where the inland-most crevasse intersects sea level
    %if the crevasses never intersect sea level
    if isempty(c) == 1;
        c = find(H<20,1,'first'); %set the calving front to a default minimum ice thickness value
    end
    %if the crevasses first interset sea level inland of the grounding line 
    if c <= gl;
        c = gl; %set the grounding line as the calving front
    end

    %Plot the annual profiles & display select variables in the command window
    if year - year_previous >= 1;
        %plot the data & keep track on the grounding line position (gl_x),
        %thickness (gl_H), & speed (gl_U)
        [x_gl,H_gl,U_gl] = flowline_model_demo_plots(year_start,year_end,year,x,hb,U,gl,c,rho_i,rho_sw,H,fig1,fig2,fig4);
        
        %uncomment the following line(s) to display select data
%         disp(['Time (yr): ',num2str(year)]);
%         disp(['Grid size (m): ',num2str(dx)]);
%         disp(['Divide Surface Elevation (m): ', num2str(nanmean(h(1:5)))]);
%         disp(['GL: position= ',num2str(x(gl))]);
%         disp(['C: position= ',num2str(x(c))]);
%         disp(' '); %leave a space between annual data in command window
    end
    
    %stop the model at a set time & save
    if year > year_end;
        disp('model complete!');
        cd(homepath)
%         save flowline_model_demo_original.mat x hb W dx H h Hf U dUdx A E beta N fwd gl c ice_end dt year;
%         yrs = year_start:1:year_end;
%         save flowline_model_demo_gltimeseries_original.mat x_gl H_gl U_gl yrs;
        figure(fig1); saveas(gcf,['flowmodel-elev',...
            '-bedrough',num2str(round(beta,2,'significant')),'bedexp',num2str(m),'-visc',num2str(round(nanmean(vm),3,'significant')),...
            '-ELA',num2str(ELA),'-fwd',num2str(fwd),'.png'],'png');
        figure(fig2); saveas(gcf,['flowmodel-vel',...
            '-bedrough',num2str(round(beta,2,'significant')),'bedexp',num2str(m),'-visc',num2str(round(nanmean(vm),3,'significant')),...
            '-ELA',num2str(ELA),'-fwd',num2str(fwd),'.png'],'png');
        figure(fig4); saveas(gcf,['flowmodel-grounding',...
            '-bedrough',num2str(round(beta,2,'significant')),'bedexp',num2str(m),'-visc',num2str(round(nanmean(vm),3,'significant')),...
            '-ELA',num2str(ELA),'-fwd',num2str(fwd),'.png'],'png');
        break
    end
    
    %refresh the plot counter
    year_previous = floor(year);
    
    %loop through
    i=i+1;
    
end

end
