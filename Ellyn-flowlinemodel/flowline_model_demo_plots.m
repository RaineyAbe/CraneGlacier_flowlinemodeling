function [x_gl,H_gl,U_gl] = flowline_model_demo_plots(start_year,end_year,year,x,hb,U,gl,c,rho_i,rho_sw,H,fig1,fig2,fig4)
%flowline_model_demo_plots 
%   Basic plotting script used for the flowline model demonstration (*demo*). 
%   Three plots will be generated at yearly model time steps using the code: 
%   elevation profile for the entire model domain (with a zoomed inset), 
%   speed profile for the entire model domain, and time series of the 
%   grounding line position and speed. The code is explained in detail in 
%   the flowline_model_demo_userguide.doc file. 
%
%   Copyright 2020 Ellyn M. Enderlin, ellynenderlin@boisestate.edu
%   $Version: 1.1 $  $Date: 30-March-2020 11:33:00 $

%create an annual colormap
cmap = colormap(parula(end_year-start_year+round(.1*(end_year-start_year))));

%define profile axes limits
xlimits = [0,55000]; xticks = [0:5000:55000]; xticklabels = [0:5:55];
ylimits_U = [0,0.5]; yticks_U = [0:0.1:0.5];
ylimits_h = [-200,1400]; yticks_h = [-200:200:1400];
xlimits_zoom = [40000,50000]; xticks_zoom = [40000:2000:50000]; xticklabels_zoom = [40:2:50];
ylimits_hzoom = [-150, 150]; yticks_hzoom = [-150:50:150];

%get the annual grounding line position, speed, and thickness
x_gl(round(year-start_year)) = x(gl);
H_gl(round(year-start_year)) = H(gl);
U_gl(round(year-start_year)) = U(gl);

%ice surface elevation
h = hb+H; %h = surface elevation (m a.s.l.)
h(gl:length(x)) = (1-rho_i/rho_sw).*H(gl:length(x)); %adjust the surface elevation of ungrounded ice to account for buoyancy

%plot the profiles 
if round(year) < end_year;
    figure(fig1);
%     clf(fig1); %uncomment to clear the figure
    
    %if the base of the terminus is below sea level, plot the ice to the calving face
    if hb(c) <= 0;
        %plot the bottom of the floating ice tongue
        zf = -(rho_i/rho_sw).*H; zf(gl) = hb(gl); %depth of the tongue below sea level (m)
        plot(x(gl:c),zf(gl:c),'color',cmap(round(year),:),'linewidth',2); hold on;

        %plot the calving face
        cx = [x(c) x(c)]; cy = [h(c) zf(c)];
        plot(cx,cy,'color',cmap(round(year),:),'linewidth',2); hold on;

        %plot the surface of the ice
        plot(x(1:c),h(1:c),'color',cmap(round(year),:),'linewidth',2); hold on; %exclude 'calved' ice from c:terminus that is needed for stability

        %plot the bed
        plot(x,hb,'k','linewidth',2); hold on;
    else
        plot(x,h,'color',cmap(round(year),:),'linewidth',2); hold on; %plot the entire grounded profile
        plot(x,hb,'k','linewidth',2); hold on; %plot the bed
    end
    plot(x,zeros(1,length(x)),'--','color',[0.5 0.5 0.5]); %plot sea level
    if round(year==1); 
        set(gca,'ylim',ylimits_h,'xlim',xlimits,'ytick',yticks_h,...
            'xtick',xticks,'xticklabel',xticklabels);
        xlabel('Distance from ice divide (km)','fontsize',12); ylabel('Elevation (m)','fontsize',12); 
    end
    set(fig1,'position',[0 100 500 400]); %set the position
    title(['year: ',num2str(round(year))],'fontsize',12); drawnow;

    %speed profile
    figure(fig2);
%     clf(fig2); %uncomment to clear the figure
    plot(x(1:c),U(1:c).*86400,'color',cmap(round(year),:),'linewidth',2); hold on;
    set(gca,'ylim',ylimits_U,'xlim',xlimits,'ytick',yticks_U,...
        'xtick',xticks,'xticklabel',xticklabels);
    xlabel('Distance from ice divide (km)','fontsize',12);
    ylabel('Speed (m d^{-1})','fontsize',12);
    set(fig2,'position',[500 100 500 400]); %set the position
%     annotation(gcf,'textbox',[0.3 0.85 0.2 0.05],'string',['year: ',num2str(round(year))],...
%         'backgroundcolor',[1 1 1],'edgecolor',[1 1 1],'fontsize',12);
    title(['year: ',num2str(round(year))],'fontsize',12); drawnow;

else %plot the last year's profiles in green
    figure(fig1);
%     clf(fig1); %uncomment to clear the figure
    
    %if the base of the terminus is below sea level, plot the ice to the calving face
    if hb(c) <= 0;
        %plot the bottom of the floating ice tongue
        zf = -(rho_i/rho_sw).*H; zf(gl) = hb(gl); %depth of the tongue below sea level (m)
        plot(x(gl:c),zf(gl:c),'g','linewidth',2); hold on;

        %plot the calving face
        cx = [x(c) x(c)]; cy = [h(c) zf(c)];
        plot(cx,cy,'g','linewidth',2); hold on;

        %plot the surface of the ice
        plot(x(1:c),h(1:c),'g','linewidth',2); hold on; %exclude 'calved' ice from c:terminus that is needed for stability

        %plot the bed
        plot(x,hb,'k','linewidth',2); hold on;
    else
        plot(x,h,'color',cmap(round(year),:),'linewidth',2); hold on; %plot the entire grounded profile
        plot(x,hb,'k','linewidth',2); hold on; %plot the bed
    end
    plot(x,zeros(1,length(x)),'--','color',[0.5 0.5 0.5]); %plot sea level
    set(gca,'ylim',ylimits_h,'xlim',xlimits,'ytick',yticks_h,...
        'xtick',xticks,'xticklabel',xticklabels);
    xlabel('Distance from ice divide (m)','fontsize',12);
    ylabel('Elevation (m)','fontsize',12);
    set(fig1,'position',[0 100 500 400]); %set the position
    title(['year: ',num2str(round(year))],'fontsize',12); drawnow;

    %speed profile
    figure(fig2);
%     clf(fig2); %uncomment to clear the figure
    plot(x(1:c),U(1:c).*86400,'g','linewidth',2); hold on;
    set(gca,'ylim',ylimits_U,'xlim',xlimits,'ytick',yticks_U,...
        'xtick',xticks,'xticklabel',xticklabels);
    xlabel('Distance from ice divide (m)','fontsize',12);
    ylabel('Speed (m d^{-1})','fontsize',12);
    set(fig2,'position',[500 100 500 400]); %set the position
%     annotation(gcf,'textbox',[0.3 0.85 0.2 0.05],'string',['year: ',num2str(round(year)-1)],...
%         'backgroundcolor',[1 1 1],'edgecolor',[1 1 1],'fontsize',12);
    title(['year: ',num2str(round(year))],'fontsize',12); drawnow;

end

%plot the profile for the zoomed inset in a dummy plot
fig3 = figure;
set(fig3,'position',[0 20 1 1]);
if hb(c) <= 0;
    %plot the bottom of the floating ice tongue
    zf = -(rho_i/rho_sw).*H; zf(gl) = hb(gl); %depth of the tongue below sea level (m)
    plot(x(gl:c),zf(gl:c),'color',cmap(round(year),:),'linewidth',2); hold on;
    
    %plot the calving face
    cx = [x(c) x(c)]; cy = [h(c) zf(c)];
    plot(cx,cy,'color',cmap(round(year),:),'linewidth',2); hold on;
    
    %plot the surface of the ice
    plot(x(1:c),h(1:c),'color',cmap(round(year),:),'linewidth',2); hold on; %exclude 'calved' ice from c:terminus that is needed for stability
    
    %plot the bed
    plot(x,hb,'k','linewidth',2); hold on;
else
    plot(x,[h',hb'],'linewidth',2); %plot the entire grounded profile & the bed
end
plot(x,zeros(1,length(x)),'--','color',[0.5 0.5 0.5]); %plot sea level
set(gca,'ylim',ylimits_hzoom,'xlim',xlimits_zoom,'ytick',yticks_hzoom,'yticklabel',[],...
        'xtick',xticks_zoom,'xticklabel',[]); %set the limits for the inset
if round(year==1);
    set(gca,'ylim',ylimits_hzoom,'xlim',xlimits_zoom,'ytick',yticks_hzoom,'yticklabel',yticks_hzoom,...
        'xtick',xticks_zoom,'xticklabel',xticklabels_zoom); %set the limits for the inset
end

%use the inset.m script to replot the zoom as an inset
[h_m h_i]=inset(fig1,fig3); %h_m is the handle for figure1, h_i is the handle for the inset
close(fig3); %close the dummy plot
figure(fig1); drawnow;
% annotation(gcf,'textbox',[0.3 0.85 0.2 0.05],'string',['year: ',num2str(round(year)-1)],...
%     'backgroundcolor',[1 1 1],'edgecolor',[1 1 1],'fontsize',12);

%time series of grounding line position & speed
figure(fig4);
subplot(2,1,1);
plot(round(year)-1,x_gl(end),'*','color',cmap(round(year),:)); hold on;
set(gca,'xlim',[0,end_year]); grid on;
xlabel('Time (yrs)','fontsize',12);
ylabel('x_{grounding line} (m from divide)','fontsize',12);

subplot(2,1,2);
plot(round(year)-1,U_gl(end).*86400,'*','color',cmap(round(year),:)); hold on;
set(gca,'xlim',[0,end_year]); grid on;
xlabel('Time (yrs)','fontsize',12);
ylabel('U_{grounding line} (m d^{-1})','fontsize',12);
set(fig4,'position',[1000 100 500 400]); %set the position
drawnow;

end