function [SMB,submelt] = flowline_model_demo_SMB(x,h,ELA,gl,homepath)
%   flowline_model_demo_SMB 
%   Basic surface mass balance and submarine melting parameterizations used
%   for the flowline model demonstration (*demo*). The code is explained 
%   in detail in the flowline_model_demo_userguide.doc file.
%
%   The SMB and submelt profiles can can be modified directly in the flowline model or 
%   in a tab-delimited ascii text file .
%
%   Copyright 2013 Ellyn M. Enderlin, ellyn.enderlin@gmail.com
%   $Version: 1.0 $  $Date: 28-March-2013 15:38:00 $

    %calculate surface mass balance as a function of elevation with respect
    %to the equilibrium line altitude (ELA)
    x_ELA=find(h>=ELA,1,'last'); %find the location of the ELA
    if isempty(x_ELA) == 1;
        disp('warning: poor ELA choice');
    end

    %SMB profile (default is a piecewise linear relationship)
    c1 = 0.95; c2 = 0.3; %slope in the accumulation (c1) & ablation (c2) zones
    lapserate = -0.0071; %air temperature lapse rate (0.71^{o}C/100m)
    SMB = zeros(1,length(x)); %set-up an empty matrix of the correct size to fill-in with values
    SMB(1:x_ELA) = c1.*lapserate.*(h(x_ELA)-h(1:x_ELA))./31536000; %accumulation balance gradient (yr^-1) = c1*0.0071
    SMB(x_ELA+1:length(x)) = c2.*lapserate.*(h(x_ELA)-h(x_ELA+1:length(x)))./31536000; %ablation balance gradient (yr^-1) = c2*0.0071
     
    %submarine melt rate under floating tongue
    cd(homepath);
    P = dlmread('flowline_model_demo_ascii_submelt.txt','\t',1,0); %read the ascii text file containing the submarine melt rate data
    x_input = P(:,1)';
    submelt_input = P(:,2)';

    %construct the submarine melt rate profile
    submelt_floating = interp1(x_input,submelt_input,x,'linear','extrap'); %submarine melt rate (m s^-1)
    submelt_inland = zeros(1,gl-1); %set-up an empty matrix of the correct size to fill-in with values where the ice is floating
    submelt = -[submelt_inland, submelt_floating(1:length(x)-(gl-1))]; %concatenate 
    
end