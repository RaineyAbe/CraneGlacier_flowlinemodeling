function [hb,W] = flowline_model_demo_geometry(L,dx,homepath)
%flowline_model_demo_geometry 
%   Basic geometry used for the flowline model 
%   demonstration (*demo*). The code is explained in detail in the 
%   flowline_model_demo_userguide.doc file.
%
%   The geometry can can be modified directly in the flowline model or 
%   in a tab-delimited ascii text file.
%
%   Copyright 2013 Ellyn M. Enderlin, ellyn.enderlin@gmail.com
%   $Version: 1.0 $  $Date: 28-March-2013 15:39:00 $

%set-up the grid using the model grid spacing
x = 0:dx:L;

%read the ascii text file containing the model geometry data
%row 1 = header
%column 1 = distance, column 2 = bed elevation, column 3 = width
cd(homepath);
P = dlmread('flowline_model_demo_ascii_geometry.txt','\t',1,0);
x_input = P(:,1)';
hb_input = P(:,2)';
W_input = P(:,3)';

%interpolate the bed elevation and width profiles to match the model grid
hb = interp1(x_input,hb_input,x,'linear','extrap');
W = interp1(x_input,W_input,x,'linear','extrap');

end