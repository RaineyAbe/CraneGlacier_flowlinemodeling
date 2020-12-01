function [main_fig, h_inset]=inset(main_handle, inset_handle,inset_size)

% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.4
% The outputs are the axes-handles of both.
% 
% An examle can found in the file: inset_example.m
% 
% Moshe Lindner, August 2010 (C).

if nargin==2
    inset_size=0.4;
end

inset_size=inset_size*.7;
% figure
new_fig=main_handle;
main_fig = findobj(main_handle,'Type','axes');
% h_main = copyobj(main_fig(end),new_fig);
% set(h_main(end),'Position',get(main_fig(end),'Position'))
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
ax=get(main_fig(end),'Position');
set(h_inset,'Position', [.7*ax(1)+ax(3)-inset_size .5*ax(2)+ax(4)-inset_size inset_size 1.1*inset_size])
