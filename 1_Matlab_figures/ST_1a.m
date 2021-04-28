%% 1. generate a figure form your data with code, as good as possible
A= peaks;
figure, 
t=annotation('textbox','String','a)', 'position',[0.01 0.99 0 0]); t.FontSize=14; t.LineStyle='None';
subplot(2,2,1), plot(A(1,:),'Displayname','B1 1µT'); axis([-Inf Inf -2 0.5]); xlabel('meas. index'); ylabel('signal S/S_0 [a.u.]');

t=annotation('textbox','String','b)', 'position',[0.47 0.99 0 0]); t.FontSize=14; t.LineStyle='None';
subplot(2,2,2), plot(A(2,:),'Displayname','B1 2µT'); axis([-Inf Inf -2 0.5]); xlabel('meas. index');ylabel('signal S/S_0 [a.u.]');

t=annotation('textbox','String','c)', 'position',[0.47 0.53 0 0]); t.FontSize=14; t.LineStyle='None';
subplot(2,2,3), plot(A(3,:)); axis([-Inf Inf -2 0.5]); xlabel('meas. index');ylabel('signal S/S_0 [a.u.]');

t=annotation('textbox','String','d)', 'position',[0.01 0.53 0 0]); t.FontSize=14; t.LineStyle='None';
subplot(2,2,4), imagesc(peaks); axis image; set(gca,'YTick',[]); set(gca,'XTick',[]); xlabel('image'); 
colorbar;

% this changes all fonts of axis to Garamond 11
child=get(gcf,'children');
for y = 1:length(child)
   set(child(y), 'fontname', 'Garamond','fontsize',12);  
end
% this changes all fonts of textboxes (lables a-d) to Garamond 16
child=findall(gcf,'Type','TextBox');
for y = 1:length(child)
   set(child(y), 'fontname', 'Garamond','fontsize',16);  
end
% this sets the size of the figure, bets to use this to adjust your final output size
set(gcf,'Position', [1200,350, 440,400]); % single column
% set(gcf,'Position', [1000,200, 840,480]); % double column

%% 2. rearrange your figure as you want

%% 3. use some codes to make subfigures same sizes again or color limites the same, etc.
% set certain size of axis ( first activate one axis in the figue by
% clicking on it
P=get(gca,'Position')
set(gca,'Position',[P(1),P(2),0.7,0.3]);

% set certain clims of image
get(gca,'Clim')
set(gca,'Clim',[-2 2])


%% 4.manipulate data 
% set interactive mode and click on one curve 
% then use  the following to extract X and Y data from the plot
x=get(gco,'XData')
y=get(gco,'YData')

% manipulate the data and plot again, copy this new line in the old plot
Y=y/max(abs(y));
figure, plot(x,Y);

% same for images, copy the image by clicking on it, STR-C, then click on
% the AXIS you want to paste in, (hitthe thin line at the edge of the
% current image or just delet the previosu image, then STR-V
C=get(gco,'CData')
c=C.*C;
figure, imagesc(c)



%% to resize plots it is sometimes useful to chaneg teh units from normalized to absolute
set(gca,'Units','centimeters');
set(gca,'Units','normalized');

%% set Latex interpreter as default
set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
set(groot, 'DefaultAxesFontName', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');

%% better subplot tools
https://www.mathworks.com/matlabcentral/fileexchange/20003-panel
https://de.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w






