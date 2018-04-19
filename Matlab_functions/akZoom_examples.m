%% a) Simple Plot
figure
plot(10:24,rand(1,15));
akZoom();

%% b) Image
figure
imagesc(magic(40));
akZoom();

%% c) Plotyy (linked axes)
figure
plotyy(1:15, rand(1,15), 1:15, rand(1,15));
akZoom();

%% d) Plotyy (independent axes)
figure
ax = plotyy(1:15, rand(1,15), 1:15, rand(1,15));
akZoom({ax(1), ax(2)});

%% e) Subplots (independent axes)
figure
for k = 1:4
  y = rand(1,15);
  subplot(2, 2, k);
  plot(y);
end
akZoom();

%% e) Subplots (linked axes)
figure
ax = NaN(4,1); % this will store the axes handles
for k = 1:4
  y = rand(1,15);
  ax(k) = subplot(2, 2, k);
  plot(y);
end
akZoom(ax);

%% f) Subplots (mixture of linked and indipendent axes)
figure
ax = NaN(4,1); % this will store the axes handles
for k = 1:4
  y = rand(1,15);
  ax(k) = subplot(2, 2, k);
  plot(y);
end
akZoom({[ax(1),ax(3)],ax(2),[ax(3),ax(4)]});

%% g) Different figures (linked)
% Note that you can simply call akZoom('all_linked') for linking all open
% figures.
figure;
imshow(imread('peppers.png'));
ax1 = gca;
figure
imshow(rgb2gray(imread('peppers.png')));
ax2 = gca;
akZoom([ax1,ax2]);

% move position of second figure a bit to make user aware of it:  
pos = get(gcf,'Position'); set(gcf, 'Position', [pos(1)+100, pos(2)-100, pos(3), pos(4)]);

%% h) All axes in all open figures (independent)
figure
plotyy(1:15, rand(1,15), 1:15, rand(1,15));
figure
plot(rand(1,15));
akZoom('all');

% move position of second figure a bit to make user aware of it:  
pos = get(gcf,'Position'); set(gcf, 'Position', [pos(1)+100, pos(2)-100, pos(3), pos(4)]);

%% i) All axes in all open figures (linked)
figure
plotyy(1:15, rand(1,15), 1:15, rand(1,15));
figure
semilogy(1:150, rand(1,150));
akZoom('all_linked');

% move position of second figure a bit to make user aware of it:  
pos = get(gcf,'Position'); set(gcf, 'Position', [pos(1)+100, pos(2)-100, pos(3), pos(4)]);

%% j) LogLog plot with out-of-bounds feedback
% Note the blinking boundaries in the plot when you try to pan/zoom in a 
% way that would make the limits exceed valid double values, 
% i.e. >1.8e308 or <2.2e-308 
figure
loglog(exp(-100:800), exp(-100:800).*rand(1,901).^2);
akZoom()

%% k) Test with patch objects
figure
xdata = [2 1; 1 8; 8 8];
ydata = [5 2; 8 4; 4 0];
patch(xdata,ydata,'r');
akZoom()

%% l) Change mapping of mouse buttons
figure
semilogx(rand(1,150));
akZoom('rlm')

%% m) Test for axes within ui-subpanel(s)
h = figure;
hp = uipanel('Title','Sub-panel', 'BackgroundColor','red', 'Position',[.1 .1 .8 .8]);
hsp = uipanel('Parent',hp,'Title','Sub-sub-panel','Position',[.3 .1 .6 .8]);
ax = axes('Parent', hsp);
plot(10:24,rand(1,15));
akZoom();

%% n) Demonstrate different MATLAB-behaviour for AspectRatio or "axis equal"
img = imread('peppers.png');

figure
imshow(img);
set(gcf, 'Position', [500,600,400,150]);
set(gcf,'Name', 'imshow()', 'NumberTitle', 'off')
akZoom

figure
set(gcf, 'Position', [500,350,400,150]);
image(img);
set(gca, 'DataAspectRatio', [1 1 1])
set(gca, 'visible', 'off')
set(gcf,'Name', 'DataAspectRatio [1 1 1] (same as imshow)', 'NumberTitle', 'off')
akZoom

figure
set(gcf, 'Position', [500,100,400,150]);
image(img);
axis equal
set(gca, 'visible', 'off')
set(gcf,'Name', 'axis equal', 'NumberTitle', 'off')
akZoom