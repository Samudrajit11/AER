function [divarea1,divarea2,divarea3,divarea4] = calc_AER(tr,FPS)

%;------
%;Calculate AER for the given FPS, for 2 particles.
%; output is the nm^2/msec vector for XX and for YY.
%; tr: is a matrice of 2 particles x and y positions (:,x1 y1 x2 y2)
%; FPS: is the FPS
%;------

%Create 4 vectors of X and Y (x1x2, y1y2, x1y2, y1x2)
xx = [tr(:,1),tr(:,3)]; %pixels
yy = [tr(:,2),tr(:,4)]; %pixels
xy = [tr(:,1),tr(:,4)]; %pixels
yx = [tr(:,2),tr(:,3)]; %pixels
% changed xx for normal modes of masses on a spring
%%xx = [tr(:,1)-tr(:,3),(tr(:,1)+tr(:,3))].*58.2e-9;
%%yy = [tr(:,2)-tr(:,4),(tr(:,2)+tr(:,4))].*58.2e-9;


%Fit to ellipse (if sig>>0 can have some effect)
%xx_ellipse = fit_ellipse(xx(:,1),xx(:,2));
%yy_ellipse = fit_ellipse(yy(:,1),yy(:,2));
%xy_ellipse = fit_ellipse(xy(:,1),xy(:,2));
%yx_ellipse = fit_ellipse(yx(:,1),yx(:,2));
%Finds ellipse centers
%center_xx = [xx_ellipse.X0_in xx_ellipse.Y0_in];
%center_yy = [yy_ellipse.X0_in yy_ellipse.Y0_in];
%center_xy = [xy_ellipse.X0_in xy_ellipse.Y0_in];
%center_yx = [yx_ellipse.X0_in yx_ellipse.Y0_in];

center_xx = mean(xx(:,1:2));
center_yy = mean(yy(:,1:2));
center_xy = mean(xy(:,1:2));
center_yx = mean(yx(:,1:2));

%Rescales the vectors for mean 0 using the centers
xx = [xx-center_xx,zeros(length(xx),1)];
yy = [yy-center_yy,zeros(length(xx),1)];
xy = [xy-center_xy,zeros(length(xx),1)];
yx = [yx-center_yx,zeros(length(xx),1)];

%Find the area of a triangle between 2 measurements using cross product
c1 = cross(xx(1:end-1,:),xx(2:end,:));
c2 = cross(yy(1:end-1,:),yy(2:end,:));
c3 = cross(xy(1:end-1,:),xy(2:end,:));
c4 = cross(yx(1:end-1,:),yx(2:end,:));

area1 = 0.5*sum(c1,2); area2 = 0.5*sum(c2,2); area3 = 0.5*sum(c3,2); area4 = 0.5*sum(c4,2);

%Adds a time vector, divides area by time
sumarea1 = cumsum([area1,ones(length(area1),1)],1);
sumarea2 = cumsum([area2,ones(length(area2),1)],1);
sumarea3 = cumsum([area3,ones(length(area3),1)],1);
sumarea4 = cumsum([area4,ones(length(area4),1)],1);

divarea1 = FPS*58.2^2*10^-3*(sumarea1(:,1)./(sumarea1(:,2))); %get pixel^2/frame and transform into nm^2./msec
divarea2 = FPS*58.2^2*10^-3*(sumarea2(:,1)./(sumarea2(:,2)));
divarea3 = FPS*58.2^2*10^-3*(sumarea3(:,1)./(sumarea3(:,2)));
divarea4 = FPS*58.2^2*10^-3*(sumarea4(:,1)./(sumarea4(:,2)));
end