function[outmat12]=output_aer_2part_fps(frame_full)

parttotal=2;
outmat12=[];

for frameind=1:2000
    frame=[];
    frame=frame_full(1:frameind:end,:);
    FPS=10000;
    FPS=FPS/frameind;      
    tr = sim_to_tr(frame,parttotal);
    st = AER_mult(tr,parttotal,FPS);
    out12=st.particle_12(end,1);
    
    %out1=(drivemat(kk)*0.01);

    outmat12=[outmat12;FPS out12];
    
end 

end

%......
function tr = sim_to_tr(frame,particle_num)
%reshapes frame for use later
L = length(frame);
tr = zeros(L*particle_num,6);
frame2 = [frame(:,1:2:end), frame(:,2:2:end)];
tr(:,1:2) = reshape(frame2,[L*particle_num,2])*10^6; %1e6 converts to um
tr(:,3:4) = ones(L*particle_num,2);

tr(:,5) = repmat(cumsum(ones(L,1)),particle_num,1);
tr_6 = ones(L,particle_num);
for i = 1:particle_num
    tr_6(:,i) = tr_6(:,i)*i;
end
tr(:,6) = tr_6(:);
end

%.....
function st = AER_mult(tr,NOF,FPS)

%;------
%; Calculate AER for the given FPS, for all particle combinations.
%; output is the a struct, containing (for every pair),
%; a [nm^2/msec] vector for XX and for YY.
%; tr: is just the regular tr for all of the particles.
%; FPS: is the FPS
%; NOF: number of particles in the tr.
%;------


%Create a matrice of all possible combinations of 2 particles
C = nchoosek(1:NOF,2);

for i=1:size(C,1)
    
    %Create the relevant tr for the AER (X,Y for the relevant particles)
    m_a = find(tr(:,6)==C(i,1)); m_b = find(tr(:,6)==C(i,2));
    tr_a = tr(m_a,1:2); tr_b = tr(m_b,1:2);
    tr2 = [tr_a tr_b];
    
    %Run the AER
    [XX,YY,XY,YX] = calc_AER(tr2,FPS);
    T = cumsum(ones(length(XX),1)).*1/FPS;
    %Create a fieldname according to the particles
    A = char(strjoin(string(C(i,:))));
    fieldname = A(~isspace(A));
    %saves into a struct, which will be outputed
    st.('particle_'+string(fieldname)) = [XX YY XY YX T];   
end


end

%....
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

divarea1 = FPS*10^3*(sumarea1(:,1)./(sumarea1(:,2))); % transform into nm^2./msec
divarea2 = FPS*10^3*(sumarea2(:,1)./(sumarea2(:,2)));
divarea3 = FPS*10^3*(sumarea3(:,1)./(sumarea3(:,2)));
divarea4 = FPS*10^3*(sumarea4(:,1)./(sumarea4(:,2)));
end