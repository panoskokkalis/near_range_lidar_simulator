function [overlap, DFO, RFOV, Atilt_max]=calculate_overlap(alt, lamda, bexp, dred, FT, DT, Dfieldstop, DTL, Atilt, TFOV, DL, height)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%caclculations according to Kamil Stelmaszczyk et al., 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alt=[0:7.5:10000]';%1.375e-01
% [overlap_far, DFO_far, RFOV_far, Atilt_max_far]=calculate_overlap_far_range(alt, 355, 10, 8, 500, 200, 0.1, 170, 0.0687, 0.25, 8, 3200);
% [overlap_near, DFO_near, RFOV_near, Atilt_max_near]=calculate_overlap_far_range(alt, 355, 10, 8, 250, 80, 0.1, 110, 0.1687, 0.25, 8, 3200);
% [overlap_depol, DFO_depol, RFOV_depol, Atilt_max_depol]=calculate_overlap_far_range(alt, 355, 10, 8, 500, 70, 0.1, 110, 0.02, 0.25, 8, 3200);
% 
% figure; plot(alt, overlap_far, alt, overlap_near, alt, overlap_depol); myfig; xLabel('Height [m, a.s.l.]'); yLabel('Overlap Function'); legend('Far Range', 'Near Range', 'Depolarization'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

igauss=3;

% Convert to m
FT=FT.*10^-3;
DT=DT.*10^-3;
Dfieldstop=Dfieldstop.*10^-3;
DTL=DTL.*10^-3;
Atilt=Atilt.*10^-3;
TFOV=TFOV.*10^-3;
DL=DL.*10^-3;
lamda=lamda.*10^-9;

% calculation of rfov
RFOV=Dfieldstop./(2.*FT); % RFOV is exported in rad

% Apply the expansion./reduction functions
if (igauss >= 2) && (bexp > 0)
    TFOV_ex = (TFOV./dred).* 0.5 .* igauss;
    DL_ex = DL.*bexp.* 0.5 .* igauss;
else
    DL_ex = DL;
end

% calculation of DFO as value
DFO=(2.*DTL+DT+DL_ex)./(2.*(RFOV-TFOV_ex+Atilt));

% Check for constraints due to Atilt
if (Atilt > RFOV-TFOV_ex);
    disp('error : Your system will not be aligned in far range');
else
    disp('ok');  
end   

Atilt_max=RFOV-TFOV_ex; % Atilt_max is exported in rad

% Check for the optimum fibers diameter
NA_fiber=0.23;
fnum_ideal=0.5.*sqrt((1./NA_fiber)^2-1);
fnum_used=FT./DT;
div_from_ideal=100.*(fnum_ideal-fnum_used)./fnum_ideal;

% make_circle_points(FT, TFOV_ex, DL_ex, DT, height, DTL, Dfieldstop, RFOV, Atilt)

% OverLap Calculation No abbeartion is incuded
% e1=(2.*TFOV_ex+(DL_ex)./alt).*(FT+abber); % The size of the image
e1=FT.*(2.*TFOV_ex+(DL_ex+DT)./alt); % The size of the image

v=FT.*(DTL-Atilt.*alt)./alt; % The dispacement from the optical axis in image space
y1=2.*acos((Dfieldstop^2+4.*v.^2-e1.^2)./(4.*v.*Dfieldstop));
y2=2.*acos((e1.^2+4.*v.^2-Dfieldstop^2.)./(4.*v.*e1));

overlap=zeros(length(alt),1);

for i=1:length(alt)
    if v(i)>(Dfieldstop + e1(i))/2;
        overlap(i)=0;
    elseif((abs(Dfieldstop-e1(i)))/2 < v(i) &  v(i) < (Dfieldstop+e1(i))/2);
        overlap(i)=((y1(i)-sin(y1(i))).*Dfieldstop^2+(y2(i)-sin(y2(i))).*e1(i)^2)./(2.*pi.*e1(i)^2);
    elseif (v(i)<(e1(i)-Dfieldstop)/2 & e1(i) > Dfieldstop);
        overlap(i)=Dfieldstop^2./e1(i)^2;
    elseif (v(i)<=(Dfieldstop-e1(i))/2 & e1(i)<=Dfieldstop);
        overlap(i)=1;   
    end
end


% v=abs(FT.*(DTL-Atilt.*alt)./alt); % The dispacement from the optical axis in image space
% y1=2.*acos((Dfieldstop^2+4.*v.^2-e1.^2)./(4.*v.*Dfieldstop));
% y2=2.*acos((e1.^2+4.*v.^2-Dfieldstop^2.)./(4.*v.*e1));
% 
% overlap=zeros(length(alt),1);
% 
% for i=1:length(alt)
%     if v(i)>(Dfieldstop + e1(i))/2;
%         overlap(i)=0;
%     elseif(((abs(Dfieldstop-e1(i)))/2 <v(i)) && (v(i)<(Dfieldstop+e1(i))/2));
%         overlap(i)=((y1(i)-sin(y1(i)))*Dfieldstop^2+(y2(i)-sin(y2(i)))*e1(i)^2)/(2*pi*e1(i)^2);
%     elseif ((v(i)<(e1(i)-Dfieldstop)/2) && (e1(i) > Dfieldstop));
%         overlap(i)=Dfieldstop^2/e1(i)^2;
%     elseif ((v(i)<(Dfieldstop-e1(i))/2) && (e1(i)<=Dfieldstop));
%         overlap(i)=1;   
%     end
% end

% figure; plot(alt, overlap); myfig; xLabel('Height [m, a.s.l.]'); yLabel('Overlap Function'); legend('Far Range');

end

function make_circle_points(FT, TFOV_ex, DL_ex, DT, height, DTL, Dfieldstop, RFOV, Atilt)

% OBJECT SPACE
r_laser_object=TFOV_ex.*height; % circle radius
x_laser_object=(DL_ex./2)+(DT./2)+DTL-Atilt.*height;
y_laser_object=0;

% if (DTL-Atilt*height)==0 % the cross point is the DTL/Atilt
%     x_laser_object=0;
%     y_laser_image=0;
% elseif (DTL-Atilt*height)<0
%     x_laser_object=-(Atilt.*height-DTL-(DL_ex./2)-(DT./2));
%     y_laser_image=(FT.*(x_laser_object)./height);
% elseif (DTL-Atilt*height)>0
%     x_laser_object=(DL_ex./2)+(DT./2)+DTL-Atilt.*height;
%     y_laser_image=(FT.*(DTL-Atilt.*height)./height);
% end

x_telescope_object=0;
y_telescope_object=0;
r_telescope_object=RFOV.*height; % circle radius


% IMAGE SPACE
r_laser_image=0.5.*FT.*(2.*TFOV_ex+(DL_ex+DT)./height); % circle radius
y_laser_image=(FT.*(DTL-Atilt.*height)./height);
x_laser_image=0;

x_telescope_image=0;
y_telescope_image=0;
r_telescope_image=Dfieldstop./2; % circle radius

hold on
th = 0:pi/500:2*pi;

xunit_laser_object=r_laser_object*cos(th) + x_laser_object;
yunit_laser_object=r_laser_object*sin(th) + y_laser_object;
xunit_telescope_object=r_telescope_object*cos(th) + x_telescope_object;
yunit_telescope_object=r_telescope_object*sin(th) + y_telescope_object;

xunit_laser_image=r_laser_image*cos(th) + x_laser_image;
yunit_laser_image=r_laser_image*sin(th) + y_laser_image;
xunit_telescope_image=r_telescope_image*cos(th) + x_telescope_image;
yunit_telescope_image=r_telescope_image*sin(th) + y_telescope_image;

figure; plot(xunit_telescope_object, yunit_telescope_object, xunit_laser_object, yunit_laser_object); myfig; xlabel('Distance Telescope-Laser axis [m]'); ylabel('Object size [m]'); legend('Telescope', 'Laser'); title({'Cross sections of Telescope and Laser at 2.258 km height'})
figure; plot(xunit_telescope_image*1000, yunit_telescope_image*1000, xunit_laser_image*1000, yunit_laser_image*1000); myfig; xlabel('Image size [mm]'); ylabel('Displacement on optical axis [mm]'); legend('Dfieldstop', 'Image Size'); title({'Cross sections of field stop and image size from 2.258 km height'})


% figure; plot(xunit_telescope, yunit_telescope); myfig

hold off
end