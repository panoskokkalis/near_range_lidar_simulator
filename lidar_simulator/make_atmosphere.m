function [export_aaer, export_baer, export_lraer, export_daer, baer_parallel, baer_cross, export_AOD]=make_atmosphere(height, alt, AOD_1, AOD_2, sm1, sm2, export, aerosol_type)

% height; the pertubation table height first layer height (1,1) up to
% height (1,2), second-top layer layer height (2,1) up to height (2,2) 
% AOD_1; the AOD at 532 of the bottom layer
% AOD_2; the AOD at 532 of the top layer
% sm1; the smoothing order of bottom layer a value 0-1 most good 0.008
% sm2; the smoothing order of top layer a value 0-1 most good 0.008
% k 0 or 1; 0 the steping finction profiles are exported while for k=1 the
% filtered 

% the bottom layer 
% load('aerosol_db\aerosol_type.mat')
load(['aerosol_db\' aerosol_type{1,1} '.mat'])

A_a_355_1=atm(1,2);A_b_355_1=atm(1,3); lr_355_1=atm(1,4); d_355_1=atm(1,5);
a_532_1=atm(2,2);b_532_1=atm(2,3); lr_532_1=atm(2,4); d_532_1=atm(2,5);
A_a_1064_1=atm(3,2);A_b_1064_1=atm(3,3); lr_1064_1=atm(3,4); d_1064_1=atm(3,5);

% the lofted layer 
% load('aerosol_db\Dust.mat')
load(['aerosol_db\' aerosol_type{2,1} '.mat'])

A_a_355_2=atm(1,2);A_b_355_2=atm(1,3); lr_355_2=atm(1,4); d_355_2=atm(1,5);
a_532_2=atm(2,2);b_532_2=atm(2,3); lr_532_2=atm(2,4); d_532_2=atm(2,5);
A_a_1064_2=atm(3,2);A_b_1064_2=atm(3,3); lr_1064_2=atm(3,4); d_1064_2=atm(3,5);

a_532_1=AOD_1/((height(1,2)-height(1,1)).*10^-3); % The extinction at 532 /km for the type and the AOD user defined
a_532_2=AOD_2/((height(2,2)-height(2,1)).*10^-3); % The extinction at 532 /km for the type and the AOD user defined

b_532_1=a_532_1/lr_532_1; % The backscatter  at 532 /km sr for the type and the LR user defined
b_532_2=a_532_2/lr_532_2; % The backscatter  at 532 /km sr for the type and the LR user defined

a_355_1=a_532_1*(350/550)^(-A_a_355_1); % The extinction at 355 for the first layer according to Aa exponent
a_355_2=a_532_2*(350/550)^(-A_a_355_2); % The extinction at 355 for the second layer according to Aa exponent

b_355_1=b_532_1*(350/550)^(-A_b_355_1); % The backscatter at 355 for the first layer according to Ab exponent
b_355_2=b_532_2*(350/550)^(-A_b_355_2); % The backscatter at 355 for the second layer according to Ab exponent

a_1064_1=a_532_1*(550/1000)^(A_a_1064_1); % The extinction at 1064 for the first layer according to Aa exponent
a_1064_2=a_532_2*(550/1000)^(A_a_1064_2); % The extinction at 1064 for the second layer according to Aa exponent

b_1064_1=b_532_1*(550/1000)^(A_b_1064_1); % The backscatter at 1064 for the first layer according to Ab exponent
b_1064_2=b_532_2*(550/1000)^(A_b_1064_2); % The backscatter at 1064 for the second layer according to Ab exponent

% array preallocation for time resuming 
baer_prof_355=zeros(size(alt));
LR_prof_355=zeros(size(alt));
aaer_prof_355=zeros(size(alt));
d_prof_355=zeros(size(alt));

for i=1:1:length(alt)
    k=alt(i);
    if (k>=height(1,1) & k<=height(1,2))  
        baer_prof_355(i,1)=b_355_1;
        LR_prof_355(i,1)=lr_355_1;
        aaer_prof_355(i,1)=a_355_1;
        d_prof_355(i,1)=d_355_1;
    elseif (k>=height(2,1) & k<=height(2,2))
        baer_prof_355(i,1)=b_355_2;
        LR_prof_355(i,1)=lr_355_2;
        aaer_prof_355(i,1)=a_355_2;
        d_prof_355(i,1)=d_355_2;
    else 
        baer_prof_355(i,1)=0;
        LR_prof_355(i,1)=0;
        aaer_prof_355(i,1)=0;
        d_prof_355(i,1)=0;
    end
end

% array preallocation for time resuming 
baer_prof_532=zeros(size(alt));
LR_prof_532=zeros(size(alt));
aaer_prof_532=zeros(size(alt));
d_prof_532=zeros(size(alt));

for i=1:1:length(alt)
    k=alt(i);
    if (k>=height(1,1) & k<=height(1,2))  
        baer_prof_532(i,1)=b_532_1;
        LR_prof_532(i,1)=lr_532_1;
        aaer_prof_532(i,1)=a_532_1;
        d_prof_532(i,1)=d_532_1;
    elseif (k>=height(2,1) & k<=height(2,2))
        baer_prof_532(i,1)=b_532_2;
        LR_prof_532(i,1)=lr_532_2;
        aaer_prof_532(i,1)=a_532_2;
        d_prof_532(i,1)=d_532_2;
    else 
        baer_prof_532(i,1)=0;
        LR_prof_532(i,1)=0;
        aaer_prof_532(i,1)=0;
        d_prof_532(i,1)=0;
    end
end

% array preallocation for time resuming 
baer_prof_1064=zeros(size(alt));
LR_prof_1064=zeros(size(alt));
aaer_prof_1064=zeros(size(alt));
d_prof_1064=zeros(size(alt));

for i=1:1:length(alt)
    k=alt(i);
    if (k>=height(1,1) & k<=height(1,2))  
        baer_prof_1064(i,1)=b_1064_1;
        LR_prof_1064(i,1)=lr_1064_1;
        aaer_prof_1064(i,1)=a_1064_1;
        d_prof_1064(i,1)=d_1064_1;
    elseif (k>=height(2,1) & k<=height(2,2))
        baer_prof_1064(i,1)=b_1064_2;
        LR_prof_1064(i,1)=lr_1064_2;
        aaer_prof_1064(i,1)=a_1064_2;
        d_prof_1064(i,1)=d_1064_2;
    else 
        baer_prof_1064(i,1)=0;
        LR_prof_1064(i,1)=0;
        aaer_prof_1064(i,1)=0;
        d_prof_1064(i,1)=0;
    end
end

aaer=([aaer_prof_355,aaer_prof_532,aaer_prof_1064]).*10^-3; baer=([baer_prof_355,baer_prof_532,baer_prof_1064]).*10^-3; 
lraer=[LR_prof_355,LR_prof_532,LR_prof_1064]; daer=[d_prof_355,d_prof_532,d_prof_1064];

y1=[0:sm1:1]; y_noise_1=abs(y1.*(1+0.55.*((-1).^randsrc(1,length(y1),[0,1])).*rand(size(y1)))); z1=abs(0.5-abs(0.5-y_noise_1));znorm1=z1./trapz(z1);
y2=[0:sm2:1]; y_noise_2=abs(y2.*(1+0.55.*((-1).^randsrc(1,length(y2),[0,1])).*rand(size(y2)))); z2=abs(0.5-abs(0.5-y_noise_2));znorm2=z2./trapz(z2);

% Break the two steping layers to two different profiles 
kk=find(alt<=height(1,2));
kkk=find(alt<=height(2,1));
aaer_1=zeros(size(aaer)); baer_1=zeros(size(baer));  
aaer_2=zeros(size(aaer)); baer_2=zeros(size(baer));  

aaer_1=[aaer(1:kk(end),:); zeros(length(aaer)-kk(end),3)]; baer_1=[baer(1:kk(end),:); zeros(length(baer)-kk(end),3)];
aaer_2=[zeros(size(aaer(1:kkk(end)-1,:))); aaer(kkk(end):end, :)]; baer_2=[zeros(size(baer(1:kkk(end)-1,:))); baer(kkk(end):end, :)];

w_aaer_1=zeros(size(aaer_1));
w_baer_1=zeros(size(baer_1));
w_lraer_1=zeros(size(lraer));
w_daer_1=zeros(size(daer));

w_aaer_2=zeros(size(aaer_2));
w_baer_2=zeros(size(baer_2));
w_lraer_2=zeros(size(lraer));
w_daer_2=zeros(size(daer));

y1=[0:sm1:1]; y_noise_1=abs(y1.*(1+0.55.*((-1).^randsrc(1,length(y1),[0,1])).*rand(size(y1)))); z1=abs(0.5-abs(0.5-y1));znorm1=z1./trapz(z1);
y2=[0:sm2:1]; y_noise_2=abs(y2.*(1+0.55.*((-1).^randsrc(1,length(y2),[0,1])).*rand(size(y2)))); z2=abs(0.5-abs(0.5-y2));znorm2=z2./trapz(z2);


% y1=[0:sm1:1]; y_noise_1=abs(y1.*(1+0.55.*((-1).^randsrc(1,length(y1),[0,1])).*rand(size(y1)))); z1=abs(0.5-abs(0.5-y_noise_1));znorm1=z1./trapz(z1);
% y2=[0:sm2:1]; y_noise_2=abs(y2.*(1+0.55.*((-1).^randsrc(1,length(y2),[0,1])).*rand(size(y2)))); z2=abs(0.5-abs(0.5-y_noise_2));znorm2=z2./trapz(z2);

for i=1:1:3;
    w_aaer_1(:,i)=conv(aaer_1(:,i), znorm1, 'same');
    w_baer_1(:,i)=conv(baer_1(:,i), znorm1, 'same');
end

for i=1:1:3;
    w_aaer_2(:,i)=conv(aaer_2(:,i), znorm2, 'same');
    w_baer_2(:,i)=conv(baer_2(:,i), znorm2, 'same');
end

w_aaer=w_aaer_1+w_aaer_2;
w_baer=w_baer_1+w_baer_2;

% lraer=w_aaer./w_baer;

% Calculate_AOD
w_AOD=zeros(size(aaer));
AOD=zeros(size(aaer));
for i=2:1:length(aaer);
    for j=1:1:3;
        w_AOD(1,j)=0;
        w_AOD(i,j)=w_AOD(i-1,j)+(0.5.*(alt(i)-alt(i-1)).*(w_aaer(i-1,j)+w_aaer(i,j)));
        AOD(1,j)=0;
        AOD(i,j)=AOD(i-1,j)+(0.5.*(alt(i)-alt(i-1)).*(aaer(i-1,j)+aaer(i,j)));
    end
end

if export==0;
    export_aaer=aaer;
    export_baer=baer;
    export_lraer=lraer;
    export_daer=daer;
    export_AOD=AOD;
elseif export==1;
    export_aaer=w_aaer;
    export_baer=w_baer;
    export_lraer=lraer;
    export_daer=daer;
    export_AOD=AOD;
else
    
end

baer_cross=(daer./100).*(export_baer./(1+daer./100));
baer_parallel=(export_baer./(1+daer./100));

% figure; subplot(1,3,1); plot(aaer(:,1), alt, aaer(:,2), alt, aaer(:,3), alt, w_aaer(:,1), alt, w_aaer(:,2), alt, w_aaer(:,3), alt); myfig; subplot(1,3,2); plot(baer(:,1), alt, baer(:,2), alt, baer(:,3), alt, w_baer(:,1), alt, w_baer(:,2), alt, w_baer(:,3), alt); myfig;...
%     subplot(1,3,3);plot(lraer(:,1), alt, lraer(:,2), alt, lraer(:,3), alt); myfig;

% figure; subplot(1,4,1); plot(aaer(:,1), alt, aaer(:,2), alt, aaer(:,3), alt, w_aaer(:,1), alt, w_aaer(:,2), alt, w_aaer(:,3), alt); myfig; subplot(1,4,2); plot(baer(:,1), alt, baer(:,2), alt, baer(:,3), alt, w_baer(:,1), alt, w_baer(:,2), alt, w_baer(:,3), alt); myfig;...
%     subplot(1,4,3);plot(lraer(:,1), alt, lraer(:,2), alt, lraer(:,3), alt); myfig;...
%     subplot(1,4,4);plot(AOD(:,1), alt, AOD(:,2), alt, AOD(:,3), alt); myfig;

% figure; subplot(1,4,1); plot(aaer(:,1), alt, aaer(:,2), alt, aaer(:,3), alt);  xlabel('\alpha aer [1/m]'); ylabel('Height [m, a.s.l.]'); myfig; subplot(1,4,2); plot(baer(:,1), alt, baer(:,2), alt, baer(:,3), alt); xlabel('\beta aer [1/m/sr]'); myfig;...
%     subplot(1,4,3);plot(lraer(:,1), alt, lraer(:,2), alt, lraer(:,3), alt); xlabel('Saer [sr]'); myfig;...
%     subplot(1,4,4);plot(daer(:,1), alt, daer(:,2), alt, daer(:,3), alt); xlabel('\delta [%]'); myfig;

% figure; subplot(1,4,1); plot(export_aaer(:,1), alt, export_aaer(:,2), alt, export_aaer(:,3), alt); xlabel('Aerosol Ext. Coeff. [1/m]'); ylabel('Height [m, a.s.l.]'); myfig;...  
%         subplot(1,4,2); plot(export_baer(:,1), alt, export_baer(:,2), alt, export_baer(:,3), alt); xlabel('Aerosol Back. Coeff. [1/m/sr]'); myfig;...
%         subplot(1,4,3); plot(export_daer(:,1), alt, export_daer(:,2), alt, export_daer(:,3), alt); xlabel('Particle Depol. [%]'); myfig;...
%         subplot(1,4,4);plot(export_AOD(:,1), alt, export_AOD(:,2), alt, export_AOD(:,3), alt); xlabel('AOD'); myfig;

end
