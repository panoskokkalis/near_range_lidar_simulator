function [bp, ap, erbpk, erapk]=klett_old(alt, files_summed, bm, am, ref_height, lr, br, nbeans, rcs, erP)

% This function solves the Klett/Fernald backward integration for retrieving backscatter coefficient 
% Copyright (C) 2015,  Kokkalis Panos, panko@noa.gr
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%======================================================================================


%alt altitude m
%files_summed signals that was corrected and then summed
%bm backscatter molecular
%ref_height reference height in m
%lr lidar ratio
%br backscatter ratio
%in the column (channel) selected for retrievals
%nbeans the number of beans above and below ref_height for the molecular
%backscater in the free particles atmoshere


index=find(alt>=ref_height); %&&&&&&&&&&&&&& signal a matrix
bm_z0(:,1)=mean(bm((index(1)-nbeans):(index(1)+nbeans),1));
am_z0(:,1)=mean(am((index(1)-nbeans):(index(1)+nbeans),1));

% rcs(:,1)=(alt(:,1).*alt(:,1)).*files_summed(:,1);
rcs_z0(:,1)=mean(rcs((index(1)-nbeans):(index(1)+nbeans),1));

%  finish=index(1)+nbeans;
 finish=index(1);
% test_rcs=rcs(finish(:,in))
int_beta_m=zeros([index(1)], 1);
int_beta_m(finish,1)=bm_z0;
for i=(finish-1):-1:1;
    int_beta_m(i,1)=int_beta_m(i+1,1)+(0.5.*(alt(i,1)-alt(i+1,1)).*(bm(i+1,1)+bm(i,1)));
end

for i=1:1:finish
    A(i,1)=rcs(i,1).*exp(-2.0.*(lr-8.37758).*int_beta_m(i,1)); %%%% Is it rcs maybe ? 
end

int_A=zeros([index(1), 1]);
int_A(finish,1)=A(finish,1);
% Q2temp(finish,in)=0.0;
for i=(finish-1):-1:1;
    int_A(i,1)=int_A(i+1)+(0.5.*(alt(i)-alt(i+1)).*(A(i+1)+A(i)));
end

 B=rcs_z0/(br.*bm_z0);
%  B=rcs(index(1),1)/(br.*bm(index(1),1));
%  B=2.696e-6/(br.*bm(index(1),in));

for i=1:1:finish
    bp(i,1)=-bm(i,1)+(A(i,1)./(B-2.*lr.*int_A(i,1)));
end
 
ap(:,1)=bp(:,1).*lr;

%Start the error caclulations according to Voler phd due to statistical
%error on signal

if isempty(erP)==0;
    
s1=lr;
s2=8.37758;
s=rcs;

step=ones(size(bp)-1);
A=ones(size(bp)-1);
Z=ones(size(bp)-1);
N=ones(size(bp)-1);
db1ar=ones(size(bp)-1);
db1=ones(size(bp)-1);
db2ar=ones(size(bp)-1);
db2=ones(size(bp)-1);
erbpk=ones(size(bp)-1);
erapk=ones(size(bp)-1);

for i=2:1:finish
    step(i,1)=alt(i)-alt(i-1);
    A(i,1)=(s1-s2).*(bm(i)+bm(i-1)).*step(i);
    Z(i,1)=(s(i)).*exp(A(i));
    N(i,1)=(s(i)./(bp(i)+bm(i))) + s1.*(s(i)+s(i-1).*exp(A(i))).*step(i,1);
    db1ar(i,1)=-(((alt(i).^2)/(bp(i)+bm(i)))+s1.*(alt(i).^2).*step(i)).*Z(i);
    db1(i,1)=db1ar(i,1).*erP(i,1)./(N(i,1).^2);
    db2ar(i,1)=(alt(i-1).^2).*exp(A(i)).*N(i)-s1.*(alt(i-1).^2).*step(i).*exp(A(i))*Z(i);
    db2(i,1)=db2ar(i,1).*erP(i-1)./(N(i,1).^2);

    erbpk(i,1)=10.*sqrt((db1(i,1).^2)+(db2(i,1).^2));
    erapk(i,1)=10.*erbpk(i,1).*s1;
end 

else
    erbpk=0;
    erapk=0;
end
end


