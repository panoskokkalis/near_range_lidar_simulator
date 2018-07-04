function [bp_near, ap_near]=near_range_klett(alt, files_summed, bmol_t, amol, ref_height, lr, br, nbeans, rcs_near, erP)

% This function solves the Klett/Fernald backward integration for retrieving backscatter coefficient 
% Copyright (C) 2015,  Kokkalis Panos, panko@noa.gr
%
% This programol is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This programol is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this programol.  If not, see <http://www.gnu.org/licenses/>.


%======================================================================================



index=find(alt>=ref_height); %&&&&&&&&&&&&&& signal a matrix
bmol_t_z0(:,1)=mean(bmol_t((index(1)-nbeans):(index(1)+nbeans),1));
amol_z0(:,1)=mean(amol((index(1)-nbeans):(index(1)+nbeans),1));

% rcs(:,1)=(alt(:,1).*alt(:,1)).*files_summed(:,1);
rcs_z0(:,1)=mean(rcs_near((index(1)-nbeans):(index(1)+nbeans),1));

%  finish=index(1)+nbeans;
 finish=index(1);
% test_rcs=rcs(finish(:,in))
int_beta_m=zeros([index(1)], 1);
int_beta_m(finish,1)=bmol_t_z0;
for i=(finish-1):-1:1;
    int_beta_m(i,1)=int_beta_m(i+1,1)+(0.5.*(alt(i,1)-alt(i+1,1)).*(bmol_t(i+1,1)+bmol_t(i,1)));
end

for i=1:1:finish
    A(i,1)=rcs_near(i,1).*exp(-2.0.*(lr-8.37758).*int_beta_m(i,1));  
end

int_A=zeros([index(1), 1]);
int_A(finish,1)=A(finish,1);
% Q2temp(finish,in)=0.0;
for i=(finish-1):-1:1;
    int_A(i,1)=int_A(i+1)+(0.5.*(alt(i)-alt(i+1)).*(A(i+1)+A(i)));
end

 B=rcs_z0/(br.*bmol_t_z0);
%  B=rcs(index(1),1)/(br.*bmol_t(index(1),1));
%  B=2.696e-6/(br.*bmol_t(index(1),in));

for i=1:1:finish
    bp_near(i,1)=-bmol_t(i,1)+(A(i,1)./(B-2.*lr.*int_A(i,1)));
end
 
ap_near(:,1)=bp_near(:,1).*lr;

%Start the error caclulations according to Voler phd due to statistical
%error on signal

if isempty(erP)==0;
    
s1=lr;
s2=8.37758;
s=rcs_near;

step=ones(size(bp_near)-1);
A=ones(size(bp_near)-1);
Z=ones(size(bp_near)-1);
N=ones(size(bp_near)-1);
db1ar=ones(size(bp_near)-1);
db1=ones(size(bp_near)-1);
db2ar=ones(size(bp_near)-1);
db2=ones(size(bp_near)-1);
erbpk=ones(size(bp_near)-1);
erapk=ones(size(bp_near)-1);

for i=2:1:finish
    step(i,1)=alt(i)-alt(i-1);
    A(i,1)=(s1-s2).*(bmol_t(i)+bmol_t(i-1)).*step(i);
    Z(i,1)=(s(i)).*exp(A(i));
    N(i,1)=(s(i)./(bp_near(i)+bmol_t(i))) + s1.*(s(i)+s(i-1).*exp(A(i))).*step(i,1);
    db1ar(i,1)=-(((alt(i).^2)/(bp_near(i)+bmol_t(i)))+s1.*(alt(i).^2).*step(i)).*Z(i);
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


