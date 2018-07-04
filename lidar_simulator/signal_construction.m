function [P, P_photons]=signal_construction(baer, ext, bmol, amol, alt, K, p_length, Eo, DT, trans_R, trans_T, qe_detector, lamda)


trans_R=trans_R./100;  trans_T=trans_T./100; qe_detector=qe_detector./100;

c=299792458; %m/s
p_length=p_length*10^-9;

photons_emitted=calculate_nphotons(lamda, Eo);

DT=DT*10^-3;

A=pi*(DT/2)^2;

% calculate the integral of exp
atotal=ext+amol;
btotal=baer+bmol;

integral=zeros(size(atotal));

for i=2:1:length(ext);
    integral(1)=0;
    integral(i)=integral(i-1)+(0.5.*(alt(i)-alt(i-1)).*(atotal(i-1)+atotal(i)));
end

P=Eo*c*(p_length/2).*A.*(btotal./alt.^2).*K.*trans_R.*trans_T.*qe_detector.*exp(-2.*integral); % signal in mJ
P_photons=photons_emitted*c*(p_length/2).*A.*(btotal./alt.^2).*K.*trans_R.*trans_T.*qe_detector.*exp(-2.*integral); % number of photons detected per pulse

end

