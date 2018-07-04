function [amol, bmol_c, bmol_t, dmol_c, dmol_t, k_c, k_t, Cs]=ray_cabannes(lamda, temp, press, humidity, xCO2)

% INPUTS: 
% lamda: is the wavelength in nm 
% temp: the temperature profile in K
% press: the pressure profile in Pa

% OUTPUTS: 
% amol: molecular extinction coefficient in 1/m 
% bmol_c: backscatter molecular coefficient (1/m/sr) for the Cabannes line 
% bmol_t: backscatter molecular coefficient (1/m/sr) for the total Rayleigh spectrum
% dmol_c: the molecular depolarization ratio for the Cabannes line 
% dmol_t: the molecular depolarization ratio for the total Rayleigh spectrum 

% Calculate the molecular number density in real temperature pressure conditions
% k=1.3806504e-23; % Boltzman constant J/K taken from VFs document Rayleigh scattering coefficients and linear depolarization ratios at several EARLINET lidar wavelegths v 1.4 pp4

k=1.3806488e-23; % Boltzman constant J/K taken from http://physics.nist.gov/cuu/Constants/index.html
N=press./(k.*temp); % in #/m3 from eq 10 of VFs document

% Calculate the refractive index of air in standard atmospheric conditions 
[Fk, nr, epsilon, mw, lamda_air]=king_factor_air_ref(lamda, 288.15, 101325, humidity, xCO2);

% Na=6.02214129e+23; % Avogadros number 1/mol taken from http://physics.nist.gov/cuu/Constants/index.html
% % The density of air for standard atmospheric conditions is (std_rho) and can be calculated from stdatmo routine but for speeding up i just give it as input 
% std_rho=1.225000018124288e+3; % in g/m3
% % Calculate the temperature-pressure and molecular number density in standard atmospheric conditions molecules/m^3
% Ns=((std_rho)/mw)*Na; %the molecular number density 1/m3 calculated for standard atmospheric conditions and variable mw
% 
% % The above two formulations have been done in order to avoide using a
% % constant Ns number from bibliography like 
Ns=2.5469514145749058e+25; %the molecular number density 1/m3 - Bucholtz first paragraph I think VFs calculations are taken the Ns value as constant something that may be wrong since its strongly depended on the airs molecular weight. Further more the dry air is more dense than the moist !!!

% Calculate the molecular backscatter for the cabannes line
amol=N.*((24*pi^3)./((lamda.*10^-9).^4.*Ns.^2)).*((nr.^2-1).^2./(nr.^2+2).^2).*Fk;

% Calculate the molecular depolarization ratio for cabannes line Eq.16 VF  
dmol_c=(3.*Fk-3)./(4.*Fk+36);

% Calculate the molecular depolarization ratio for total Rayleigh spectrum Eq.15 VF  
dmol_t=(3.*Fk-3)./(4.*Fk+6);

% Calculate the conversion factors for total Rayeleigh and cabannes lines Eq.20 && 22 VF  
k_t=(1+2.*dmol_t)./(1+dmol_t); % Eq 20
k_c=(1+2.*dmol_t)./(1-(3/4).*dmol_t); % Eq 22

% Calculate the molecular backscatter for total Rayleigh and cabannes lines Eq 18 && 20 && 21
bmol_t=amol./(k_t.*8*pi./3);
bmol_c=amol./(k_c.*8*pi./3);

Cs=(1./k).*((24*pi^3)./((lamda.*10^-9).^4.*Ns.^2)).*((nr.^2-1).^2./(nr.^2+2).^2).*Fk;

end
