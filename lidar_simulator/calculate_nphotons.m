function photons_emitted=calculate_nphotons(lamda, Eo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%caclulates the number of photons emited from the laser beam 
%https://www.wyzant.com/resources/answers/1065/how_many_photons_are_produce
%d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS 
% lamda (nm); Eo is the laser pulse energy (mJ)

c=299792458; % speed of light (in m/sec)

h=6.63*10^-34; % Planck's constant (in J*s)

lamda=lamda.*10^-9; %conversion of wavelength from nm to m

E=h.*c./lamda; % The Energy of a single photon at wavelength lamda (in J)

Eo=Eo.*10^-3; % The initial energy of the entire laser pulse is converted from mJ to J  

photons_emitted=Eo./E; % Division of the energy of the whole pulse by the energy per photon (at a wavelength) to get the number of photons

end