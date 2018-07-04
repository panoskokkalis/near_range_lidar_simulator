function [N_solar, N_solar_photons]=calculate_solar_noise_simple(qe_detector, diffused_irrad, sampling_time_bin, Tmax, IF_eq_bdw, RFOV, DT, trans_R, lamda)

sampling_time_bin=sampling_time_bin.*10^-9; % sec
DT=DT.*10^-3; %m
trans_R=trans_R./100; %unitless
qe_detector=qe_detector./100; %unitless
RFOV=RFOV; % Is already calculated in rad !!!
Tmax=Tmax./100;

IF_eq_bdw=IF_eq_bdw*10^-9; %m

diffused_irrad=diffused_irrad.*10^9; %mW/m2/m

N_solar=qe_detector.*(diffused_irrad./pi).*(pi.*(DT).^2./4).*((pi*(RFOV^2))).*sampling_time_bin.*trans_R.*IF_eq_bdw.*Tmax;

N_solar_photons=calculate_nphotons(lamda, N_solar);

end
