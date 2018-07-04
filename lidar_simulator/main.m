clc;
clear all;

%% import the inputs from the inputs.xlsx spreadsheet

[mydata, myheader] = xlsread('near_range.xlsx');
for i = 4:length(myheader)
      % compose a command to assign each column to a variable with the same 
      % name as the header
      if isnan(mydata(i-3,1))==1;
         commandExec = [myheader{i,3}, ' = ', 'mydata(' num2str(i-3), ',[2:end]);'];
         % execute the composed command to actually create the variable
         evalin('base', commandExec);
      else
           commandExec_const=[myheader{i,3}, ' = ', 'mydata(' num2str(i-3), ',1);'];
           evalin('base', commandExec_const);

      end
end

%% Some additional inputs for the atmosphere and the solar background conditions 
alt=[0:7.5:10000]';

Pzero=101325; % Pa
Tzero=288.15; % K

humidity=50; % the relative humidity in the atmosphere (0-100%) 0
xCO2=360; % the concentration of CO2 in the atmosphere in ppmv  385

AOD_1=0.1; % the AOD at 532 nm at the bottom layer  
AOD_2=0.2; % the AOD at 532 nm at the lofted layer
height=[0, 1200; 3000, 4000]; % the height and the depth of the layers 
sm1=0.008; % smooth filter for the bottom layer 
sm2=0.008; % smooth filter for the lofted layer
export=1; % keeping the step (0) or filtered (1) and export the corresponding profiles

aerosol_type={'Industrial'; 'Dust'};

diffused_irrad=616; % worst case scenario with intense forwared scattering from water clouds at low altitude (4km) with COD=2
                    % the input parameter should be in mW/m2/nm. 616.1 mW/m2/nm @ 355 nm - 1177.8 mW/m2/nm @ 532 nm and 402.9 mW/m2/nm @ 1064 nm

                     
snr=50; % constant in the atmospheric profile, white gaussian noise

ref_height=5800; %the reference height in m for the lidar signal inversion
lr=45; %the lidar ratio in sr for the lidar signal inversion 57.785 is themean value of a industrial layer (70.91 sr) and a dust (51.8 sr) at 532 nm
% For a reason that i cannot understand i have to user lower LR value than
% the mean of the profile otherwise i take negative values in between of
% the layers (where molecular exist)

br=1.0; % the backscattering ratio for the lidar signal inversion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% call the appropriate routines for the lidar signal construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Forwared problem-From the atmsophere and systems characteristics estimate the final recorded signal

%% call for overlap function
[overlap, DFO, RFOV, Atilt_max]=calculate_overlap(alt, lamda, bexp, dred, FT, DT, Dfieldstop, DTL, Atilt, TFOV, DL, []);

%% US standard atmosphere function
[temp,press,rho]=stdatmo(alt,Tzero);

%% molecular atmosphere
[amol, bmol_c, bmol_t, dmol_c, dmol_t, k_c, k_t, Cs]=ray_cabannes(lamda, temp, press, humidity, xCO2);

%% pertubation-aerosol layers
[aaer, baer, lraer, daer, baer_parallel, baer_cross, AOD]=make_atmosphere(height, alt, AOD_1, AOD_2, sm1, sm2, export, aerosol_type);

%% Calculate the bg signal introduced by downward diffused irradiance
[N_solar, N_solar_photons]=calculate_solar_noise_simple(qe_detector, diffused_irrad, t_sampling, IFF_Tmax, IFF_eq_bdw, RFOV, DT, trans_R, lamda);

%% Construct signal
[P, P_photons]=signal_construction(baer(:,2), aaer(:,2), bmol_t, amol, alt, overlap, p_length, Eo, DT, trans_R, trans_T, qe_detector, lamda);
%P_final=P_photons;%+N_solar_photons; % add the background 
% If i add the solar photons the solution becomes totally unstable for a
% reason that i cannot understant

P_photons(isnan(P_photons))=0;
P_final=awgn(P_photons,snr,'measured');%+N_solar_photons; 

%% Inverse problem - From the Singal go back to the aerosol profiles

%% calculate the rcs signal
[RCS]=calculate_RCS(alt, P_final);
[bp, ap, erbpk, erapk]=klett_old(alt, [], bmol_t, amol, ref_height, lr, br, 100, RCS, []);

