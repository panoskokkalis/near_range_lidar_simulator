function [Fk, nr, epsilon, mw, lamda_air]=king_factor_air_ref(lamda, temp, press, humidity, xCO2)

% INPUTS:
% lamda: the laser wavelength in vacuum in nm
% temp: the temperature profile in Kelvin
% pressure: the pressure profile in Pa
% humidity: relative humidity %
% xCO2: the volume concentration of CO2 in the atmosphere in ppmvthe CO2 concentration in the atmosphere a value of 385 to 450 ppmv
% concidered to be realistic 

% OUTPUTS:
% Fk the King Factor regarding anisotropy of the molecule, depended only on
% lamda
% nr the real part of refractive index of air 
% mw is the air molecular wheigth (in g/mol) depends on relative humidity and xCO2
% lamda_air the laser wavelength in air (in nm) since it is storngly depended on
% the refractive index and number desnity of air 

% Caluclation of the real part of refractive index of air according to
% Ciddor et al 2002 and Tomasi et al 2005 follwing the set of equations 10-19 


x1=temp-273.15; x2=press; x3=humidity;
S=1e6/lamda^2;
T=temp;

% IAPWS Check the humidity conversion formulas from vaisala pdf !!! 
% W. Wagner and A. Pruß:" The IAPWS Formulation 1995 for the
% Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use ",
% Journal of Physical and Chemical Reference Data, June 2002 ,Volume 31, Issue 2, pp.
% 387535
% http://emtoolbox.nist.gov/Wavelength/Documentation.asp#EdlenorCiddor
K1=1.16705214528e3;
K2=-7.24213167032e5;
K3=-1.70738469401e1;
K4=1.20208247025e4;
K5=-3.23255503223e6;
K6=1.49151086135e1;
K7=-4.82326573616e3;
K8=4.05113405421e5;
K9=-2.38555575678e-1;
K10=6.50175348448e2;
Omega=T+K9./(T-K10);
A=Omega.^2+K1*Omega+K2;
B=K3*Omega.^2+K4*Omega+K5;
C=K6*Omega.^2+K7*Omega+K8;
X=-B+sqrt(B.^2-4*A.*C);
psv1=1e6*(2*C./X).^4;

% Over Ice
A1=-13.928169;
A2=34.7078238;
Theta=T/273.16;
Y=A1*(1-Theta.^-1.5)+A2*(1-Theta.^-1.25);
psv2=611.657*exp(Y);

psv=psv1.*(x1>=0)+psv2.*(x1<0);

% Convert humidity to mole fraction for Ciddor
alpha=1.00062;
beta=3.14e-8;
gamma=5.60e-7;
fpt=alpha+beta*x2+gamma*x1.^2;
xv=(x3/100).*fpt.*psv./x2;

pres_vapour=(x3/100).*psv; % the partial water vapour pressure 

% Ciddor Equation
w0=295.235;w1=2.6422;w2=-0.03238;w3=0.004028;
k0=238.0185;k1=5792105;k2=57.362;k3=167917;
a0=1.58123e-6;a1=-2.9331e-8;a2=1.1043e-10;
b0=5.707e-6;b1=-2.051e-8;
c0=1.9898e-4;c1=-2.376e-6;
d=1.83e-11;e=-0.765e-8;
pR1=101325;TR1=288.15;
Za=0.9995922115;
rhovs=0.00985938;
R=8.314472;Mv=0.018015;
ras=1e-8*(k1/(k0-S)+k3/(k2-S));
rvs=1.022e-8*(w0+w1*S+w2*S^2+w3*S^3);
Ma=0.0289635+1.2011e-8*(xCO2-400);
raxs=ras*(1+5.34e-7*(xCO2-450));
Zm=1-(x2./T).*(a0+a1*x1+a2*x1.^2+(b0+b1*x1).*xv+ ... 
    (c0+c1*x1).*xv.^2)+(x2./T).^2.*(d+e*xv.^2);
rhoaxs=pR1*Ma/(Za*R*TR1);
rhov=xv.*x2.*Mv./(Zm*R.*T);
rhoa=(1-xv).*x2.*Ma./(Zm*R.*T);

nr=1+(rhoa./rhoaxs)*raxs+(rhov/rhovs)*rvs;

% Kings refractive index according to Eq 20 found in Tomasi et al 2005
lamda_micro_meter=lamda./1000;
F1=1.034+3.17*10^-4*lamda_micro_meter^-2; C1=0.78084; % dry air and N2 molecules contribution
F2=1.096+1.385*10^-3*lamda_micro_meter^-2+1.448*10^-4*lamda_micro_meter^-4; C2=0.20946; % dry air and O2 molecules contribution
F3=1; C3=0.00934; % dry air and Ar molecules contribution
F4=1.15; C4=10^-6*xCO2;% dry air and CO2 molecules contribution-C is the CO2 concentration in ppmv
F5=1.001; % various moisture conditions 

Fk=(C1*F1+C2*F2+C3*F3+C4*F4+(pres_vapour./x2).*F5)./(C1+C2+C3+C4+(pres_vapour./x2));

% Calculate epsilon parameter involved in Fk for total molecular extinction
% coeff in equation 14 of VF pp 5
epsilon=(45*Fk-45)/10;

% Calulate the molecular weight of air as the sum of molecular wheights of the individual components
% http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
% http://www.engineeringtoolbox.com/molecular-weight-gas-vapor-d_1156.html
% http://renewableenergy.wikia.com/wiki/Molecular_weight_of_dry_air
mw=0.78084*28.02+0.20946*32+0.00934*39.94+xCO2*10^-6*44.01; % in g/mol the H2O is not included 

% caclulation of the lamda air according to
% http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/pyasl_wvlconv.html#PyAstronomy.pyasl.airtovac
% and Atomic data for resonance absorption lines. I - Wavelengths longward
% of the Lyman limit [ Erratum: 1992ApJS...81..883M] Morton, Donald C
lamda_air=lamda./nr;

end
