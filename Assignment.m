clear all;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. These are used by all Nasa Functions. 
global Runiv Pref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions. 
%  For the final assignment take the ones from the specific case you are supposed to do.                  
v1=200;Tamb=250;P3overP2=7;Pamb=55*kPa;mfurate=0.68*kg/s;AF=85.65;             % These are the ones from the book
cFuel='CH4';                                                           % Pick Gasoline as the fuel (other choices check Sp.Name)
%% Select species for the case at hand
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79];                                                   % Order is important. Note that these are molefractions
MAir = Xair*Mi';                                                            % Row times Column = inner product 
Yair = Xair.*Mi/MAir;                                                       % Vector. times vector is Matlab's way of making an elementwise multiplication
%% Fuel composition
Yfuel = [1 0 0 0 0];                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:5000];NTR=length(TR);
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR 
    hia(:,i) = HNasa(TR,SpS(i));                                            % hia is a NTR by 5 matrix
    sia(:,i) = SNasa(TR,SpS(i));                                            % sia is a NTR by 5 matrix
    uia(:,i) =  UNasa(TR,SpS(i));
end
hair_a= Yair*hia';                                                          % Matlab 'inner product': 1x5 times 5xNTR matrix muliplication, 1xNTR resulT -> enthalpy of air for range of T 
sair_a= Yair*sia';                                                          % same but this thermal part of entropy of air for range of T
uair_a=Yair*uia';
% whos hia sia hair_a sair_a                                                  % Shows dimensions of arrays on commandline
%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the interpolation method
% Bisection is in the next 'cell'
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using INTERPOLATION
cMethod = 'Interpolation Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
T2 = interp1(hair_a,TR,h2);                                                 % Interpolate h2 on h2air_a to approximate T2. Pretty accurate
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';                                                        % Single value (1x5 times 5x1). Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral part of th eentropy)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total specific entropy
S2  = s2thermal - Rg*log(P2/Pref);
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
T2int = T2;

%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the Bisection method
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using bisection (https://en.wikipedia.org/wiki/Bisection_method)
cMethod = 'Bisection Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
TL = T1;
TH = 1000;                                                                  % A guess for the TH (must be too high)
iter = 0;
while abs(TH-TL) > 0.01
    iter = iter+1;
    Ti = (TL+TH)/2;
    for i=1:NSp
        hi2(i)    = HNasa(Ti,SpS(i));
    end
    h2i = Yair*hi2';                                                        % Single value (1x5 times 5x1). Intermediate value
    if h2i > h2
        TH = Ti; % new right boundary
    else
        TL = Ti; % new left boundary
    end
end
T2 = (TH+TL)/2;
T2bis = T2;
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total entropy stage 1
S2  = s2thermal - Rg*log(P2/Pref);                                          % Total entropy stage 2
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
%% Difference between two approaches: so close but not identical
fprintf('----------------------------------------------\n%8s| %9.4f %9.4f  [K]\n----------------------------------------------\n','T2-int vs T2-bis',T2int,T2bis);
%% Here starts your part (compressor,combustor,turbine and nozzle). ...
% Make a choice for which type of solution method you want to use.





%Compressor part



function [P3, T3, h3, Wcomp] = compressor(P2, T2, s2, P3overP2, Yair, SpS, mfurate, AF)

% Inputs:
%   P2, T2    - inlet pressure (Pa) and temperature (K)
%   s2        - inlet entropy (J/kg-K)
%   P3overP2  - pressure ratio
%   Yair      - species mass fractions (1 x NSp)
%   SpS       - species identifiers (used by HNasa, SNasa)
%   mfurate   - fuel mass flow (kg/s)
%   AF        - air-fuel ratio

% Outputs:
%   P3, T3    - outlet pressure (Pa), outlet temperature (K)
%   h3        - outlet enthalpy (J/kg)
%   Wcomp     - compressor work (W)

% Pressure ratio
P3 = P2 * P3overP2;

% Solve T3 from isentropic relation (s3 = s2)
entropy_residual = @(T) SNasa(T, SpS, Yair, P3) - s2;
Tguess = T2 * (P3overP2)^0.28; % crude gamma-based guess
T3 = fzero(entropy_residual, Tguess);

% Outlet enthalpy from NASA polynomials
NSp = length(SpS);
hi3 = zeros(1, NSp);
for i = 1:NSp
    hi3(i) = HNasa(T3, SpS(i));
end
h3 = Yair * hi3'; % mixture enthalpy

% Inlet enthalpy (needed for Δh)
hi2 = zeros(1, NSp);
for i = 1:NSp
    hi2(i) = HNasa(T2, SpS(i));
end
h2 = Yair * hi2';

% Compressor work
mflow2 = AF * mfurate;        % air mass flow rate
Wcomp  = mflow2 * (h3 - h2);

end













%combustor part


%Combustion process
%  CH4+2O2=CO2+2H2O
a=1;b=2;c=1;d=2;
%Calculate Equivulence ratio
AFstoich=((b/a)/Xair(2))*MAir/Mi(1); %Stoichiometric A/F ratio
equ=AFstoich/AF; %Equivulence ratio

%Calculate u before ignition 
Yfuel=1/(1+AF)
YO2=Yair(2)*(AF/(1+AF))
YN2=Yair(5)*(AF/(1+AF))

Xfuel=(Yfuel/Mi(1))/(YO2/Mi(2)+Yfuel/Mi(1)+YN2/Mi(5))
XO2=(YO2/Mi(2))/(YO2/Mi(2)+Yfuel/Mi(1)+YN2/Mi(5))
XN2=(YN2/Mi(5))/(YO2/Mi(2)+Yfuel/Mi(1)+YN2/Mi(5))
U_beforecomb=Yfuel*UNasa(T2,SpS(1))+YN2*UNasa(T2,SpS(5))+YO2*UNasa(T2,SpS(2));







%Calculate u after ignition 

%set n_air=1
n_fuel=Xfuel
% moles changed (per reaction stoichiometry)
dXO2  = (b/a) * Xfuel; % O2 consumed
dXCO2 = (c/a) * Xfuel; % CO2 produced
dXH2O = (d/a) * Xfuel; % H2O produced




X_aftercomb = [0 XO2-dXO2 dXCO2 dXH2O XN2]
M_aftercomb = X_aftercomb*Mi';                                                            % Row times Column = inner product 
Y_aftercomb = X_aftercomb.*Mi/M_aftercomb; 
uair_aftercomb=Y;
%U_beforecomb=U_aftercomb find T after comustion
U_aftercomb=U_beforecomb;
T2 = interp1(uair_aftercomb,TR,U_aftercomb); 
%Use ideal gas 





%turbine part






%nozzle part