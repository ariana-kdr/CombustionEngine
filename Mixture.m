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
cFuel='Gasoline';  % Initial fuel
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Select species for the case at hand
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2','C2H5OH'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79 0];                                             % Order is important (check myFind). Note that these are molefractions
MAir = Xair*Mi';                                                            % Row times Column = inner product 
Yair = Xair.*Mi/MAir;                                                       % Vector. times vector is Matlab's way of making an elementwise multiplication

%% Percentage ethanol
prompt = 'Enter percentage of ethanol (0.0-1.0): ';
dlgtitle = 'Requested input';
dlg_default = {'0'};
dlg = inputdlg(prompt, dlgtitle, 1 , dlg_default);     %Fraction ranging from 0.0 to 1.0
e = str2double(dlg(1));
%% Mass and volume comparison
V_e = e;                                           %volume percentage ethanol
V_g = 1-e;                                         %volume percentage gasoline
rho_g = 748.9;                                        %density gasoline [kg/m^3]
rho_e = 789;                                          %density ethanol [kg/m^3]
Yfuel = ([(V_g*rho_g) (V_e*rho_e) 0 0 0 0])/...   %mass fraction of total fuel
    (sum([(V_g*rho_g) (V_e*rho_e) 0 0 0 0]));
%% Number of moles and fractions
Xfuel = (Yfuel./Mi)/sum(Yfuel./Mi); 
C = SpS(1).Elcomp(3);                                   %Carbon atoms in gasoline
H = SpS(1).Elcomp(2);                                   %Hydrogen atoms in gasoline
alpha = e/(1-e);
%Reaction equation to compute masses
a = Xfuel(1)+alpha*Xfuel(2);
c = a*(C + alpha*2);
d = a*(H+alpha*6)/2;
b = c + 0.5*d - 0.5*a*alpha;
e = 3.76*b;

m_f = a*(Mi(1) + alpha*Mi(6));
m_a = b*Mi(2) + e*Mi(5);

AF = m_a / m_f;