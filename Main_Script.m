%% General

clear all;
close all;
clc;

%Add different path folders
relativepath_to_generalfolder = 'General';
addpath(relativepath_to_generalfolder);
relativepath_to_nasafolder = 'Nasa';
addpath(relativepath_to_nasafolder);

%Load Nasa Database
TdataBase = fullfile('Nasa','NasaThermalDatabase.mat');
load(TdataBase);

% Constants 
global Runiv Pref Tref                                                      %Set global variables (so they cannot be changed)
Runiv = 8.314472;                                                           %Universal Gas Constant[J/K*mol]
Pref = 1.01235e5;                                                           %Reference Pressure (1atm) [Pa]
Tref = 298.15;                                                              %Reference Temperature [K]

%Starting with only Gasoline
Species = Sp(myfind({Sp.Name},{'Gasoline','O2','CO2','H2O','N2'}));         %Subselection of species in database
Mi = [Species.Mass];                                                        %Array of molecular masses for each component
n_Nitr = (25*0.79)/(0.21);                                                  %No. moles nitrogen using air composition (79%)

%Mole Fractions using balanced chemical equation
X_reac = [2/(2+25+n_Nitr) 25/(2+25+n_Nitr) 0 0 n_Nitr/(2+25+n_Nitr)];
X_prod = [0 0 16/(16+18+n_Nitr) 18/(16+18+n_Nitr) n_Nitr/(16+18+n_Nitr)];

%Convert to mass fractions
Y_reac = (X_reac.*Mi)/(X_reac*Mi');
Y_prod = (X_prod.*Mi)/(X_prod*Mi');

%Molar mass
M_reac = X_reac*Mi';
M_prod = X_prod*Mi';

%Gas constants
R_reac = Runiv/M_reac;
R_prod = Runiv/M_prod;

%Dimensions (using dimensions found online)
B = 68e-3;                                                                  %Bore [m]
r = 27e-3;                                                                  %Radius of crank shaft [m]
S = 54e-3;                                                                  %Stroke [m]
L = 84.7e-3;                                                                %Length of rod [m]
r_c = 8.5;                                                                  %Compression Ratio [m]
V_d = pi/4*B^2*S;                                                           %Displacement Volume [m^3]
V_c = V_d/(r_c-1);                                                          %Clearance volume [m^3]
%% Idealized Cycle
% Assumptions:
% * Both intake and exhaust gases are _ideal gases_ and are at _atmospheric pressure_
% * _Adiabatic_ compression
% * _Isochoric_ and _Instantaneous_ combustion
% * _Adiabatic_ expansion
% * _Isochoric_ exhaust 
% * _Isothermal/Isobaric_ exhaust/intake stroke

%Create empty matrices for each state property
p_ideal = [];                                                               %Ideal cycle pressure array [Pa]                                                       
V_ideal = [];                                                               %Ideal cycle volume array [m^3]
T_ideal = [];                                                               %Ideal cycle temperature array [K]
m_ideal = [];                                                               %Ideal cycle mass array [kg]
theta_crank = [];                                                           %Ideal cycle crank angle array [rad]

%Time Values
dt = 0.0001;                                                                %Time step to solve differential equations [s]
dthetadt = 3000/60*2*pi;                                                    %Convert speed of engine to [rad/s]
t_end = 2*(60/3000);                                                        %End time [s] (one rotation is 0.01 sec however, one cycle is two revolutions)
time = [0:dt:t_end];                                                        %Time array [s]

%Create vector array of crank angle and volume as a function of time
i = 0;                                                                      %Index Value
for t = [0:dt:t_end]
    i = i+1;                                                                %Update Index
    theta_crank(i) = t*dthetadt;                                            %Calculate angle for each time instance
    V_ideal(i) = Vcyl(theta_crank(i),B,r,L,V_c);                            %Calculate volume for angle using volume function
end
theta_crank_degrees = theta_crank.*(180/pi);                                %Convert crank angle to degrees for easy life

%Set the angle for the start position of each of the stages
intake_stroke_start_theta = 0;
compression_start_theta = 180;
ignition_start_theta = 360;                                                 %These will have to be changed for real cycle
expansion_start_theta = 360;
valve_exhaust_start_theta = 540;
exhaust_stroke_start_theta = 540;

%Find the index in the crank angle vector array for which it is at this angle
intake_stroke_index = max(find(theta_crank_degrees==intake_stroke_start_theta));
compression_start_index = max(find(theta_crank_degrees<compression_start_theta));
ignition_start_index = max(find(theta_crank_degrees<ignition_start_theta));
expansion_start_index = max(find(theta_crank_degrees<expansion_start_theta));
valve_exhaust_start_index = max(find(theta_crank_degrees<valve_exhaust_start_theta));
exhaust_stroke_start_index = max(find(theta_crank_degrees<exhaust_stroke_start_theta));

%Find the time instant for each start of the strokes
intake_stroke_start_time = time(intake_stroke_index);
compression_start_time = time(compression_start_index);
ignition_start_time = time(ignition_start_index);
expansion_start_time = time(expansion_start_index);
valve_exhaust_start_time = time(valve_exhaust_start_index);
exhaust_stroke_start_time = time(exhaust_stroke_start_index);

%% Ideal Intake Stroke

i = 0;                                                                      %Index Value
for t = [intake_stroke_start_time:dt:compression_start_time]                %Start to end time
    i = i+1;
    p_ideal(i) = Pref;                                                      %Isobaric
    T_ideal(i) = Tref;                                                      %Isothermal
    m_ideal(i) = (p_ideal(i)*V_ideal(i))/(T_ideal(i)*R_reac);               %Mass intake calculated using ideal gas law
end

P1 = p_ideal(i);
V1 = V_ideal(i);
T1 = T_ideal(i);
m1 = m_ideal(i);

%% Ideal Compression

for t = [compression_start_time+dt:dt:ignition_start_time]
    i = i+1;                                                                %Keep updating same index as last time

    %Calculate the volumetric heat capacity for each molecule
    for j = [1:length(Species)]
        cvj(j) = CvNasa(T_ideal(i-1),Species(j));
    end
    
    Cv = Y_reac*cvj';                                                       %Calculate total heat capacity using mass fractions

    m_ideal(i) = m_ideal(i-1);                                              %Mass stays constant since valves are closed

    %Euler Forward Method:
    dTdt = -(R_reac*T_ideal(i-1))/(V_ideal(i-1)*Cv)*dVdt(theta_crank(i),dthetadt,B,r,L);
    T_ideal(i)= T_ideal(i-1)+dt*dTdt;

    p_ideal(i)= (m_ideal(i)*R_reac*T_ideal(i))/(V_ideal(i));                %Ideal Gas Law
end

plot(V_ideal(1:i),p_ideal)

%% Ideal Ignition


%% Ideal Power Stroke 


%% Ideal Valve Exhaust


%% Ideal Exhaust Stroke


%% Ideal Intake Stroke


%% Functions

function [Vcyl] = Vcyl(theta,B,r,L,V_c)
%This function is used to calculate the volume of the cylinder 
%   See lecture for derivation of formula
    x = r*cos(theta)+sqrt(L^2-r^2*sin(theta)^2);                            %Distance from center of crankshaft to center of piston head [m]
    d = L+r-x;                                                              %Distance of top of piston head to its maximum height [m]

    Vcyl = pi/4*B^2*d+V_c;                                                  %Cylinder Volume [m^3]
end

function [dVdt] = dVdt(theta,dthetadt,B,r,L)
%Calculate derivative of volume with respect to time for given theta
    dVdt = pi/4*B^2*(dthetadt*r*sin(theta)+0.5*(dthetadt*r^2*sin(2*theta))/(sqrt(L^2-0.5*r^2*(1-cos(2*theta)))));
end
