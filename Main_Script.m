% General

%Known Issues:
% - Indexing in for loop is not perfect yet
% - For loop expans each state property matrix increasing run time
% - Since volume is a function of time ignition and valve exhaust are not
% instantaneous

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

w = 50;                                                                     %Angular Veloctiy [hz]                                                                      
%% Idealized Cycle
% Assumptions:
% - Both intake and exhaust gases are _ideal gases_ and are at _atmospheric 
%pressure_
% - _Adiabatic_ compression
% - _Isochoric_ and _Instantaneous_ combustion
% - _Adiabatic_ expansion
% - _Isochoric_ exhaust 
% - _Isothermal/Isobaric_ exhaust/intake stroke

% - Angular velocity is _constant_ 

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
theta_intake_stroke = 0;
theta_compression = 180;
theta_ignition = 360;                                                       %These will have to be changed for real cycle
theta_expansion = 360;
theta_valve_exhaust = 540;
theta_exhaust_stroke = 540;

%Find the index in the crank angle vector array for which it is at this
%angle
index_intake_stroke = max(find(theta_crank_degrees==theta_intake_stroke));
%Since angle is never perfectly 180 take max index where it is lower and
%add one (cheap workaround)
index_compression = max(find(theta_crank_degrees<theta_compression))+1;
index_ignition = max(find(theta_crank_degrees<theta_ignition))+1;
index_expansion = max(find(theta_crank_degrees<theta_expansion))+1;
index_valve_exhaust = max(find(theta_crank_degrees<theta_valve_exhaust))+1;
index_exhaust_stroke = max(find(theta_crank_degrees<theta_exhaust_stroke))+1;

%Find the time instant for each start of the strokes
time_intake_stroke = time(index_intake_stroke);
time_compression = time(index_compression);
time_ignition = time(index_ignition);
time_expansion = time(index_expansion);
time_valve_exhaust = time(index_valve_exhaust);
time_exhaust_stroke = time(index_exhaust_stroke);

%% Ideal Intake Stroke

i = 0;                                                                      %Index Value
for t = [time_intake_stroke:dt:time_compression]                            %Start to end time
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

for t = [time_compression+dt:dt:time_ignition]
    i = i+1;                                                                %Keep updating same index as last time
    
    %Calculate the volumetric heat capacity for each molecule
    for j = [1:length(Species)]
        Cvj(j) = CvNasa(T_ideal(i-1),Species(j));
    end
    Cv = Y_reac*Cvj';                                                       %Calculate total heat capacity using mass fractions

    m_ideal(i) = m_ideal(i-1);                                              %Mass stays constant since valves are closed

    %Euler Forward Method:
    dTdt = -(R_reac*T_ideal(i-1))/(V_ideal(i-1)*Cv)*dVdt(theta_crank(i),dthetadt,B,r,L);
    T_ideal(i)= T_ideal(i-1)+dt*dTdt;

    p_ideal(i)= (m_ideal(i)*R_reac*T_ideal(i))/(V_ideal(i));                %Ideal Gas Law
end

P2 = p_ideal(i);
V2 = V_ideal(i);
T2 = T_ideal(i);
m2 = m_ideal(i);
%% Ideal Ignition

%Calculate lower heating value
for j = [1:length(Species)]
    h_0(j) = HNasa(Tref,Species(j));                                        %Calculate enthalpy of formation for each species
    Cv_2i(j) = CvNasa(T2,Species(j));                                       %Volumetric Heat Capacity at T2 for each species
end
Cv_2 = Y_reac*Cv_2i';                                                       %Volumetric Heat Capacity of Mixture [J/KgK]

%INDICES WILL HAVE TO BE CHANGED ONCE ETHANOL ADDED 
Q_lhv = (Y_reac)/(Y_reac(1))*h_0'-Y_prod/(Y_reac(1))*h_0';                  %Lower heating value of gasoline [J/kg] 
                                                                            %Found to be 4.3416e6 [J/kg] online
                                                                            
V3  = V2;                                                                   %Isochoric
Q23 = Q_lhv*m2*Y_reac(1);                                                   %Caclulate total heat gain during combustion
T3  = T2+Q23/(m2*Cv_2);                                                     %Assume instant heating up 
m3  = m2;                                                                    

P3 = (m3*R_reac*T3)/(V3);                                                   %Ideal Gas Law

i = i+1;                                                                    %Move forward one dt, so not instantaneous

p_ideal(i) = P3;
T_ideal(i) = T3;
m_ideal(i) = m3;
%% Ideal Power Stroke 

for t = [time_expansion+2*dt:dt:time_valve_exhaust]
    i = i+1;                                                                %Keep updating same index as last time
    
    %Calculate the volumetric heat capacity for each molecule
    for j = [1:length(Species)]
        Cvj(j) = CvNasa(T_ideal(i-1),Species(j));
    end
    Cv = Y_prod*Cvj';                                                       %Calculate total heat capacity using mass fractions

    m_ideal(i) = m_ideal(i-1);                                              %Mass stays constant since valves are closed

    %Euler Forward Method:
    dTdt = -(R_prod*T_ideal(i-1))/(V_ideal(i-1)*Cv)*dVdt(theta_crank(i),dthetadt,B,r,L);
    T_ideal(i)= T_ideal(i-1)+dt*dTdt;

    p_ideal(i)= (m_ideal(i)*R_prod*T_ideal(i))/(V_ideal(i));                %Ideal Gas Law
end

P4 = p_ideal(i);
V4 = V_ideal(i);
T4 = T_ideal(i);
m4 = m_ideal(i);

%% Ideal Valve Exhaust

i = i+1;
p_ideal(i) = Pref;
T_ideal(i) = Tref;
m_ideal(i) = (p_ideal(i)*V_ideal(i))/(R_prod*T_ideal(i));

Q45 = (Tref - T4)*Cv*m4;                                            %Caclulate total heat loss during exhaust
%% Ideal Exhaust Stroke

for t = [time_exhaust_stroke+2*dt:dt:t_end]
    i = i+1;
    p_ideal(i) = Pref;                                                      %Isobaric
    T_ideal(i) = Tref;                                                      %Isothermal
    m_ideal(i) = (p_ideal(i)*V_ideal(i))/(T_ideal(i)*R_prod);               %Mass intake calculated using ideal gas law
end

%Work, Power & Efficiency Calculations
W_index = (cumtrapz(V_ideal,p_ideal));                                      %Index of work done throughout cycle [J]
W = W_index(end);                                                           %Total work done per cyycle [J]


P = W*w;                                                                    %Power of engine [W]

n_W = abs(W/Q23)                                                            %Efficiency using work and heat intake
n_Q = 1 - abs(Q45/Q23)                                                      %Efficiency using heat intake and outtake
n_C = 1 - abs(Tref/T3)                                                      %Maximum carnot efficiency 

%% Figures

figure 
hold on 
set(gca,'FontSize',18) 
title('Idealized Cycle')
xlabel('Volume [m^3]')
ylabel('Pressure [bar]')
xlim([0 inf])
ylim([0 inf])
grid on

plot(V_ideal,p_ideal.*1e-5)

hold off

figure 
semilogy(V_ideal,p_ideal)
hold on 
set(gca,'FontSize',18) 
title('Idealized Cycle')
xlabel('Volume [m^3]')
ylabel('Pressure [Pa]')

xlim([0 inf])
ylim([0 inf])

grid on 

hold off

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
