%% General

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
%% Idealized Cycle
% Assumptions:
% - Both intake and exhaust gases are _ideal gases_ and are at _atmospheric 
%pressure_
% - _Adiabatic_ compression
% - _Isochoric_ and _Instantaneous_ combustion
% - _Adiabatic_ expansion
% - _Isochoric_ exhaust 
% - _Isothermal/Isobaric_ exhaust/intake stroke

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
                                                                            
V3 = V2;                                                                    %Isochoric
T3 = T2+(Q_lhv*m2*Y_reac(1))/(m2*Cv_2);                                     %Assume instant heating up 
m3 = m2;                                                                    

p3 = (m3*R_prod*T3)/(V3);                                                   %Ideal Gas Law

i = i+1;                                                                    %Move forward one dt, so not instantaneous

p_ideal(i) = p3;
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
%% Ideal Exhaust Stroke

for t = [time_exhaust_stroke+2*dt:dt:t_end]
    i = i+1;
    p_ideal(i) = Pref;                                                      %Isobaric
    T_ideal(i) = Tref;                                                      %Isothermal
    m_ideal(i) = (p_ideal(i)*V_ideal(i))/(T_ideal(i)*R_prod);               %Mass intake calculated using ideal gas law
end

%% Figures

figure 
hold on 
set(gca,'FontSize',18) 
title('Idealized Cycle')
xlabel('Volume [m^3]')
ylabel('Pressure [bar]')
xlim([0 inf])
ylim([0 150])
grid on

plot(V_ideal,p_ideal.*1e-5)

hold off

%% Real Cycle 
%Currently implemented:
% - Heat loss during: compression, ignition and expansion
% - Finite rate of combustion
% - Valve timing for intake

%Remaining:
% - Heat loss in all other cycles
% - Friction
% - Pumping/Negative Cycle
% - Valve timing for exhaust

%Create empty matrices for each state property
p_real = [];                                                                %Real cycle pressure array [Pa]                                                       
V_real = V_ideal;                                                           %Real cycle volume array [m^3]
T_real = [];                                                                %Real cycle temperature array [K]
m_real = [];                                                                %Real cycle mass array [kg]

%Parameters 
method_heat_coeff = 'Honenbe';                                              %Use either Woschni or Honenberg relation
t_w = 30e-3;                                                                %Estimate of thickness wall surrounding cylinder [m]
k = 14.4;                                                                   %Conductive heat transfer coefficient [W/mK]
ignition_offset = 25;                                                       %Degrees before TDC where ignition starts

%Set the angle for the start position of each of the stages
theta_intake_stroke = 0;
theta_compression = 240;                                                    %Using measurement data
theta_ignition = 360-ignition_offset;
theta_expansion = 360;
theta_valve_exhaust = 540;
theta_exhaust_stroke = 540;

%Find the index in the crank angle vector array for which it is at this
%angle
index_intake_stroke = max(find(theta_crank_degrees==theta_intake_stroke));
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

%% Real Intake Stroke
%Still the same as ideal, have to add negative work

i = 0;                                                                      %Index Value
for t = [time_intake_stroke:dt:time_compression]                            %Start to end time
    i = i+1;
    p_real(i) = Pref;                                                       %Isobaric
    T_real(i) = Tref;                                                       %Isothermal
    m_real(i) = (p_real(i)*V_real(i))/(T_real(i)*R_reac);                   %Mass intake calculated using ideal gas law
end
%% Real Compression
%No longer adiabatic, heat loss to surroundings is implemented

for t = [time_compression+dt:dt:time_ignition]
    i = i+1;                                                                %Keep updating same index as last time
    
    %Calculate the volumetric heat capacity for each molecule
    for j = [1:length(Species)]
        Cvj(j) = CvNasa(T_real(i-1),Species(j));
    end
    Cv = Y_reac*Cvj';                                                       %Calculate total heat capacity using mass fractions

    m_real(i) = m_real(i-1);                                                %Mass stays constant since valves are closed

    %Calculate thermal resistance depending on method selected
    if method_heat_coeff == 'Woschni'
        alpha = 1/(1/Woschni(B,p_real(i-1),T_real(i-1),S)+t_w/k);           %Thermal Resistance
    else 
        alpha = 1/(1/Honenberg(p_real(i-1),T_real(i-1),V_real(i-1),S)+t_w/k); %Total Thermal Resistance
    end
    
    %Euler Forward Method:   
    dTdt = -(alpha*Acyl(theta_crank(i-1),B,r,L))/(m_real(i)*Cv)*(T_real(i-1)-Tref)-(R_reac*T_real(i-1))/(V_real(i-1)*Cv)*dVdt(theta_crank(i),dthetadt,B,r,L);
    T_real(i)= T_real(i-1)+dt*dTdt;

    p_real(i)= (m_real(i)*R_reac*T_real(i))/(V_real(i));                    %Ideal Gas Law
end
%% Real Ignition

%INDICES WILL HAVE TO BE CHANGED ONCE ETHANOL ADDED 
Q_lhv = (Y_reac)/(Y_reac(1))*h_0'-Y_prod/(Y_reac(1))*h_0';                  %Lower heating value of gasoline [J/kg] 
                                                                            %Found to be 4.3416e6 [J/kg] online
                                                                            
n = 3;                                                                      %Wiebe Form Factor (experimental)
a = 5;                                                                      %Wiebe Efficiency Factor (expirimental)
theta_s = theta_ignition*pi/180;                                            %Start angle of combustion [rad]
theta_d = ignition_offset*pi/180;                                           %Combustion duration [rad]
                                                                            
for t = [time_ignition+dt:dt:time_expansion]
    i = i+1;
    %Calculate lower heating value
    for j = [1:length(Species)]
        h_0(j) = HNasa(Tref,Species(j));                                    %Calculate enthalpy of formation for each species
        Cv_2i(j) = CvNasa(T_real(i-1),Species(j));                          %Volumetric Heat Capacity at T for each species
    end
    Cv_2_reac = Y_reac*Cv_2i';                                              %Volumetric Heat Capacity of Reactants [J/KgK]
    Cv_2_prod = Y_prod*Cv_2i';                                              %Volumetric Heat Capacity of Products [J/KgK]
    
    %Calculate the Heat Capacity and gas constant of the mixture using the
    %fraction of mass that has combusted (Wiebe Function)
    Cv_2 = (1-Wiebe(n,a,theta_crank(i-1),theta_s,theta_d))*Cv_2_reac+Wiebe(n,a,theta_crank(i-1),theta_s,theta_d)*Cv_2_prod;
    R = (1-Wiebe(n,a,theta_crank(i-1),theta_s,theta_d))*R_reac+Wiebe(n,a,theta_crank(i-1),theta_s,theta_d)*R_prod;
    
    m_real(i) = m_real(i-1);                                                %Mass stays constant since valves are closed
    m_f = Y_reac(1)*m_real(i);                                              %Mass of fuel using mass composition [kg]
    
     %Calculate thermal resistance depending on method selected
    if method_heat_coeff == 'Woschni'
        alpha = 1/(1/Woschni(B,p_real(i-1),T_real(i-1),S)+t_w/k);
    else 
        alpha = 1/(1/Honenberg(p_real(i-1),T_real(i-1),V_real(i-1),S)+t_w/k);
    end
    
    %Euler Forward Method
    dTdt = (Q_lhv)/(m_real(i-1)*Cv_2)*dmdtheta(m_f,n,a,theta_crank(i-1),theta_s,theta_d)*dthetadt-(alpha*Acyl(theta_crank(i-1),B,r,L))/(m_real(i)*Cv)*(T_real(i-1)-Tref)-(R*T_real(i-1))/(V_real(i-1)*Cv)*dVdt(theta_crank(i),dthetadt,B,r,L);
    T_real(i)= T_real(i-1)+dt*dTdt;
    
    p_real(i)= (m_real(i)*R*T_real(i))/(V_real(i));                         %Ideal Gas Law
end
%% Real Power Stroke

for t = [time_expansion+dt:dt:time_valve_exhaust]
    i = i+1;                                                                %Keep updating same index as last time
    
    %Calculate the volumetric heat capacity for each molecule
    for j = [1:length(Species)]
        Cvj(j) = CvNasa(T_real(i-1),Species(j));
    end
    Cv = Y_prod*Cvj';                                                       %Calculate total heat capacity using mass fractions

    m_real(i) = m_real(i-1);                                                %Mass stays constant since valves are closed
    
    %Calculate thermal resistance depending on method selected
    if method_heat_coeff == 'Woschni'
        alpha = 1/(1/Woschni(B,p_real(i-1),T_real(i-1),S)+t_w/k);
    else 
        alpha = 1/(1/Honenberg(p_real(i-1),T_real(i-1),V_real(i-1),S)+t_w/k);
    end

    %Euler Forward Method:
    dTdt = -(alpha*Acyl(theta_crank(i-1),B,r,L))/(m_real(i)*Cv)*(T_real(i-1)-Tref)-(R_prod*T_real(i-1))/(V_real(i-1)*Cv)*dVdt(theta_crank(i),dthetadt,B,r,L);
    T_real(i)= T_real(i-1)+dt*dTdt;

    p_real(i)= (m_real(i)*R_prod*T_real(i))/(V_real(i));                    %Ideal Gas Law
end

i_power = i;
%% Real Valve Exhaust
%Still needs to be changed for real cycle

i = i+1;
p_real(i) = Pref;
T_real(i) = Tref;
m_real(i) = (p_real(i)*V_real(i))/(R_prod*T_real(i));
%% Real Exhaust Stroke
%Still needs to be changed for real cycle

for t = [time_exhaust_stroke+2*dt:dt:t_end]
    i = i+1;
    p_real(i) = Pref;                                                       %Isobaric
    T_real(i) = Tref;                                                       %Isothermal
    m_real(i) = (p_real(i)*V_real(i))/(T_real(i)*R_prod);                   %Mass intake calculated using ideal gas law
end
%% Figures

figure 
hold on 
set(gca,'FontSize',18) 
title('Real Cycle')
xlabel('Volume [m^3]')
ylabel('Pressure [bar]')
xlim([0 inf])
ylim([0 100])
grid on

plot(V_real,p_real.*1e-5)

hold off

figure 
hold on 
set(gca,'FontSize',18) 
title('Comparison')
xlabel('Volume [m^3]')
ylabel('Pressure [bar]')
xlim([0 inf])
ylim([0 150])
grid on

plot(V_ideal,p_ideal.*1e-5)
plot(V_real,p_real.*1e-5)

legend('Ideal Cycle', 'Real Cycle')

hold off
%% Functions

function [Vcyl] = Vcyl(theta,B,r,L,V_c)
%This function is used to calculate the volume of the cylinder 
%   See lecture for derivation of formula
    x = r*cos(theta)+sqrt(L^2-r^2*sin(theta)^2);                            %Distance from center of crankshaft to center of piston head [m]
    d = L+r-x;                                                              %Distance of top of piston head to its maximum height [m]

    Vcyl = pi/4*B^2*d+V_c;                                                  %Cylinder Volume [m^3]
end

function [Acyl] = Acyl(theta,B,r,L)
%Calculate cylinder area 
    x = r*cos(theta)+sqrt(L^2-r^2*sin(theta)^2);                            %Distance from center of crankshaft to center of piston head [m]
    d = L+r-x;                                                              %Distance of top of piston head to its maximum height [m]
    
    Acyl = pi*B*d+2*pi/4*B^2;
end

function [dVdt] = dVdt(theta,dthetadt,B,r,L)
%Calculate derivative of volume with respect to time for given theta
    dVdt = pi/4*B^2*(dthetadt*r*sin(theta)+0.5*(dthetadt*r^2*sin(2*theta))/(sqrt(L^2-0.5*r^2*(1-cos(2*theta)))));
end

function [alpha] = Woschni(B,p,T,S)
%Calculate heat transfer coefficient using Woschni Relation
    v_p = 2*(3000)/(60)*S;                                                  %Average piston speed (used in head coeff relations)
    alpha = 3.26*B^(-0.2)*p^(0.8)*T^(-0.55)*(2.28*v_p)^(0.8);
end

function [alpha] = Honenberg(p,T,V,S)
%Calculate heat transfer coefficient using Honenberg Relation
    v_p = 2*(3000)/(60)*S;                                                  %Average piston speed (used in head coeff relations)
    alpha = 3.26*p^(0.8)*T^(-0.4)*V^(-0.06)*(v_p+1.4)^0.8;
end

function [x_frac] = Wiebe(n,a,theta,theta_s,theta_d)
%Fraction of mass combusted using Wiebe Function
    x_frac = 1-exp(-a*((theta-theta_s)/theta_d)^n);
end

function [dmdtheta] = dmdtheta(m_f,n,a,theta,theta_s,theta_d)
%Wiebe function from project manual 
    dmdtheta = m_f*n*a*(((theta-theta_s)/theta_d)^(n-1))/(theta_d)*exp(-a*((theta-theta_s)/theta_d)^n);
end
