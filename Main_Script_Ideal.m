%% General


% Known Issues:
% - Indexing in for loop is not perfect yet
% - Since volume is a function of time ignition and valve exhaust are not
% instantaneous
% - Q_lhv behaves strangely when ethanol is added to the blend, and in
% general it has to be reconsidered due to the signs in the code not
% matching the signs in the slides 


clear all;
close all;
clc;


% Add directory of Nasa routines to Matlab-path
addpath('General/Nasa');    


% Load Nasa Database
TdataBase = fullfile('General\Nasa','NasaThermalDatabase.mat');
load(TdataBase);


% Constants 
global Runiv Pref Tref                                                      %Set global variables (so they cannot be changed)
Runiv = 8.314472;                                                           %Universal Gas Constant[J/K*mol]
Pref = 1.01235e5;                                                           %Reference Pressure (1atm) [Pa]
Tref = 298.15;                                                              %Reference Temperature [K]


% Loading chemical properties 
Fuel1 = 'Gasoline';                                                          % Make fuel its own variable for easy changes
Fuel2 = 'C2H5OH'; 
Species = Sp(myfind({Sp.Name},{Fuel1, Fuel2,'O2','CO2','H2O','N2'}));               % Subselection of species in database
NSpecies = length(Species);                                                 % Create variable for length of Species array for later use 
Mi = [Species.Mass];                                                        % Array of molecular masses for each component


%Dimensions (using dimensions found online)
B = 68e-3;                                                                  %Bore [m]
r = 27e-3;                                                                  %Radius of crank shaft [m]
S = 54e-3;                                                                  %Stroke [m]
L = 84.7e-3;                                                                %Length of rod [m]
r_c = 8.5;                                                                  %Compression Ratio [m]
V_d = pi/4*B^2*S;                                                           %Displacement Volume [m^3]
V_c = V_d/(r_c-1);                                                          %Clearance volume [m^3]
%% ideal Cycle 
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

steps = 4000;                                                                           %Number of steps to be used for cycle. 10000 seems to be a good balance between precision and speed

theta_crank = zeros(1,steps);                                                          %ideal cycle crank angle array [rad]
p_ideal = zeros(1,steps);                                                               %ideal cycle pressure array [Pa]                                                       
V_ideal = zeros(1,steps);                                                               %ideal cycle volume array [m^3]
T_ideal = zeros(1,steps);                                                               %ideal cycle temperature array [K]
m_ideal = zeros(1,steps);                                                               %ideal cycle mass array [kg]
s_ideal = zeros(1,steps);                                                               %ideal cycle entropy 
Q_in_ideal = zeros(1,steps);                                                           %Total heat in in ideal cycle [J]

%Create struct to combine all the results of the model

eth = [0,0.05,0.10,0.15,0.2,0.4,0.6,0.8];                                                                  % ethanol fractions (can be any string of numbers between 0 and 1)
n_eth = length(eth);                                                               % number of ethanol blends considered


Results = struct('Blend',[],'Percentage', [], 'Volume',[],'Pressure',[],'Temperature',[],'Entropy',[],'Work',[],'Heat',[], 'Efficiency', [], 'Power', []); 

for t = 1:n_eth
    Results(t).Blend = 100*eth(t) + "% ethanol";
end


                                                                          
%Time Values
t_end = 2*(60/3000);                                                        %End time [s] (one rotation is 0.01 sec however, one cycle is two revolutions)
dt = t_end/(steps-1);                                                       %Time step to solve differential equations [s]
dthetadt = 3000/60*2*pi;                                                    %Convert speed of engine to [rad/s]
time = [0:dt:t_end];                                                        %Time array [s]


%Parameters 
method_heat_coeff = 'Honenbe';                                              %Use either Woschni or Honenberg relation
t_w = 30e-3;                                                                %Estimate of thickness wall surrounding cylinder [m]
k = 14.4;                                                                   %Conductive heat transfer coefficient [W/mK]
ignition_offset = 25;                                                       %Degrees before TDC where ignition starts
ignition_duration = 40;


%Create vector array of crank angle and volume as a function of time
i = 0;                                                                      %Index Value
for t = [0:dt:t_end]
    i = i+1;                                                                %Update Index
    theta_crank(i) = t*dthetadt;                                            %Calculate angle for each time instance
    V_ideal(i) = Vcyl(theta_crank(i),B,r,L,V_c);                             %Calculate volume for angle using volume function
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
index_intake_stroke = find(theta_crank_degrees==theta_intake_stroke, 1, 'last' );
%Since angle is never perfectly 180 take max index where it is lower and
%add one (cheap workaround)
index_compression = find(theta_crank_degrees<theta_compression, 1, 'last' )+1;
index_ignition = find(theta_crank_degrees<theta_ignition, 1, 'last' )+1;
index_expansion = find(theta_crank_degrees<theta_expansion, 1, 'last' )+1;
index_valve_exhaust = find(theta_crank_degrees<theta_valve_exhaust, 1, 'last' )+1;
index_exhaust_stroke = find(theta_crank_degrees<theta_exhaust_stroke, 1, 'last' )+1;

%Find the time instant for each start of the strokes
time_intake_stroke = time(index_intake_stroke);
time_compression = time(index_compression);
time_ignition = time(index_ignition);
time_expansion = time(index_expansion);
time_valve_exhaust = time(index_valve_exhaust);
time_exhaust_stroke = time(index_exhaust_stroke);

p = 0;
for e = eth
p = p + 1;

%Calculate Fuel Properties
x = (1 - e)*Species(1).Elcomp(3) + e*Species(2).Elcomp(3);                  % Summing the moles of given ethanol/gasoline ratio for carbon
y = (1 - e)*Species(1).Elcomp(2) + e*Species(2).Elcomp(2);                  % Summing the moles of given ethanol/gasoline ratio for hydrogen
z = (1 - e)*Species(1).Elcomp(1) + e*Species(2).Elcomp(1);                  % Summing the moles of given ethanol/gasoline ratio for oxygen

a = x + y/4 - z/2;


% Number of moles in chemical reaction
N_reac = [1-e e a 0 0 a*3.76];
N_prod = [0 0 0 x y/2 a*3.76];


% Mole Fractions using balanced chemical equation
X_reac = N_reac/sum(N_reac);
X_prod = N_prod/sum(N_prod);


%Convert to mass fractions
Y_reac = (X_reac.*Mi)/(X_reac*Mi');
Y_prod = (X_prod.*Mi)/(X_prod*Mi');


%Molar mass
M_reac = X_reac*Mi';
M_prod = X_prod*Mi';


%Gas constants
R_reac = Runiv/M_reac;
R_prod = Runiv/M_prod;
%% ideal Intake Stroke
%Still the same as ideal, have to add negative work

%% Ideal Intake Stroke

i = 0;                                                                      %Index Value
for t = [time_intake_stroke:dt:time_compression]                            %Start to end time
    i = i+1;
    p_ideal(i) = Pref;                                                      %Isobaric
    T_ideal(i) = Tref;                                                      %Isothermal
    
    %Calculate entropy
    for j = [1:length(Species)]
        sj(j)= SNasa(T_ideal(i),Species(j));
    end
    s_ideal(i)= Y_reac*sj' - R_reac * (log(p_ideal(i))-log(Pref));

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

    %Calculate entropy
    for j = [1:length(Species)]
        sj(j)= SNasa(T_ideal(i),Species(j));
    end
    s_ideal(i) = Y_reac*sj' - R_reac*(log(p_ideal(i))-log(Pref));
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

%Lower heating value of gasoline + ethanol[J/kg] 
%Found to be 43.416e6 [J/kg] online for gasoline     
Q_lhv = Y_reac/(Y_reac(1) + Y_reac(2))*h_0' - Y_prod/(Y_reac(1) + Y_reac(2))*h_0';  
                         
V3 = V2;                                                                    %Isochoric
Q23 = Q_lhv*m2*(Y_reac(1));                                     %Calculate total heat gain during combustion
T3 = T2 + Q23/(m2*Cv_2);                                                    %Assume instant heating up 
m3 = m2;                                                                    

p3 = (m3*R_reac*T3)/(V3);                                                   %Ideal Gas Law

i = i+1;                                                                    %Move forward one dt, so not instantaneous

%Calculate entropy

for j = [1:length(Species)]
        sj(j)= SNasa(T_ideal(i),Species(j));
end
s_ideal(i) = Y_reac*sj' - R_reac * (log(p_ideal(i))-log(Pref));

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

    %Calculate entropy
    for j = [1:length(Species)]
        sj(j)= SNasa(T_ideal(i),Species(j));
    end
    s_ideal(i)= Y_prod*sj' - R_prod * (log(p_ideal(i))-log(Pref));
end

P4 = p_ideal(i);
V4 = V_ideal(i);
T4 = T_ideal(i);
m4 = m_ideal(i);


%% Ideal Valve Exhaust

i = i+1;
p_ideal(i) = Pref;
T_ideal(i) = Tref;

%Calculate entropy
for j = [1:length(Species)]
    sj(j)= SNasa(T_ideal(i),Species(j));
end
s_ideal(i)= Y_prod*sj' - R_prod * (log(p_ideal(i))-log(Pref));

m_ideal(i) = (p_ideal(i)*V_ideal(i))/(R_prod*T_ideal(i));
Q45 = (Tref - T4)*Cv*m4;                                                    %Calculate total heat loss during exhaust

%% Ideal Exhaust Stroke

for t = [time_exhaust_stroke+2*dt:dt:t_end]
    i = i+1;
    p_ideal(i) = Pref;                                                      %Isobaric
    T_ideal(i) = Tref;                                                      %Isothermal

    %Calculate entropy
    for j = [1:length(Species)]
        sj(j)= SNasa(T_ideal(i),Species(j));
    end
    s_ideal(i)= Y_prod*sj' - R_prod * (log(p_ideal(i))-log(Pref));

    m_ideal(i) = (p_ideal(i)*V_ideal(i))/(T_ideal(i)*R_prod);               %Mass intake calculated using ideal gas law
end

%% Work, Power & Efficiency Calculations
W_index = (cumtrapz(V_ideal,p_ideal));                                        %Index of work done throughout cycle [J]
W = W_index(end);                                                           %Total work done per cyycle [J]




P = W*3000/120;                                                                    %Power of engine [W]

Q_in = Q23;                                                      %Heat added per cycle
n_W = abs(W/Q_in);                                                           %Efficiency using work and heat intake                                                  
n_C = 1 - abs(min(T_ideal)/max(T_ideal));                                                     %Maximum carnot efficiency

%% Record results in struc
Results(p).Percentage = e;
Results(p).Volume= V_ideal;
Results(p).Pressure = p_ideal;
Results(p).Temperature = T_ideal;
Results(p).Entropy = s_ideal;
Results(p).Work = W;
Results(p).Heat = Q_in;
Results(p).Efficiency = n_W;
Results(p).Power = P;



end




%% Figures


figure 
hold on 
set(gca,'FontSize',18) 
title('P-v Diagram for the ideal Cycle')
xlabel('Volume [m^3]')
ylabel('Pressure [bar]')

xlim([0 3e-4])
ylim([0 inf])

grid on

for p = 1:n_eth
plot(Results(p).Volume,Results(p).Pressure.*1e-5)
end

legend(100*eth + "% Ethanol")

hold off


figure
hold on
set(gca,'FontSize',18) 
title('Efficiency over Ethanol Percentage')
xlabel('Ethanol (%)')
ylabel('Thermal Efficiency (%)')
ylim([0.4 0.7])
plot([Results.Percentage], [Results.Efficiency])
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
