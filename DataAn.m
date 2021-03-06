clear all;close all;clc
%% add general to matlab path

addpath('General');
%%
DataDir='C:\Users\20202587\Documents\TU Eindhoven\Jaar 2\4GB10 Combustion Engine\Experiment 2\Second attempt\Training 2';  % The directory with the files. (Make sure you have backups!)
ColumnOrder={'time','Encoder','Sensor'};
%% Constants

rpm = 3000; % RPM
t_r = 60/rpm; % Period
omega = (2*pi*rpm)/60; % Angular velocity
%% Converting and Importing the Data

Dim_grid = 3;                                                               % set size of subplot grid
FDir= dir(DataDir); 
%%Files = FDir(4:4,:);
Files = FDir(3:2+Dim_grid^2,:);                                                      % remove '.' and '..', limit the amount of plots to 9
nFiles=length(Files);                                                       % dir gives a directory listing, only *.txt files in this case
h=waitbar(0,'Converting');
for i=1:nFiles
    waitbar(i/nFiles,h);
    fname       = Files(i).name;                                            % Take a name from the list
    curfilename = fullfile(DataDir,fname);                                  % Create the full name
    %comma2point_overwrite(curfilename);                                    % Convert comma to dot decimal separator, this is extremly fast!!
end
delete(h)

for iFiles=1:nFiles
    fname       = Files(iFiles).name;                                       % Take a name from the list
    curfilename = fullfile(DataDir,fname);                                  % Create the full name
    Data        = ImportData4GB10(curfilename,ColumnOrder);                 % Read the data. Type help ImportData4GB10
    % store it for each case. Yes a struct that contains a struct and other
    % things. Why not?
    Case(iFiles).Data     = Data;
    Case(iFiles).filename = fname;
    Case(iFiles).DataDir  = DataDir;
    % preamble, put data in easy to use arrays
    t      = Data.t;
    p      = Data.pulse;
    V      = Data.Volt;
    RevEnd = Data.RevEnds;
    V_s = 5; % Maximum voltage
    P = ((V./V_s)-0.115)/(0.00385*4); % Conversion from voltage to pressure

    % Input of mean values
    if contains(fname,'no') == true
        Q = 75;
    elseif contains(fname,'half') == true
        Q = 60;
    elseif contains(fname,'full') == true
        Q = 50;
    else
        Q = 60;
    end
    P = movmean(P,Q); % Remove noise
    nCyc = (length(RevEnd) - 1)/2; % Number of cycles

    % Input of amount of wanted cycles
    out = 10;
    for i = 1:out
        Q = nCyc-i;
    end

    %Drift correction
    P_amb = 1.012;                                                          %Ambient pressure according to wheather report
    %%Ca_off = 1.4779;                                                        %Crank angle offset found according to ignition
    Ca_off = 1.445; %New TDC
    P_exh = P(RevEnd - round(Ca_off/(4 * pi)*mean(diff(RevEnd))));          %Exhaust pressure at end of exhaust
    T_exh = [ones(length(t(RevEnd)),1) t(RevEnd)];
    Reg = T_exh\P_exh;                                                      %Regression between pressure and time

    P_drift = P + (P_amb - Reg(1,:)) - Reg(2,:) * t;

    subplot(Dim_grid,Dim_grid,iFiles)              %Create a 3x3 subplot
    title(fname(1:end-4));
    xlabel('Volume [m^3]');
    ylabel('Pressure [bar]');
    %legend('Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5','Cycle 6','Cycle 7','Cycle 8','Cycle 9','Cycle 10');
    hold on
    for iCyc = 1:nCyc-Q
        Cycle_s = 2*iCyc; % Start of one cycle
        Cycle_e = Cycle_s+2; % End of one cycle

        t_c = t(RevEnd(Cycle_s):RevEnd(Cycle_e)); % Time of one cycle
        t_c = t_c - t_c(1);

        P_c = P_drift(RevEnd(Cycle_s):RevEnd(Cycle_e)); % Pressure of the specific cycle(s)
        omega_m = (4*pi)/max(t_c);              % Determined angular velocity for Ca

        Ca = omega_m*t_c - Ca_off;                       % Crank angle

        [~,P_max(iCyc)] = max(P(RevEnd(Cycle_s):RevEnd(Cycle_e)));
        P_max_rel(iCyc) = P_max(iCyc)/length(P(RevEnd(Cycle_s):RevEnd(Cycle_e)));

        V_c = V_cyl(Ca); % Volume using the function script
        plot(V_c,P_c);
    end
end

%%
%% Work, Power & Efficiency Calculations
%W_index = (cumtrapz(V_ci_total, P_ci_total));                                        %Index of work done throughout cycle [J]
%W = W_index(end);                                                           %Total work done per cyycle [J]

W_total = V_ci_total .* P_ci_total;
W_sum = sum(W_total(2:end,:));

P = W_sum*3000/60;                                                                    %Power of engine [W]

return

        %W = trapz(V_c,P_c)*10^5;
        Q_lhv = [4.2894e7 4.3525e7 4.4227e7 4.5012e7];
    if contains(fname,'0 eth') == true
        Q = Q_lhv(1,1);
    elseif contains(fname,'5 eth') == true
        Q = Q_lhv(1,2);
    elseif contains(fname,'10 eth') == true
        Q = Q_lhv(1,3);
    elseif contains(fname,'15 eth') == true
        Q = Q_lhv(1,4);
    end
    m_dot = importdata('C:\Users\20192807\Desktop\Combustion\Triaining\matlab\Data\massflow.txt');
    if contains(fname,'no') == true
        m_flow = m_dot(:,1);  % No load flows
    elseif contains(fname,'half') == true
        m_flow = m_dot(:,2);
    elseif contains(fname,'full') == true
        m_flow = m_dot(:,3);
    end
    eta = W./(Q_lhv*m_flow);
