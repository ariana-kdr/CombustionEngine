clear all;close all;clc
%% add general to matlab path 
addpath('General');
%%
DataDir='C:\Users\20202832\Downloads\Training';  % The directory with the files. (Make sure you have backups!)
ColumnOrder={'time','Encoder','Sensor'};
%% Constants
rpm = 3000; % RPM
t_r = 60/rpm; % Period
omega = (2*pi*rpm)/60; % Angular velocity
%% Converting and Importing the Data
Files=dir(fullfile(DataDir,'0 eth no load.txt'));
nFiles=length(Files);                                                       % dir gives a directory listing, only *.txt files in this case
h=waitbar(0,'Converting');
for i=1:nFiles
    waitbar(i/nFiles,h);
    fname       = Files(i).name;                                            % Take a name from the list
    curfilename = fullfile(DataDir,fname);                                  % Create the full name
    comma2point_overwrite( curfilename);                                    % Convert comma to dot decimal separator, this is extremly fast!!
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
    P = ((V./V_s)-0.115)/0.00385; % Conversion from voltage to pressure
    
    % Input of mean values
    prompt2 = 'Data no- [1], half- [2] or full-load [3]? ';
    dlgtitle2 = 'Requested input';
    dlg_default2 = {'1'};
    dlg2 = inputdlg(prompt2, dlgtitle2, 1 , dlg_default2);     
    in = str2double(dlg2(1));
        if in == 1
            sub = 75;
        elseif in == 2
            sub = 60;
        elseif in == 3
            sub = 50;    
    end
    P = movmean(P,sub); % Remove noise
    nCyc = (length(RevEnd - 1))/2; % Number of cycles
        
    % Input of amount of wanted cycles
    prompt2 = 'How many cycles do you want? [1-10] ';
    dlgtitle2 = 'Requested input';
    dlg_default2 = {'1'};
    dlg2 = inputdlg(prompt2, dlgtitle2, 1 , dlg_default2);     
    out = str2double(dlg2(1));
    for i = 1:out
        sub = nCyc-i;
    end

    %Drift correction
    P_amb = 1.012;      %Ambient pressure according to wheather report
    P_exh = P(RevEnd);      %Exhaust pressure at RevEnd
    T_exh = [ones(length(t(RevEnd)),1) t(RevEnd)];
    Reg = T_exh\P_exh;      %Regression between pressure and time

    P_drift = P + (P_amb - Reg(1,:)) - Reg(2,:) * t;

    figure()
    hold on
for iCyc = 1:nCyc-sub
    Cycle_s = 2*iCyc; % Start of one cycle
    Cycle_e = Cycle_s+2; % End of one cycle

    t_c = t(RevEnd(Cycle_s):RevEnd(Cycle_e)); % Time of one cycle
    t_c = t_c - t_c(1);
    
    P_c = P_drift(RevEnd(Cycle_s):RevEnd(Cycle_e)); % Pressure of the specific cycle(s)
    omega_m = (4*pi)/max(t_c);              % Determined angular velocity for Ca
    Ca_off = 1.4476;                        %Crank angle offset found according to ignition
    
    Ca = omega_m*t_c - Ca_off;                       % Crank angle

    [~,P_max(iCyc)] = max(P_drift(RevEnd(Cycle_s):RevEnd(Cycle_e)));
    P_max_rel(iCyc) = P_max(iCyc)/length(P_drift(RevEnd(Cycle_s):RevEnd(Cycle_e)));

    V_c = V_cyl(Ca); % Volume using the function script
    plot(V_c,P_c);
    legend('Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5','Cycle 6','Cycle 7','Cycle 8','Cycle 9','Cycle 10');
    title('pV-diagram of experimental data of engine in full-load');
    xlabel('Volume [m^3]');
    ylabel('Pressure [bar]');
    
    W_index = cumtrapz(V_c,P_c);
    W = W_index(end)*10^5;
end
end
