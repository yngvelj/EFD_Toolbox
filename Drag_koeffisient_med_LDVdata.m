clear all;
Forcedata = load('ForceData.mat');
LDVdata = load('LDVdata.mat');
water = water_properties;
T = 20; %benchmark reynoldsnumber temperature
h = 60e-3; %Heigh of bump
rho = water.rho(T); %Density: 997.7998
u = 0.5:0.1:1; %Velocity vector
Rn = u*h/water.nu(T); %Benchmark reynoldsnumber
Ap = 24021.851e-6; %Projected area of the bump


for i = 20:7:length(Forcedata.Data)
    if(i ~= length(Forcedata.Data))
       for j = 1:7
            FData(i+j).name = Forcedata.Data(i+j).name;
            FData(i+j).data = Forcedata.Data(i+j).data.a2;
       end
    end
end

%7600:Force without bump
Data.F7600 = Forcedata.Data(7:16);
F7600.temp = 20;
%LDV Velocity
Data.LDV7600 = LDVdata.LDVdata(77:86);
LDV7600.vel = [Data.LDV7600.stats]; %temporarly variable for the mean velocity from LDV
F7600.Rn = [LDV7600.vel.mean]*h./water.nu(F7600.temp);%Calculate reynoldsnumber

%1000:Decay test
Data.Decay = Forcedata.Data(1:6);


%7000:Force measurements
Data.F7000 = Forcedata.Data(17:20);
% %LDV Velocity ... Finner ikke LDV-data...

%7430:Force
Data.F7430 = Forcedata.Data(21:27);
F7430.temp = 20.9;
%LDV Velocity
Data.LDV7430 = LDVdata.LDVdata(14:20);
LDV7430.vel = [Data.LDV7430.stats]; %temporarly variable for the mean velocity from LDV
F7430.Rn = [LDV7430.vel.mean]*h./water.nu(F7430.temp);%Calculate reynoldsnumber


%7440:Force
Data.F7440 = Forcedata.Data(28:34);
F7440.temp = 20.9;
%LDV Velocity
Data.LDV7440 = LDVdata.LDVdata(21:27);
LDV7440.vel = [Data.LDV7440.stats]; %temporarly variable for the mean velocity from LDV
F7440.Rn = [LDV7440.vel.mean]*h./water.nu(F7440.temp);%Calculate reynoldsnumber

%7450:Force
Data.F7450 = Forcedata.Data(35:41);
F7450.temp = 20.9;
%LDV Velocity
Data.LDV7450 = LDVdata.LDVdata(28:34);
LDV7450.vel = [Data.LDV7450.stats]; %temporarly variable for the mean velocity from LDV
F7450.Rn = [LDV7450.vel.mean]*h./water.nu(F7450.temp);%Calculate reynoldsnumber

%7460:Force
Data.F7460 = Forcedata.Data(42:48);
F7460.temp = 20.5;
%LDV Velocity
Data.LDV7460 = LDVdata.LDVdata(35:41);
LDV7460.vel =[Data.LDV7460.stats]; %temporarly variable for the mean velocity from LDV
F7460.Rn = [LDV7460.vel.mean]*h./water.nu(F7460.temp);%Calculate reynoldsnumber

%7470:Force
Data.F7470 = Forcedata.Data(49:55);
F7470.temp = 20.5;
%LDV Velocity
Data.LDV7470 = LDVdata.LDVdata(42:48);
LDV7470.vel = [Data.LDV7470.stats]; %temporarly variable for the velocitydata from LDV
F7470.Rn = [LDV7470.vel.mean]*h./water.nu(F7470.temp);%Calculate reynoldsnumber

%7480:Force
Data.F7480 = Forcedata.Data(56:62);
F7480.temp = 20.5;
%LDV Velocity
Data.LDV7480 = LDVdata.LDVdata(49:55);
LDV7480.vel = [Data.LDV7480.stats]; %temporarly variable for the mean velocity from LDV
F7480.Rn = [LDV7480.vel.mean]*h./water.nu(F7480.temp);%Calculate reynoldsnumber

%7490:Force
Data.F7490 = Forcedata.Data(63:69);
F7490.temp = 19.8;
%LDV Velocity
Data.LDV7490 = LDVdata.LDVdata(56:62);
LDV7490.vel = [Data.LDV7490.stats]; %temporarly variable for the mean velocity from LDV
F7490.Rn = [LDV7490.vel.mean]*h./water.nu(F7490.temp);%Calculate reynoldsnumber

%7500:Force
Data.F7500 = Forcedata.Data(70:76);
F7500.temp = 19.8;
%LDV Velocity
Data.LDV7500 = LDVdata.LDVdata(63:69);
LDV7500.vel = [Data.LDV7500.stats]; %temporarly variable for the mean velocity from LDV
F7500.Rn = [LDV7500.vel.mean]*h./water.nu(F7500.temp);%Calculate reynoldsnumber

%7510:Force
Data.F7510 = Forcedata.Data(77:83);
F7510.temp = 19.8;
%LDV Velocity
Data.LDV7510 = LDVdata.LDVdata(70:76);
LDV7510.vel = [Data.LDV7510.stats]; %temporarly variable for the mean velocity from LDV
F7510.Rn = [LDV7510.vel.mean]*h./water.nu(F7510.temp);%Calculate reynoldsnumber




%% Drag on the flat plate 

for i = 1:length(Data.F7600)
   F7600_mean(i) = mean(Data.F7600(i).data.a2(16).data);
   Vel7600_mean(i) = Data.LDV7600(i).stats.mean;
end

for i = 1:length(Data.F7600)
   F7600_mean(i) = mean(Data.F7600(i).data.a2(16).data);
   Vel7600_mean(i) = Data.LDV7600(i).stats.mean;
end
for i = 1:length(Data.F7600)
   F7600_mean(i) = mean(Data.F7600(i).data.a2(16).data);
   Vel7600_mean(i) = Data.LDV7600(i).stats.mean;
end
for i = 1:length(Data.F7600)
   F7600_mean(i) = mean(Data.F7600(i).data.a2(16).data);
   Vel7600_mean(i) = Data.LDV7600(i).stats.mean;
end
for i = 1:length(Data.F7600)
   F7600_mean(i) = mean(Data.F7600(i).data.a2(16).data);
   Vel7600_mean(i) = Data.LDV7600(i).stats.mean;
end


 %% ########################Drag of plate###########################
 %---------------------------------------------------------------
 %Calculate and interpolate the drag of the plate as a function of
 %Reynoldsnumber
 %Calculate the friction coefficient
B = 0.56;
L = 1.3;
S = 2*L*B;
t = 6e-3;
nu = water.nu(T); %Viscosity
Rex = Vel7600_mean.*L/nu; %Reynoldsnumber (length of plate)
Ap_plate = B*t;
%Drag from the struts:
H = 295e-3; %Height of the strut
t = 21e-3;
c = 100e-3;
Cf_struts = 0.075./(log10(Vel7600_mean*c/nu)+2).^2;
% DCf_struts = (111*(150*Vel7600_mean).^0.21-404).*Cf_struts.^2; %Roughness hadde lite eller ingen ting å si for resultatet
Cd_struts = 2*(Cf_struts)*(1+2*(t/c)+ 60*(t/c)^4); %2D Drag coefficient of the struts

Cd_plate = 1; %pressure drag coefficient for a square in infinite fluid Hoerner
F_plate_basedrag = Cd_plate*1/2*rho*Ap_plate*Vel7600_mean.^2;%Induced viscouse pressure drag on plate

F_struts = 1/2*rho*Vel7600_mean.^2*c.*Cd_struts*H*4; %Calculated drag from the four struts

F7600_plate = F7600_mean-F_struts-F_plate_basedrag; %Friction drag on plate
F7600_mean_intp = interp1([F7600.Rn],F7600_mean,Rn); %Drag from struts and plate without bump

%Friction coefficients
Cf = (F7600_plate./S)./(1/2*rho*Vel7600_mean.^2); %Calculated from experiment
Cf_blasius = 0.664./(sqrt(Rex)); %Blasius
Cf_Prandtl = 0.027./Rex.^1/7; %Prandtl
Cf_ittc = 0.075./(log10(Rex)-2).^2; %ITTC57

figure()
hold on;
ylabel('Friction coefficient')
xlabel('Rn [-]');
plot(Rex,Cf);
plot(Rex,Cf_blasius);
plot(Rex,Cf_Prandtl);
plot(Rex,Cf_ittc);
legend('C_f experiment','C_f Blasius','C_f Prandtl','C_f ITTC57');



%% ####################Drag analysis###########################
%------------------------------------------------------------
%Calculates the drag coefficient for the bump at 0 degree orientation
%Run numbers: 7430, 7470, 7480, 7510


% for i = 1:4 %Finner ikke LDV-data(har ikke inflow velocity..)
%     F7000_mean(i) = mean(Data.F7000(i).data.a2(16).data); 
% end
for i = 2:7
    F7430_mean(i) = mean(Data.F7430(i).data.a2(16).data);
    F7440_mean(i) = mean(Data.F7440(i).data.a2(16).data);
    F7470_mean(i) = mean(Data.F7470(i).data.a2(16).data);
    F7480_mean(i) = mean(Data.F7480(i).data.a2(16).data);
    F7510_mean(i) = mean(Data.F7510(i).data.a2(16).data);
    
end


%Run 7430 0degrees rotation drag
F7430_mean_intp = interp1([F7430.Rn],F7430_mean,Rn);
F7440_mean_intp = interp1([F7440.Rn],F7440_mean,Rn);
F7470_mean_intp = interp1([F7470.Rn],F7470_mean,Rn);
F7480_mean_intp = interp1([F7480.Rn],F7480_mean,Rn);
F7510_mean_intp = interp1([F7510.Rn],F7510_mean,Rn);

figure()
hold on;
title('Force vs Reynolds number');
xlabel('Rn [-]');
ylabel('Force [N]');
plot(F7430.Rn,F7430_mean);
plot(Rn,F7430_mean_intp);

plot(F7440.Rn,F7440_mean);
plot(Rn,F7440_mean_intp);

plot(F7470.Rn,F7470_mean);
plot(Rn,F7470_mean_intp);

plot(F7480.Rn,F7480_mean);
plot(Rn,F7480_mean_intp);

plot(F7510.Rn,F7510_mean);
plot(Rn,F7510_mean_intp);
legend('F7430 mean', 'F7430 intrp','F7440 mean', 'F7440 intrp','F7470 mean', 'F7470 intrp','F7480 mean', 'F7480 intrp','F7510 mean', 'F7510 intrp');

figure()
hold on;
title('Drag vs velocity');
xlabel('Inflow velocity [m/s]');
ylabel('Drag [N]');
plot([LDV7430.vel.mean],F7430_mean);
plot(u,F7430_mean_intp);

plot([LDV7440.vel.mean],F7440_mean);
plot(u,F7440_mean_intp);

plot([LDV7470.vel.mean],F7470_mean);
plot(u,F7470_mean_intp);

plot([LDV7480.vel.mean],F7480_mean);
plot(u,F7480_mean_intp);

plot([LDV7510.vel.mean],F7510_mean);
plot(u,F7510_mean_intp);
legend('F7430 mean', 'F7430 intrp','F7440 mean', 'F7440 intrp','F7470 mean', 'F7470 intrp','F7480 mean', 'F7480 intrp','F7510 mean', 'F7510 intrp');

figure()
hold on;
xlabel('runnr')
ylabel('velocity');
run = [0 1 2 3 4 5 6];
plot(run,[LDV7430.vel.mean]);
plot(run,[LDV7440.vel.mean]);
plot(run,[LDV7470.vel.mean]);
plot(run,[LDV7480.vel.mean]);
plot(run,[LDV7510.vel.mean]);
legend('7430','7470','7480','7510');


%Drag on the bump
F7430_bump = F7430_mean_intp-F7600_mean_intp;
F7440_bump = F7440_mean_intp-F7600_mean_intp;
F7470_bump = F7470_mean_intp-F7600_mean_intp;
F7480_bump = F7480_mean_intp-F7600_mean_intp;
F7510_bump = F7510_mean_intp-F7600_mean_intp;

%Drag coefficients
Cd_7430_bump = F7430_bump./(1/2*rho*u.^2*Ap); %Drag coefficient for one bump
Cd_7440_bump = F7440_bump./(1/2*rho*u.^2*Ap);
Cd_7470_bump = F7470_bump./(1/2*rho*u.^2*Ap);
Cd_7480_bump = F7480_bump./(1/2*rho*u.^2*Ap);
Cd_7510_bump = F7510_bump./(1/2*rho*u.^2*Ap);


figure()
hold on;
title('Drag coefficients BUMP 0 [deg] rotation');
ylabel('C_d [-]');
xlabel('Reynolds number [-]');
plot(Rn,Cd_7430_bump);
plot(Rn,Cd_7440_bump);
plot(Rn,Cd_7470_bump);
plot(Rn,Cd_7480_bump);
plot(Rn,Cd_7510_bump);
legend('7430','7440','7470','7480','7510');






%% ################### LIFT ANALYSIS ##############################
%---------------------------------------------------------------




for i = 1:7
    FY7450_mean(i) = mean(Data.F7450(i).data.a2(16).data);
    FY7460_mean(i) = mean(Data.F7460(i).data.a2(16).data);
    FY7500_mean(i) = mean(Data.F7500(i).data.a2(16).data);
    FY7510_mean(i) = mean(Data.F7510(i).data.a2(16).data);
    
end






%% Z-direction analysis
for i = 1:7
    FY7450_mean(i) = mean(Data.F7450(i).data.a2(16).data);
    FY7460_mean(i) = mean(Data.F7460(i).data.a2(16).data);
    FY7500_mean(i) = mean(Data.F7500(i).data.a2(16).data);
    FY7510_mean(i) = mean(Data.F7510(i).data.a2(16).data);
    
end














