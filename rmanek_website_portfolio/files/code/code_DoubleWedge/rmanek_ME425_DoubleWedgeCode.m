%% Double Wedge Lab
% Author: Rushan Manek
% BUID: U89490767
% Date Created: 14/OCT/22
% Date of Last Revision: 15/OCT/22

%% Clear Previous Variables
clc
clf
clear
%close all
%% Initial Conditiions
gamma = 1.4;
region(1).Region = 1;

%% Asking User for Input Information
region(1).M = input("Enter the Upstream Mach #: ");
region(1).P = input("Enter the Upstream Static Pressure: ");
AoA = input("Enter the Angle of Attack (degrees): ") * (pi/180);
theta_f = input("Enter the Forward Wedge Half Angle (degrees): ") * (pi/180);
theta_a = input("Etner the Rearward Wedge Half Angle (degrees): ") * (pi/180);



%% Creating the Points for the Double Wedge

c = 1;
x0 = 0;
y0 = 0;
xp = (c * tan(theta_a)) / (tan(theta_f) + tan(theta_a));
yp = xp * tan(theta_f);


%% Creating the Double Wedge
edge1 = [x0 y0; xp yp];
edge2 = [x0 y0; xp -yp];
edge3 = [xp yp; c 0];
edge4 = [xp -yp; c 0];

%% Plotting the Double Wedge
airfoil = figure(1);
hold on
plot(edge1(:,1),edge1(:,2),Color='#15a322',LineWidth=1.5,DisplayName="Region 2");
plot(edge2(:,1),edge2(:,2),'b',LineWidth=1.5,DisplayName="Region 3");
plot(edge3(:,1),edge3(:,2),'k',LineWidth=1.5,DisplayName="Region 4");
plot(edge4(:,1),edge4(:,2),'r',LineWidth=1.5,DisplayName="Region 5");
xlabel("Airfoil Length (m)")
ylabel("Airfoil Height (m)")
title("Airfoil Cross Section")
legend

%Adjust scaling to look like airfoil
xlim([x0-(c*0.3), c+(c*0.3)])
ylim([-yp-(yp*1), yp+(yp*1)])

%Adjust position and size to centre of screen
%airfoil.Position = [100 100 700 400];
%movegui('center');

%% Calculate Angles
theta2 = theta_f - AoA;
theta3 = theta_f + AoA;
theta4 = theta_f + theta_a;
theta5 = theta_f + theta_a;


%% Calculate Mach # & Pressure for Each Region
if theta2 > 0
    region(2) = calcMnMP_oblique(2,region(1).M,theta2,gamma,region(1).P);
    region(3) = calcMnMP_oblique(3,region(1).M,theta3,gamma,region(1).P);
    region(4) = calcPM_expansion(4,region(2).M,region(2).P,theta4,gamma);
    region(5) = calcPM_expansion(5,region(3).M,region(3).P,theta5,gamma);
elseif theta2 < 0
    region(2) = calcPM_expansion(2,region(1).M,region(1).P,theta2,gamma);
    region(3) = calcMnMP_oblique(3,region(1).M,theta3,gamma,region(1).P);
    region(4) = calcPM_expansion(4,region(2).M,region(2).P,theta4,gamma);
    region(5) = calcPM_expansion(5,region(3).M,region(3).P,theta5,gamma);
end

%% Calculate Force for Each Region
region(1).Fx = 0;
region(1).Fy = 0;
for i = 2:5
    [region(i).Fx, region(i).Fy] = calcForce(i,region(i).P,xp,yp,c);
end

%% Calculate Total Force, Lift & Drag
Fx_total = sum([region.Fx]);
Fy_total = sum([region.Fy]);

L = (Fy_total * cos(AoA)) - (Fx_total * sin(AoA));
D = (Fx_total * cos(AoA)) + (Fy_total * sin(AoA));

Cl = (2 * L) / (gamma * region(1).P * (region(1).M)^2 * c);
Cd = (2 * D) / (gamma * region(1).P * (region(1).M)^2 * c);
clc

%% Print Information

%Input Table Information 
fprintf("\nInput Information\n")
inputTable = table(region(1).M,region(1).P,AoA*(180/pi),theta_f*(180/pi),theta_a*(180/pi));
inputTableTitles = ["M1 ","P1 (kPa)",'α (deg)', "α_f (deg)","α_a (deg)"];
inputTable.Properties.VariableNames = inputTableTitles;
disp(inputTable)

disp('_________________________________________________________________')


%Region Information
fprintf("\nRegion Information\n")
regionTable = struct2table(region);
regionTableTitles = ["Region", "Mach #", "Pressure (kPa)", "Fx (kN)", "Fy (kN)"];
regionTable.Properties.VariableNames = regionTableTitles;
disp(regionTable)

disp('_________________________________________________________________')


%Output Information
fprintf("\nOutput Information\n")
outputTable = table(Fx_total*1000, Fy_total*1000, L*1000, D*1000, Cl, Cd);
outputTableTitles = ["Total Fx (N)", "Total Fy (N)", "Lift Force (N)", "Drag Force (N)","Cl","Cd"];
outputTable.Properties.VariableNames = outputTableTitles;
disp(outputTable)
fprintf("\n\n")

%% Plot Annotation

annotationLocation = [x0-(c*0.25), 0;
                      x0-(c*0.15), yp;
                      x0-(c*0.15), -yp;
                      c, yp/2;
                      c, -yp/2];
for i = 1:5
    txt = plotText(i,round(region(i).M,1),round(region(i).P,1));
    text(annotationLocation(i,1),annotationLocation(i,2),txt,EdgeColor='k')
end








