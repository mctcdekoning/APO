% Optimization file problem 1.1
clear all
close all
clc

%% Add Cplex paths
warning('off','MATLAB:dispatcher:pathWarning')
% Karel
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64')
% Mathias
addpath ('/Users/mathiasdekoning/Applications/IBM/ILOG/CPLEX_Studio1271/cplex/matlab/x86-64_osx')
% Gael


%% Determine input

% select input file and sheet
data    = [pwd '/AE4423_Datasheet_20.xlsx'];

% number of EU airports
NodesEU   = input('Number of EU airports (between 5 and 20): ');
while ismember(NodesEU,5:20) == 0
    NodesEU = input('Invalid input. Number of EU airports (between 5 and 20): ');
end

% number of US airports
NodesUS   = input('Number of US airports (between 0 and 4): ');
while ismember(NodesUS,0:4) == 0
    NodesEU = input('Invalid input. Number of US airports (between 0 and 4): ');
end

Nodes = NodesEU + NodesUS;

% number of AC types
ACtype = 5;

% number of AC per AC type
disp('Number of AC per type: ')       
AC1  = input('type 1: '); 
AC2  = input('type 2: ');
AC3  = input('type 3: ');
AC4  = input('type 4: ');
AC5  = input('type 5: ');       

% total number of AC
AC = AC1 + AC2 + AC3 + AC4 + AC5;

%% import airport data

AirportcellsEU      = ['C4:', char(66+NodesEU), '5'];
[~,AirportsEU,~]    = xlsread(data,1,AirportcellsEU);
LatcellsEU          = ['C6:', char(66+NodesEU), '6'];
LongcellsEU         = ['C7:', char(66+NodesEU), '7'];
RunwaycellsEU       = ['C8:', char(66+NodesEU), '8'];
DemandcellsEU       = ['C16:', char(66+NodesEU), int2str(15+NodesEU)];
if NodesUS ~= 0
    AirportcellsUS      = ['W4:', char(86+NodesUS), '5'];
    [~,AirportsUS,~]    = xlsread(data,1,AirportcellsUS);
    Airports            = [AirportsEU, AirportsUS];
    LatcellsUS          = ['W6:', char(86+NodesUS), '6'];
    Lat_deg             = [xlsread(data,1,LatcellsEU),xlsread(data,1,LatcellsUS)]; %[deg]
    LongcellsUS         = ['W7:', char(86+NodesUS), '7'];
    Long_deg            = [xlsread(data,1,LongcellsEU),xlsread(data,1,LongcellsUS)];       %[deg]
    RunwaycellsUS       = ['W8:', char(86+NodesUS), '8'];
    RunwayAir           = [xlsread(data,1,RunwaycellsEU),xlsread(data,1,RunwaycellsUS)];
    DemandcellsUSi      = ['C36:', char(66+NodesEU), int2str(35+NodesUS)];
    DemandcellsUSj      = ['W16:', char(86+NodesUS), int2str(15+NodesEU)];
    DemandcellsUSij     = ['W36:', char(86+NodesUS), int2str(35+NodesUS)];
    Demand              = [xlsread(data,1,DemandcellsEU),xlsread(data,1,DemandcellsUSj);
                           xlsread(data,1,DemandcellsUSi),xlsread(data,1,DemandcellsUSij)];
else
    [~,Airports,~]  = AirportsEU;
    Lat_deg         = xlsread(data,1,LatcellsEU);     %[deg]
    Long_deg        = xlsread(data,1,LongcellsEU);    %[deg]
    RunwayAir       = xlsread(data,1,RunwaycellsEU);
    Demand          = xlsread(data,1,DemandcellsEU);
end

Lat       = deg2rad(Lat_deg);              %[rad]
Long      = deg2rad(Long_deg);             %[rad]

disp('Your selected airports are '), disp(Airports)

%% import aircraft data

ac           = 5;
v_ac         = xlsread(data,1,'AC7:AG7');       % [km/hr]
s_ac         = xlsread(data,1,'AC8:AG8');       % [-]
TAT_ac       = xlsread(data,1,'AC9:AG9');       % [min]
r_ac         = xlsread(data,1,'AC10:AG10');     % [km]
Runway_ac    = xlsread(data,1,'AC11:AG11');     % [m]

% Cost Parameters
Cost_Lease_ac  = xlsread(data,1,'AC13:AG13');   % [Euro]
Cost_Fixed_ac  = xlsread(data,1,'AC14:AG14');   % [Euro]
Cost_Time_ac   = xlsread(data,1,'AC14:AG14');   % [Euro/hr]
Cost_Fuel_ac   = xlsread(data,1,'AC15:AG15');   % [Euro/km]

%% Calculate the distance between the nodes:

dd_ij = zeros(Nodes,Nodes);
Re = 6371;  % [km]

for i = 1:Nodes;
    for j = 1:Nodes;
        ddd_ij = 2*asin(sqrt(((sin((Lat(i)-Lat(j))/2))^2)+(cos(Lat(i))*cos(Lat(j))*((sin((Long(i)-Long(j))/2))^2)))); 
        dd_ij (i,j) = Re * ddd_ij;
    end
end     % verified Londen-Paris = 350 [km]

% Convert to correct matrix
d_ij = reshape (dd_ij, Nodes*Nodes, 1);


%% Calculate the yield

YY_ij = zeros(Nodes,Nodes);

for i = 1:NodesEU
    for j = 1:NodesEU
        YY_ij(i,j) = (5.9.*dd_ij(i,j).^(-0.67))+0.043;
    end
    for j = NodesEU+1:Nodes
        YY_ij(i,j) = 0.05*dd_ij(i,j);
    end
end
for i = NodesEU+1:Nodes
    for j = 1:Nodes
        YY_ij(i,j) = 0.05*dd_ij(i,j);
    end
end

Y_ij = reshape(YY_ij,Nodes*Nodes, 1);


%% Initiate CPLEX

% Decision variables
DV  = Nodes * Nodes * ACtype;      % Number of DV = Nodes (i) * Nodes (j) * ACtype (k)

% function out = varindex(m, n, p)
%     out = (m-1) * Nodes + n + Nodes*Nodes*(p-1)
% end


