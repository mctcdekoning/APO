function Multicommodity ()
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
data = [pwd '/AE4423_Datasheet_20.xlsx'];

% number of EU airports
NodesEU = input('Number of EU airports (between 5 and 20): ');
while ismember(NodesEU,5:20) == 0
    NodesEU = input('Invalid input. Number of EU airports (between 5 and 20): ');
end

% number of US airports
NodesUS = input('Number of US airports (between 0 and 4): ');
while ismember(NodesUS,0:4) == 0
    NodesEU = input('Invalid input. Number of US airports (between 0 and 4): ');
end

Nodes = NodesEU + NodesUS;

% number of AC types and AC per AC type
disp('Number of AC per type: ')
ACtype = 0;
AC1  = input('type 1: '); 
if AC1 ~= 0
    ACtype = ACtype + 1;
end
AC2  = input('type 2: ');
if AC2 ~= 0
    ACtype = ACtype + 1;
end
AC3  = input('type 3: ');
if AC3 ~= 0
    ACtype = ACtype + 1;
end
AC4  = input('type 4: ');
if AC4 ~= 0
    ACtype = ACtype + 1;
end
AC5  = input('type 5: ');
if AC5 ~= 0
    ACtype = ACtype + 1;
end

% total number of AC
AC_total = AC1 + AC2 + AC3 + AC4 + AC5;
AC = [AC1, AC2, AC3, AC4, AC5];

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
    Airports        = AirportsEU;
    Lat_deg         = xlsread(data,1,LatcellsEU);     %[deg]
    Long_deg        = xlsread(data,1,LongcellsEU);    %[deg]
    RunwayAir       = xlsread(data,1,RunwaycellsEU);
    Demand          = xlsread(data,1,DemandcellsEU);
end

Lat       = deg2rad(Lat_deg);              %[rad]
Long      = deg2rad(Long_deg);             %[rad]

disp('Your selected airports are '), disp(Airports)

% Transfer passengers from departing from the hub are not allowed
Demand_G = Demand;
Demand_G(5,:) = zeros(1,Nodes);
Demand_G(:,5) = zeros(Nodes,1);


%% import aircraft data

actypes      = [AC1~=0,AC2~=0,AC3~=0,AC4~=0,AC5~=0];
v_ac         = xlsread(data,1,'AC7:AG7');       % [km/hr]
v_ac         = v_ac(v_ac.*actypes~=0);          % [km/hr]
s_ac         = xlsread(data,1,'AC8:AG8');       % [-]
s_ac         = s_ac(s_ac.*actypes~=0);          % [-]
TAT_ac       = xlsread(data,1,'AC9:AG9');       % [min]
TAT_ac       = TAT_ac(TAT_ac.*actypes~=0);      % [min]
r_ac         = xlsread(data,1,'AC10:AG10');     % [km]
r_ac         = r_ac(r_ac.*actypes~=0);          % [km]
Runway_ac    = xlsread(data,1,'AC11:AG11');     % [m]
Runway_ac    = Runway_ac(Runway_ac.*actypes~=0);% [m]
LF_EU        = O.75                             % [-]
LF_US        = O.85                             % [-]

% Cost Parameters
Cost_Lease_ac  = xlsread(data,1,'AC13:AG13');   % [Euro]
Cost_Lease_ac  = Cost_Lease_ac(Cost_Lease_ac.*actypes~=0);
Cost_Fixed_ac  = xlsread(data,1,'AC14:AG14');   % [Euro]
Cost_Fixed_ac  = Cost_Fixed_ac(Cost_Fixed_ac.*actypes~=0);
Cost_Time_ac   = xlsread(data,1,'AC14:AG14');   % [Euro/hr]
Cost_Time_ac   = Cost_Time_ac(Cost_Time_ac.*actypes~=0);
Cost_Fuel_ac   = xlsread(data,1,'AC15:AG15');   % [Euro/km]
Cost_Fuel_ac   = Cost_Fuel_ac(Cost_Fuel_ac.*actypes~=0);

%% Calculate the distance between the nodes:

dd_ij = zeros(Nodes,Nodes);
Re = 6371;  % [km]

for i = 1:Nodes
    for j = 1:Nodes
        ddd_ij = 2*asin(sqrt(((sin((Lat(i)-Lat(j))/2))^2)+(cos(Lat(i))*cos(Lat(j))*((sin((Long(i)-Long(j))/2))^2)))); 
        dd_ij (i,j) = Re * ddd_ij;
    end
end     % verified Londen-Paris = 350 [km]

% Convert to correct matrix
d_ij = reshape (dd_ij, Nodes*Nodes, 1);

% Define matrix a_ijk which is 0 is range of aircraft is not sufficient for distance airport i to j
a_ijk=zeros(Nodes,Nodes,ACtype);
for i = 1 : Nodes
    for j = 1:Nodes
        for k = 1: ACtype;
            if dd_ij(i,j)<r_ac(k);
                a_ijk(i,j,k)=1000;
            end
        end
    end
end

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

%% Calculate the LF

LF = zeros(Nodes,Nodes);

for i = 1:NodesEU
    for j = 1:NodesEU
        LF(i,j) = LF_EU;
    end
    for j = NodesEU+1:Nodes
        LF(i,j) = LF_US;
    end
end
for i = NodesEU+1:Nodes
    for j = 1:Nodes
        LF(i,j) = LF_US
    end
end

%% Calculated TAT

TAT_ac

%% Initiate CPLEX

% Decision variables
% Number of DV = Nodes (i) * Nodes (j) * ACtype (k)
DV = (Nodes * Nodes) + (Nodes * Nodes) + ACtype + (Nodes * Nodes) * ACtype;

% Initialize the CPLEX object
model               = 'Cplex_model';
cplex               = Cplex('model');
cplex.Model.sense   = 'maximize';

%% Objective function

% Variable operational costs
for k = 1:ACtype;
    % multiple dimension matrix of fixed cost AC
    C_fix(:,:,k) = ones(Nodes*Nodes,1)*Cost_Fixed_ac(k);
    % multiple dimension matrix of time-based cost AC
    C_time(:,:,k)  =   Cost_Time_ac(k).*(d_ij./v_ac(k));
    % multiple dimension matrix of fuel cost AC
    C_fuel(:,:,k)  =   Cost_Fuel_ac(k).*d_ij;
end

Cost_X_ij   =   Y_ij.*d_ij;
Cost_W_ij   =   Y_ij.*d_ij;
Cost_Z_ij   =   C_fix + C_time + C_fuel;

Cost_X      =   reshape(Cost_X_ij, numel(Cost_X_ij), 1);
Cost_W      =   reshape(Cost_W_ij, numel(Cost_W_ij), 1);
Cost_Z      =   reshape(Cost_Z_ij, numel(Cost_Z_ij), 1);

Cost_AC      =   transpose(Cost_Lease_ac(:,1:ACtype));

        obj     =      [Cost_X ; Cost_W ; -Cost_AC ; -Cost_Z];
        lb      =      zeros(DV,1);
        ub      =      inf(DV,1);
        ctype   =      char(ones(1, (DV)) * ('I')); 
        
% Array with DV names
        l = 1; 
        for i = 1:Nodes;
            for j = 1:Nodes;    % of the X_{ij} variables
                NameDV (l,:)  = ['X_' num2str(i,'%02d') ',' num2str(j,'%02d') '_' num2str(0,'%02d')];
                l = l + 1;
            end;
        end;
        
        for i = 1:Nodes;
            for j = 1:Nodes;    % of the W_{ij} variables
                NameDV (l,:)  = ['W_' num2str(i,'%02d') ',' num2str(j,'%02d') '_' num2str(0,'%02d')] ;
                l = l + 1;
            end;
        end;
        
        for k = 1:ACtype        % of the AC^k variables
            NameDV (l,:)  = ['A_' num2str(0,'%02d') ',' num2str(0,'%02d') '_' num2str(k,'%02d')];
            l = l + 1;           
        end
        
        for k = 1:ACtype        % of the Z_{ij}^k variables
            for i = 1:Nodes
                for j = 1:Nodes
                    NameDV (l,:)  = ['Z_' num2str(i,'%02d') ',' num2str(j,'%02d') '_' num2str(k,'%02d')];
                    l = l + 1;
                end
            end
        end

cplex.addCols(obj, [], lb, ub, ctype, NameDV);

%% Constraints

% All flow each airport leaves the airport:
for i = 1:Nodes
        for j = 1:Nodes
            C1 = zeros(1, DV);                          %Setting coefficient matrix with zeros
            C1(Xindex(i,j)) = 1;                        % Passenger 
            C1(Windex(i,j)) = 1;                        %Link getting OUT the node
            cplex.addRows(-inf, C1, Demand(i,j), sprintf('FlowBalanceNode%d_%d',i,j))
        end
end

%Transfer passengers aircraft only are only if hub is not origin or destination:
for i = 1:Nodes
    for j = 1:Nodes
        C2 = zeros(1, DV);
        C2(Windex(i,j)) = 1;
        cplex.addRows(-inf, C2, Demand_G(i,j), sprintf('TransferPassengerNotHub%d_%d',i,j))
    end
end

%Capacity verification in each flight leg:


%Balance incoming outgoing flight
for i = 1:Nodes;
    for k = 1:ACtype;
        C4 = zeros(1, DV);
        for j = Nodes
            C4(Zindex(i,j,k))    = 1;
            C4(Zindex(j,i,k))    = -1;
        end
        cplex.addRows(0,C4,0,sprintf('FlowBalanceNode%d_%d',i,j,k)) ;
    end
end

%Use aircraft limited to block hours
for k = 1:ACtype;
    C5 = zeros(1,DV);
    for i = 1:Nodes;
        for j = 1:Nodes;
            C5(Zindex(i,j,k)) = (dd_ij(i,j)/v_ac(k))+TAT_ac(i,j,k);
        end
    end
    C5(ACindex(k)) = -(70*AC(k));
    cplex.addRows(-inf,C5,0,sprintf('FlowBalanceNode%d_%d',i,j,k));
end
%       !!!!!TAT_ac(i,j,k) needs to be build!!!!!!!

% number of aircraft type k equals the fleet number:
for k = 1:ACtype
    C6 = zeros(1, DV);
    C6(ACindex(k)) = n(k);
    C6(ACindex(k)) = AC(k);
    cplex.addRows(0,C6,0,sprintf('Fleetnumber%d_%d',k));
end

% Aircraft flying to airport is able to fly the distance, 
for i = 1:Nodes
    for j = 1:Nodes
        for k = 1:ACtype
            C7 = zeros(1,DV);
            C7(Zindex(i,j,k)) = 1          
            cplex.addRows(-inf,C7,a_ijk(i,j,k),sprintf('MaxRange%d_%d',i,j,k));
        end
    end
end


%% Execute model

% cplex.solve();
% cplex.writeModel([model '.lp'])


%% 
    function out = Xindex(m, n)
        out = (m - 1) * Nodes + n;  % Function given the variable index for each X(i,j) [=(m,n)]  
              %column          %row   
    end

    function out = Windex(m, n)
        out = Nodes*Nodes + (m - 1) * Nodes + n;  % Function given the variable index for each W(i,j) [=(m,n)]  
              %X counter     %column          %row   
    end
    
    function out = ACindex(p)
        out = p + 2*Nodes*Nodes;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
              %column %X/Wcounter 
    end
    
    function out = Zindex(m, n, p)
        out = 2*Nodes*Nodes + ACtype + (m - 1) * Nodes + n + Nodes*Nodes*(p-1);  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
              %X/W/L Counter            %column          %row %parallel matrixes (k=1 & k=2)
    end
end

