function [nextX, macroPowerVector, CRE, time] = dynamicProgramming(X, currentTraffic, sec, blockingThreshold, RupperObjective)
tic
%debug
% clc,clear all, close all
% X = [0 1 0 0];
% currentTraffic = 5;
% sec = 1;
% blockingThreshold = 0.2;
% RupperObjective = 250000;

%% OFFLINE Phase

% LR data
data = load('dataLR');

%Data
thetaBlockade = data.lrStruct.blockade.thetaBlockade;
muBlockade = data.lrStruct.blockade.muBlockade;
sigmaBlockade = data.lrStruct.blockade.sigmaBlockade;

thetaPico = data.lrStruct.pico.thetaPico;
muPico = data.lrStruct.pico.muPico;
sigmaPico = data.lrStruct.pico.sigmaPico;

thetaMacro = data.lrStruct.macro.thetaMacro;
muMacro = data.lrStruct.macro.muMacro;
sigmaMacro = data.lrStruct.macro.sigmaMacro;

% Activation Cost
betta = 0.1; %[activation]
picoPowerBase = 56; %[W]
activationCost = picoPowerBase*betta;

stages = 24;
maxUEtoSector = 300;
averageDailyTraffic = [0.125 0.1 0.075 0.05 0.03 0.024 0.024 0.027 0.04 0.07 0.085 0.0964 0.1 0.105 0.11 0.115 0.12 0.1225 0.125 0.13 0.14 0.16 0.155 0.15];
t = maxUEtoSector*averageDailyTraffic;
nPicos = length(X);

load('meanABSCRE');
% Aquí estan guardados los arrays "promedios" y "armsValues". Los utilizo
% para saber los valores óptimos esperados de ABS y CRE en función del
% tráfico y el estado del sistema

% Consumption matrix
if exist('MatrixC.mat', 'file') %If the file exists, it is loaded
    load('MatrixC');
end

%If not, we calculate it
if ~exist('MatrixC.mat', 'file') || (nPicos+1) ~= size(C,1) 
    C = zeros(nPicos+1,length(t));
    display('Calculating Consumption Matrix...')
    for i=1:nPicos+1
        s = ones(1,i-1); 
        s2 = zeros(1,nPicos - (i-1));
        arrayX = [s2 s];       
        for j=1:length(t)
            partial = find(promedios(:,1) == (i-1));
            partial = promedios(partial,:);
            [~, indexArm] = min(abs(partial(:,2)-t(j)));
            arm = armsValues(partial(indexArm,3),:)
            nABSsubframes = arm(1);
            CRE = arm(2);
            
            C(i,j) = calculateConsumption2(arrayX,t(j), thetaBlockade, muBlockade, sigmaBlockade, thetaPico, muPico, sigmaPico, thetaMacro, muMacro, sigmaMacro, blockingThreshold, nABSsubframes, CRE);
        end
    end
    save('MatrixC','C');
end

% Cost Activation Matrix (Don't Delete!!!)
if exist('activationMatrix.mat', 'file') %If the file exists, it is loaded
    load('activationMatrix');
end

%If not, we calculate it
if ~exist('activationMatrix.mat', 'file') || (nPicos+1) ~= size(activationMatrix,1)
    activationMatrix = zeros(nPicos+1,nPicos+1);
    display('Calculating Activation Cost Matrix...')
    for i=1:nPicos+1
        activationMatrix(i,:) = calculateActivationCostMatrix(i-1,nPicos);
    end
    save('activationMatrix','activationMatrix');
end

% Creación de la tabla o matriz de la Programación Dinámica. Para ello se
% hace uso de la matriz de activación y la matriz de consumo antes
% calculadas. 

% Dynamic Programming
if exist('MatrixJ.mat', 'file') %If the file exists, it is loaded
    load('MatrixJ');
end

%If not, we calculate it
if ~exist('MatrixJ.mat', 'file') || (nPicos+1) ~= size(J,1)
    J = zeros(nPicos+1,stages);
    display('Calculating J Matrix...')
    for k = stages:-1:1
        validStates = find(C(:,k)~= -1);
        for s = 1:length(validStates)
            %Activation Cost
            activationArray = activationMatrix(validStates(s),:)*activationCost;
            
            if k == stages
                validActions = find(C(:,1)~= -1);   % It is valid too: validActions = C(1,find(C(1,:)~= -1));
                validValues = C(validActions,1);
                costs = validValues + activationArray(validActions)';
            else
                validActions = find(C(:,k+1)~= -1);
                validValues = C(validActions,k+1);
                costs = validValues + activationArray(validActions)' + J(validActions,k+1); 
            end 
            J(validStates(s),k) = min(costs); 
        end
    end
    save('MatrixJ','J');
end

%% ONLINE Phase

% La fase anterior solo se realizará la primera vez que se ejecute el
% simulador, posteriormente se accederá directamente a esta fase. Se
% determina la acción óptima para la siguiente etapa en funión del estado
% actual y del tráfico actual

% currentTraffic is rounded to the nearest column of J
[~, k] = min(abs(t-currentTraffic));
s = sum(X == 1) + 1; % current State -> row of J

activationArray = activationMatrix(s,:)*activationCost;

if k == stages
    validActions = find(C(:,1)~= -1);   % It is valid too: validActions = C(1,find(C(1,:)~= -1));
    validValues = C(validActions,1);
    costs = validValues + activationArray(validActions)' + J(validActions,1);
else
    validActions = find(C(:,k+1)~= -1);   % It is valid too: validActions = C(1,find(C(1,:)~= -1));
    validValues = C(validActions,k+1);
    costs = validValues + activationArray(validActions)' + J(validActions,k+1); 
end

[~,rowControl] = min(costs); % if there is a tie, it is chosen the control with less active picos
controlActivePicos = validActions(rowControl) - 1;

% Las estaciones pico que se deben apagar o encender se escogen de manera
% aleatoria

% Randomly
currentActivePicos = sum(X == 1);
nextX = X; 

if currentActivePicos < controlActivePicos % Switch on picos is needed
    nOnPicos = controlActivePicos - currentActivePicos;
    auxArray = find(X == 0);

    for i=1:nOnPicos
        picoChosen = randi(length(auxArray));
        nextX(auxArray(picoChosen)) = 1;
        auxArray(picoChosen) = [];
    end
else   % Switch off picos is needed
    nOffPicos = currentActivePicos - controlActivePicos;
    auxArray = find(X == 1);

    for i=1:nOffPicos
        picoChosen = randi(length(auxArray));
        nextX(auxArray(picoChosen)) = 0;
        auxArray(picoChosen) = [];
    end
end

% Turn off nearest pico to macro BS
% nextX = zeros(1,nPicos);
% distancePico2Macro = sqrt(sum((picosPos - repmat(macroPos,nPicos,1)) .^ 2,2));
% sortedDistances = sort(distancePico2Macro,'descend');
% 
% for i=1:controlActivePicos
%     picoChosen = find(distancePico2Macro == sortedDistances(i));
%     nextX(picoChosen) = 1;
% end

%% ABS + CRE OPTIMIZATION Phase

% Finalmente, aplicación del UCB para determinar los valores de ABS y CRE
% óptimos para la siguiente etapa
[optimalABS, optimalCRE] = UCB(nextX,sec, currentTraffic, RupperObjective);
CRE = optimalCRE;
% ABS
nSubframes = 8;
macroPowerVector = [zeros(1,nSubframes-optimalABS) ones(1,optimalABS)];

time = toc;
end