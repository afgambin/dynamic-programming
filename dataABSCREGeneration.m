
% Este script lo utilicé para determinar los valores óptimos esperados de ABS y CRE en cada momento del día.
% Es decir, en función del tráfico y el número de picos activas que
% configuración era la mejor

%debug
clc,clear all, close all
nStructs = 10;
display('Calculating data...')

ISD = 500; % Inter-Site Distance between BSs
r = ISD/3;
apothem = sqrt(3)/2 * r;
nMaxUE =  100; % Maximum number of UE associated to a (macro or pico) BS.
plotting = 0;   % Activation of graphs
sec = 1;
RupperObjective = 250000;

% Si hay tráfico nulo siempre se obtendrá un reward de 1, el UCB debería
% escoger la situación donde la macro menos radia por ejemplo
t = [0 5 10 15 20 25 30 35 40 45 50 55];

nPicos = 10; % Number of picos per sector
matrixStates = zeros(nPicos+1, nPicos);

for i=2:size(matrixStates,1)
    aux = ones(1,i-1);
    aux2 = zeros(1,nPicos - (i-1));
    matrixStates (i,:)=[aux2 aux];
end

% Stats
data = [];
dataArms = [];

ABSvalues = 1:7;
CREvalues = [0 6 9 12 18];
armsValues = [];

for i=1:length(ABSvalues)
    aux = [repmat(ABSvalues(i),length(CREvalues),1) CREvalues'];
    armsValues = [armsValues; aux];
end

for index=1:nStructs
    
    display(index);
    
    [scheme, macroPositions, simulatedSectorCenters] = inicialization(ISD, nPicos, nMaxUE, plotting);
    [scheme] = picoCellGeneration(apothem, nPicos, simulatedSectorCenters, macroPositions, scheme, plotting);
    scheme.nPicos = nPicos;
    scheme.simulatedSectorCenters = simulatedSectorCenters;
    scheme.macroPositions = macroPositions;
    save('schemeData','scheme','simulatedSectorCenters','macroPositions');
 
    for i=1:size(matrixStates,1)
        for j=1:length(t)
            X = matrixStates(i,:);
            arrayArms = zeros(length(armsValues),1);
            for k=1:length(armsValues)
                nABSsubframes = armsValues(k,1);
                CRE = armsValues(k,2);
                reward = calculateReward(X, sec, nABSsubframes, CRE, t(j), RupperObjective);
                arrayAux = [sum(X==1) t(j) nABSsubframes CRE reward];
                data = [data; arrayAux];
                arrayArms(k) = reward;
            end
            [~,bestArm] = max(arrayArms);
            arrayAux2 = [sum(X==1) t(j) bestArm];
            dataArms = [dataArms; arrayAux2];
        end
    end
end

save('dataABSCRE','dataArms','data');

%% Ordenamieno de los datos y visualización
clc, clear all, close all;
load('dataABSCRE');
datos = dataArms;

t = [0 5 10 15 20 25 30 35 40 45 50 55];
nPicos = 10;
promedios = zeros(length(t)*(nPicos+1),3);
cont = 1;

% Stats
data = [];
dataArms = [];

ABSvalues = 1:7;
CREvalues = [0 6 9 12 18];
armsValues = [];

for i=1:length(ABSvalues)
    aux = [repmat(ABSvalues(i),length(CREvalues),1) CREvalues'];
    armsValues = [armsValues; aux];
end

for i=1:(nPicos + 1)
    for j=1:length(t)
        indexes = datos(:,1) == (i-1) & datos(:,2) == t(j);
        valores = datos(indexes,3);
        bestArm = mode(valores);
        promedios(cont,:) = [i-1 t(j) bestArm]; 
        cont =  cont + 1;
    end
end

promedios
save('meanABSCRE','promedios','armsValues');

