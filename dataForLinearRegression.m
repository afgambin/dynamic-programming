% Script to obtain data for multiple linear regression

% Creo diferentes estructuras y varío el número de estaciones pico activas,
% tráfico, ABS y CRE para obtener datos.
% Los datos obtenidos los guardo directamente en la carpeta de regresion
% lineal como dataConsumption.mat

%debug
clc,clear all, close all
nStructs = 5;
display('Calculating data...')

[stat,path] = fileattrib;
pathCurrent = path.Name;
pathFolder = [pathCurrent '\regresion lineal'];

ISD = 500; % Inter-Site Distance between BSs
r = ISD/3;
apothem = sqrt(3)/2 * r;
nMaxUE =  100; % Maximum number of UE associated to a (macro or pico) BS.
plotting = 0;   % Activation of graphs

t = [0 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 65];

ABSvalues = 1:7;
CREvalues = [0 6 9 12 18];
armsValues = [];

for i=1:length(ABSvalues)
    aux = [repmat(ABSvalues(i),length(CREvalues),1) CREvalues'];
    armsValues = [armsValues; aux];
end

dataConsumptionMacroBlockade = [];
dataConsumptionPico = [];
nPicos = 12; % Number of picos per sector
matrixStates = zeros(nPicos+1, nPicos);

for i=2:size(matrixStates,1)
    aux = ones(1,i-1);
    aux2 = zeros(1,nPicos - (i-1));
    matrixStates (i,:)=[aux2 aux];
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
        X = matrixStates(i,:);
        for j=1:length(t)
            for k=1:length(armsValues)
                nABSsubframes = armsValues(k,1);
                CRE = armsValues(k,2);
                [returnDataMacroBlockade, returnDataPico] = consumptionData(X, t(j), nABSsubframes, CRE);
                dataConsumptionMacroBlockade = [dataConsumptionMacroBlockade; returnDataMacroBlockade];
                dataConsumptionPico = [dataConsumptionPico; returnDataPico];
            end
        end
    end
end

save([pathFolder '\dataConsumption.mat'],'dataConsumptionMacroBlockade','dataConsumptionPico');
display('Stored Data.')
