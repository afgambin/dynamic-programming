% INITIALIZATION PARAMETERS
clc,clear all, close all;

% Cada time slot de duración 1 min se simula una trama LTE.
% El bucle while externo controla el intervalo de confianza de los
% resultados.
% Cada 60 time slots se aplica el algoritmo propuesto (DP), el benchmark (Bench) o no se
% hace nada (NS) en función de la variable "type"
% El esquema detallado del simulador lo puedes ver en el pdf del TFM mejor
% explicado

nTimeSlots = 1440; % Every time slot has a duration of 1 min
nPicos = 4; % Number of picos per sector
%type = 'NS';
%type = 'Bench';
type = 'DP';
store = 0;
storePareto = 0;
storeTiming = 0;
blockingThreshold = 0.2;

relativeTol = 0.05;
objectiveQuality = 0.90;
followWhile = 1;

nSubframes = 8;
macroPowerVector = [1 1 1 1 0 0 0 0]; % ABS
associationBias = 18; % CRE Bias = 18 dB

ISD = 500; % Inter-Site Distance between BSs
nMaxUE =  100; % Maximum number of UE associated to a (macro or pico) BS.

r = ISD/3;
apothem = sqrt(3)/2 * r;
plotting = 1;   % Activation of graphs

%Pico Values TR 36.887
macroPower = 10^((43 - 30)/10); % 43 dBm --> 13 dBW --> 20 W
picoPower = 6.3;% 6.3 W
macroPowerSleep = 75;
picoPowerSleep = 39;
macroPowerBase = 130;
picoPowerBase = 56;
macroSlope = 4.7;
picoSlope = 2.6;
betta = 0.1; %[activation]

% BS Gains
antennaMacroGain = 10^(14/10); % 14 dBi
antennaPicoGain = 10^(5/10); % 5 dBi

%Pico Values TR 36.814
% macroPower = (10^((46 - 30)/10)); % 46 dBm --> 16 dBW --> 39.81 W
% picoPower = (10^((30 - 30)/10)); % 30 dBm --> 0 dBW --> 1W

W = 10e6; % Channel bandwidth (10 MHz) -> Same bandwidth for both types of BSs

F = 10^0.5; % Noise figure = 5 dB
T_O = 290;
K_B = 1.3806504e-23;
BW = 10e6;
N_O = F*T_O*K_B;
n0 = N_O*BW; % White noise at the receiver -> Negligible compared to interference

% Shadowing variance
shadowingVar = 10; % dB

% Initialization of macro and pico positions and simulator struct
if exist('schemeData.mat', 'file') %If the file exists, it is loaded
    display('Loading simulation scenario');
    load('schemeData');
end

%If not, we execute the initialization and save it in a file
if ~exist('schemeData.mat', 'file') || scheme.nPicos ~= nPicos
    display('Creating new simulation scenario');
    [scheme, macroPositions, simulatedSectorCenters] = inicialization(ISD, nPicos, nMaxUE, plotting);
    [scheme] = picoCellGeneration(apothem, nPicos, simulatedSectorCenters, macroPositions, scheme, plotting);
    scheme.nPicos = nPicos;
    scheme.simulatedSectorCenters = simulatedSectorCenters;
    scheme.macroPositions = macroPositions;
    save('schemeData','scheme','simulatedSectorCenters','macroPositions');
end

% Minimun distance constraints
minDistanceUE2macro = 35;
minDistanceUE2pico = 10;
nInterferingMacros =  4;
pFarUE = .1;

nSectors = size(simulatedSectorCenters,1);
picoCoverageR = 40;

% Average daily data traffic profile taken as a reference for a European
% country. This vector represents the percentage of active subscribers per
% hour
averageDailyTraffic = [0.125 0.1 0.075 0.05 0.03 0.024 0.024 0.027 0.04 0.07 0.085 0.0964 0.1 0.105 0.11 0.115 0.12 0.1225 0.125 0.13 0.14 0.16 0.155 0.15];
averagePoly = polyfit(0:23,averageDailyTraffic,7); % Degree polynomial = 7
xAxis = linspace(0,23,1440);
averageValues = polyval(averagePoly,xAxis); % Regression to obtain more points (24 hours * 60 min every hour).
figure, plot(xAxis,averageValues, 'LineWidth',2)
title('Average daily traffic')
title('Tráfico Diario Medio')
xlabel('Tiempo (horas)')
ylabel('% de Usuarios Activos')
grid on
axis tight;
% 
% % Alternative daily data Traffic
% averageDailyTraffic = ones(1,24);
% for i=1:length(averageDailyTraffic)
%     averageDailyTraffic(i) = rand*0.16;
% end
% %load('averageDailyTraffic')
% averagePoly = polyfit(0:23,averageDailyTraffic,7); % Degree polynomial = 7
% xAxis = linspace(0,23,1440);
% averageValues2 = polyval(averagePoly,xAxis); % Regression to obtain more points (24 hours * 60 min every hour).
% figure, plot(xAxis,averageValues2,'r','LineWidth',2)
% hold on
% plot(xAxis, averageValues, 'LineWidth',2)
% title('Patrones de tráfico')
% legend('Patrón inusual', 'Patrón típico')
% xlabel('Tiempo (horas)')
% ylabel('% de Usuarios Activos')
% grid on
% axis tight;

% x = [0:0.01667:24];
% averageValues = 0.12*sin(4*x)+0.12;
% figure, plot(x,averageValues)
% title('Average daily traffic')
% xlabel('Time (hour)')
% ylabel('% Active USers')
% axis tight

counter = 1;
newUsers = 0;
maxUEtoSector = 300;

% Capacity parameters [bps]
Robjective = 250000;
RupperObjective = 250000;
Rmin = 1384000*2;
Rmax = 36696000*2;
nRBs = 100; % 10MHz / 180 KHz per each RB -> 100 RBs???

% Statistics
meanConsumptionMacro = 0;
meanConsumptionPico = 0;
meanConsumptionPerSector = 0;

dailyConsumptionAcum = zeros(nTimeSlots,1);
dailyConsumptionSquareAcum = zeros(nTimeSlots,1);
dailyBlockingProbAcum = zeros(nTimeSlots,1);
dailySamples = 0;

dailyTraffic = zeros(nTimeSlots,1);

nMacroUsersAcum = 0;
nMacroUsersSamples = 0;
nPicoUsersAcum = 0;
nPicoUsersSamples = 0;

computingTimeAcum = zeros(24,1);
computingTimeSamples = 0;
firstTime = 0;
counterTime = 1;

blockAcum = 0;
nBlockSamples = 0;
contadorDay = 1;

%debug
cont2 = 0;
capacityAcum = 0;
capacitySamples = 0;
trafficAcum = zeros(nTimeSlots,1);
trafficSamples = 0;

while followWhile
    fprintf('Simulación día: %f \n', contadorDay);
    contadorDay = contadorDay + 1;
    for slot = 1:nTimeSlots
        
        % Updating  average daily traffic per minute
        counter = mod(slot,length(averageValues));
        if counter == 0
            counter = length(averageValues);
        end
        
        % Traffic Intensity
        rho = maxUEtoSector * averageValues(counter);
        counter = counter + 1;
        
        for sec = 1:1 %DEBUG!!!!!nSectors % For each sector
            secStr = ['S' num2str(sec)];
            
            % State Vectors per sector
            X = scheme.(secStr).X;
            X_previous = scheme.(secStr).X_previous;
            
            % Number of active Pico BS
            nActivePicos = sum(X == 1);
            indexesActivePicos = find(X);
            activePicosPos = scheme.(secStr).macro.picosPos(indexesActivePicos,:);
            
            %debug!!!!
            %         X = [1 0 1 0];
            %         X_previous = [0 0 0 0];
            %         nActivePicos = sum(X == 1);
            %         indexesActivePicos = find(X);
            %         activePicosPos = scheme.(secStr).macro.picosPos(indexesActivePicos,:);
            
            newUsers = round(random('poiss', rho));
            
            dailyTraffic(slot) = dailyTraffic(slot) + newUsers;
            trafficAcum(slot) =  newUsers;
            
            % Random user arrival
            newUEpos = zeros(newUsers,2);
            newUEgains = zeros(newUsers,nActivePicos+4); % +4 Interfering Macros
            
            for newUEindex = 1:abs(newUsers)
                p = randi(nActivePicos + 1); % Index of BS
                if p > nActivePicos % Macro Station UE Generation
                    
                    while 1 % this loop assures the minimum distance constraint
                        distance = rand * apothem;
                        angle = rand * 2*pi;
                        UEPosAux = scheme.(secStr).macro.sectCenter + ...
                            [distance*cos(angle) distance*sin(angle)];
                        
                        UE2macroDistance = sqrt(sum((UEPosAux - scheme.(secStr).macro.macroPos) .^ 2,2));
                        UE2picoDistance = sqrt(sum((repmat(UEPosAux, nPicos,1) - scheme.(secStr).macro.picosPos) .^ 2,2));
                        
                        if UE2macroDistance > minDistanceUE2macro && min(UE2picoDistance) > minDistanceUE2pico
                            break % If minimum distance constraints hold, the loop is breaked
                        end
                        
                    end
                    newUEpos(newUEindex,:) = UEPosAux;
                    
                else % Pico UE generation
                    
                    while 1 % this loop assures the minimum distance constraint
                        if rand < pFarUE % with probability pFarUE a new UE is generated farer than the recomended distance
                            distance = 50;
                        else
                            distance = rand * (picoCoverageR - minDistanceUE2pico) + minDistanceUE2pico; % 10 < distance < 40
                        end
                        
                        angle = rand * 2*pi;
                        UEPosAux = scheme.(secStr).macro.picosPos(p,:) + ...
                            [distance*cos(angle) distance*sin(angle)];
                        
                        UE2macroDistance = sqrt(sum((UEPosAux - scheme.(secStr).macro.macroPos) .^ 2,2));
                        UE2picoDistance = sqrt(sum((repmat(UEPosAux, nPicos,1) - scheme.(secStr).macro.picosPos) .^ 2,2));
                        
                        if UE2macroDistance > minDistanceUE2macro && min(UE2picoDistance) > minDistanceUE2pico
                            break % If minimum distance constraints hold, the loop is breaked
                        end
                    end
                    newUEpos(newUEindex,:) = UEPosAux;
                end
                
                % Macros gain
                InterfMacroIndex = scheme.(secStr).macro.InterfMacroIndex;
                freeSpaceGain = channelModel(repmat(UEPosAux,nInterferingMacros,1), macroPositions(InterfMacroIndex,:), 'm2ue');
                antennaGains =  antennaGain(InterfMacroIndex, macroPositions, UEPosAux);
                newUEgains(newUEindex,1:nInterferingMacros) = (freeSpaceGain .* antennaGains)'*antennaMacroGain;
                
                % Picos gain
                if nActivePicos > 0
                    newUEgains(newUEindex,(nInterferingMacros+1):end) = channelModel(repmat(UEPosAux,nActivePicos,1), activePicosPos, 'p2ue')*antennaPicoGain;
                end
                
                % Adding shadow fading to calculated gain
                newUEgains(newUEindex,:) = newUEgains(newUEindex,:).* 10.^(shadowingVar * randn(1,nInterferingMacros + nActivePicos) / 10);
                
                % UE Association
                % In Macro case, UE can only be associated with the macro
                % of his sector
                bestMacroGain = newUEgains(newUEindex, scheme.(secStr).macro.autoIndex);
                macroRSRP = bestMacroGain * macroPower; % RSRP:Reference Signal Receive Power (LTE)
                
                if nActivePicos > 0
                    [bestPicoGain, picoIndex] = max(newUEgains(newUEindex,(nInterferingMacros+1):end));
                    bestPicoIndex = indexesActivePicos(picoIndex);
                    picoRSRP = bestPicoGain * picoPower;
                    
                    % The UE is associated to the station with the best
                    % RSRP (CRE model).
                    if macroRSRP > picoRSRP && macroRSRP > (picoRSRP * 10^(associationBias/10))
                        winner = 'macro'; % UE associated with macro BS
                        color = 'b';
                    elseif macroRSRP > picoRSRP && macroRSRP < (picoRSRP * 10^(associationBias/10))
                        winner = ['P' num2str(bestPicoIndex)];      % UE associated with pico BS in CRE
                        color = 'r';
                    else
                        winner = ['P' num2str(bestPicoIndex)];      % UE associated with pico BS in central region
                        color = 'r';
                    end
                else
                    winner = 'macro';
                    color = 'b';
                end
                
                index = find(scheme.(secStr).(winner).UEpos(:,1) == -1,1);
                if isempty(index)
                    error('The struct is full of UEs. Increase the ''nMaxUE'' variable.');
                end
                
                % Once the station is selected, the vectors are stored in the appropriate part of the struct.
                scheme.(secStr).(winner).UEpos(index, :) = newUEpos(newUEindex,:);
                scheme.(secStr).(winner).UEgains(index,:) = newUEgains(newUEindex,:);
                %plot(newUEpos(newUEindex,1), newUEpos(newUEindex,2),'MarkerFaceColor',color,'LineStyle','o', 'MarkerSize', 3)
            end
            
            % Scheduling loop: BS Capacity applying ABS
            for p = 1:(nActivePicos+1) % For each active BS
                
                if p > nActivePicos
                    dataCapacity = zeros(sum(macroPowerVector),1);
                    nUEsBlocking = zeros(sum(macroPowerVector),1);
                else
                    dataCapacity = zeros(nSubframes,1);
                    nUEsBlocking = zeros(nSubframes,1);
                end
                
                for subf = 1:nSubframes % For each subframe
                    
                    SubframeMacroPower = macroPower * macroPowerVector(subf);   % ABS
                    
                    if p > nActivePicos % Macro station case
                        cell = 'macro';
                        cellTxPw = SubframeMacroPower;
                        indexBias = scheme.(['S' num2str(sec)]).(cell).autoIndex;
                        boolean = ~(isempty(find(scheme.(['S' num2str(sec)]).(cell).UEpos(:,1) ~= -1,1))) && SubframeMacroPower ~= 0;
                    else % Pico station cases
                        cell = ['P' num2str(indexesActivePicos(p))];
                        cellTxPw = picoPower;
                        indexBias = nInterferingMacros + p;
                        boolean = ~(isempty(find(scheme.(['S' num2str(sec)]).(cell).UEpos(:,1) ~= -1,1)));
                    end
                    
                    if boolean % If there is at least one UE and the macro power in macro scheduling is greater than 0
                        auxGains = scheme.(['S' num2str(sec)]).(cell).UEgains;
                        signalGain = auxGains(:, indexBias); % Vector of gains of usable signal
                        auxGains(:, indexBias) = 0;
                        nActiveUEs = size(auxGains,1);
                        noisePower = [repmat(SubframeMacroPower,nActiveUEs,nInterferingMacros) repmat(picoPower,nActiveUEs,nActivePicos)] .* auxGains;
                        SINR_vector = cellTxPw * signalGain ./ ( n0 + sum(noisePower,2) ); % We calculate the SINR vector of all associated UEs
                        UEcapacity = W * log2(1 + SINR_vector); % Capacity Vector
                        
                        %debug
                        capacityAcum = capacityAcum + mean(UEcapacity);
                        capacitySamples = capacitySamples + 1;
                        
                        % Algorithm to check how many RBs are being used
                        Rusers = zeros(nActiveUEs,1);
                        UEsatisfied = zeros(nActiveUEs,1);
                        UEcapacityPerRB = UEcapacity/nRBs;
                        follow = 1;
                        k = 0;
                        prueba = zeros(nActiveUEs,1);
                        while follow
                            indexes = (Robjective - Rusers).* UEcapacityPerRB;
                            [~,index] = max(indexes);
                            Rusers(index) = Rusers(index) + UEcapacityPerRB(index);
                            prueba(index) = prueba(index) + 1;
                            k = k + 1;
                            UEsatisfied(index) = Rusers(index) >= RupperObjective;
                            UEcapacityPerRB(index) = (1-UEsatisfied(index))*UEcapacityPerRB(index);
                            
                            if (sum(UEsatisfied) == nActiveUEs || k == nRBs)
                                follow = 0;
                            end
                        end
                        
                        nUEsBlocking(subf) = sum(UEsatisfied ~= 1);
                        dataCapacity(subf) = k/nRBs;
                    end
                end
                % Set BS capacity and Blocking Prob.
                scheme.(['S' num2str(sec)]).(cell).BScapacity = sum(dataCapacity)/length(dataCapacity);
                scheme.(['S' num2str(sec)]).(cell).nUEsBlocking = sum(nUEsBlocking)/length(nUEsBlocking);
            end
            
            % Calculating power Consumption per BS and sector Blocking Probability
            consumptionPerSector = 0;
            nUEsBlockingSum = 0;
            for p = 1:(nPicos+1);
                if p > nPicos % Macro station case
                    cell = 'macro';
                    BScapacity = scheme.(['S' num2str(sec)]).(cell).BScapacity;
                    consumptionBS = macroPowerBase + macroSlope*BScapacity*macroPower;
                    meanConsumptionMacro = meanConsumptionMacro + consumptionBS; % stats
                    nUEsBlockingSum = nUEsBlockingSum + scheme.(['S' num2str(sec)]).(cell).nUEsBlocking;
                else % Pico station cases
                    cell = ['P' num2str(p)];
                    BScapacity = scheme.(['S' num2str(sec)]).(cell).BScapacity;
                    
                    if X(p) == 1 && X_previous(p)== 0
                        slopeActivation = betta*picoPowerBase;
                    else
                        slopeActivation = 0;
                    end
                    
                    consumptionBS = X(p)*(picoPowerBase + picoSlope*BScapacity*picoPower) + (1-X(p))*picoPowerSleep + slopeActivation;
                    meanConsumptionPico = meanConsumptionPico + consumptionBS; % stats
                    nUEsBlockingSum = nUEsBlockingSum + scheme.(['S' num2str(sec)]).(cell).nUEsBlocking;
                end
                consumptionPerSector = consumptionPerSector + consumptionBS;
                
                nActiveUEs = sum(scheme.(['S' num2str(sec)]).(cell).UEpos(:,1) ~= -1);
                %Cleaning active UEs
                if nActiveUEs > 0
                    scheme.(['S' num2str(sec)]).(cell).UEpos(:, :) = -1;
                    scheme.(['S' num2str(sec)]).(cell) = rmfield(scheme.(['S' num2str(sec)]).(cell),'UEgains');
                    scheme.(['S' num2str(sec)]).(cell).BScapacity = 0;
                    scheme.(['S' num2str(sec)]).(cell).nUEsBlocking = 0;
                end
            end
            
            % Sector Blocking Probability at this time (current traffic)
            if newUsers == 0
                nBlockSamples = nBlockSamples + 1;
            else
                blockAcum = blockAcum + (nUEsBlockingSum/newUsers);
                nBlockSamples = nBlockSamples + 1;
                dailyBlockingProbAcum(slot) = dailyBlockingProbAcum(slot) + (nUEsBlockingSum/newUsers);
            end
            
            % Activation cost must be added only once
            if X_previous ~= -1
                scheme.(secStr).X_previous = (-1)*ones(1,nPicos);
            end
             
            % DP
            if (mod(slot,60) == 0 && ~strcmp(type,'NS'))   % The process on / off takes place every hour
                if strcmp(type,'DP')
                    display('Applying Dynamic Programming...')
                    [nextX,  newMacroPowerVector, CRE, time] = dynamicProgramming(X,newUsers, sec, blockingThreshold, RupperObjective)
                    macroPowerVector = newMacroPowerVector;
                    associationBias = CRE;
                elseif strcmp(type,'Bench')
                    display('Applying Benchmark Algorithm...')
                    [nextX, time] = benchmark(newUsers, sec, RupperObjective)
                end
                if computingTimeAcum == 0
                    firstTime = time;
                    cont2 = cont2 + 1;
                end
                computingTimeAcum(counterTime) = computingTimeAcum(counterTime) + time;
                counterTime =  counterTime + 1;
                if counterTime == 25
                    counterTime = 1;
                end
                
                scheme.(secStr).X = nextX;
                scheme.(secStr).X_previous = X;
            end
        end
        dailyConsumptionAcum(slot) = dailyConsumptionAcum(slot) + consumptionPerSector;
        dailyConsumptionSquareAcum(slot) = dailyConsumptionSquareAcum(slot) + consumptionPerSector^2;
        meanConsumptionPerSector = meanConsumptionPerSector + consumptionPerSector;
    end
    
    dailySamples = dailySamples + 1;
    computingTimeSamples = computingTimeSamples + 1;
    % Checking if it meets the confidence interval
    if dailySamples > 0
        followWhile = 0;
%         [oneMinusAlpha] = quality(relativeTol, dailySamples, dailyConsumptionAcum, dailyConsumptionSquareAcum);
%         if oneMinusAlpha >= objectiveQuality
%             followWhile = 0;
%         end
    end
    
end

dailyConsumption = dailyConsumptionAcum/dailySamples;
dailyBlockingProb = dailyBlockingProbAcum/dailySamples;
display(mean(dailyConsumption))
display(mean(dailyBlockingProb))

computingTime = computingTimeAcum/computingTimeSamples;

% Store Data
if store == 1
    storeData(nPicos, type, dailyConsumption, dailyBlockingProb);
end

if storePareto == 1
    storeDataFrentePareto(dailyConsumption, dailyBlockingProb, blockingThreshold);
end

if storeTiming == 1
    storeDataTiming(nPicos, type, computingTime, firstTime);
end
