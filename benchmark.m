function [optimumX, time] = benchmark(t,sec, RupperObjective)
% Cell Zooming Algorithm [12] -> Extracted from Genetic Algorithm Reference
%debug
% clc, clear all, close all
% t = 30;
% sec = 1;

tic
% Function arguments
secStr = ['S' num2str(sec)];
newUsers = round(random('poiss', t)); % Traffic: number of UEs

ISD = 500; % Inter-Site Distance between BSs
r = ISD/3;
apothem = sqrt(3)/2 * r;

nSubframes = 8;
macroPowerVector = [1 1 1 1 0 0 0 0]; % ABS

%Pico Values TR 36.887
macroPower = 10^((43 - 30)/10); % 43 dBm --> 13 dBW --> 20 W
picoPower = 6.3;% 6.3 W

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

% Loading struct
load('schemeData');
nPicos = scheme.nPicos;

% Minimun distance constraints
minDistanceUE2macro = 35;
minDistanceUE2pico = 10;
nInterferingMacros =  4;
pFarUE = .1;
picoCoverageR = 40;

% CRE
associationBias = 18; % CRE Bias = 18 dB

% Capacity parameters [bps]
Robjective = 250000;
%RupperObjective = 250000;
nRBs = 50; % 10MHz / 180 KHz per each RB -> 100 RBs???
bandwithRB = 180000;

    
% State Vectors per sector
X = ones(1,nPicos); % All BS are activated

% Number of active Pico BS
nActivePicos = sum(X == 1);
indexesActivePicos = find(X);
activePicosPos = scheme.(secStr).macro.picosPos(indexesActivePicos,:);

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
        elseif macroRSRP > picoRSRP && macroRSRP < (picoRSRP * 10^(associationBias/10))
            winner = ['P' num2str(bestPicoIndex)];      % UE associated with pico BS in CRE
        else
            winner = ['P' num2str(bestPicoIndex)];      % UE associated with pico BS in central region
        end
    else
        winner = 'macro';
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
    else
        dataCapacity = zeros(nSubframes,1);
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
            
            % Algorithm to check how many RBs are being used
            Rusers = zeros(nActiveUEs,1);
            UEsatisfied = zeros(nActiveUEs,1);
            UEcapacityPerRB = UEcapacity/nRBs;
            follow = 1;
            k = 0;
            
            while follow
                indexes = (Robjective - Rusers).* UEcapacityPerRB;
                [~,index] = max(indexes);
                Rusers(index) = Rusers(index) + UEcapacityPerRB(index);
                k = k + 1;
                UEsatisfied(index) = Rusers(index) >= RupperObjective;
                UEcapacityPerRB(index) = (1-UEsatisfied(index))*UEcapacityPerRB(index);
                
                if (sum(UEsatisfied) == nActiveUEs || k == nRBs)
                    follow = 0;
                end
            end
            
            dataCapacity(subf) = k/nRBs;
        end
    end
    % Set BS capacity and Blocking Prob.
    scheme.(['S' num2str(sec)]).(cell).BScapacity = sum(dataCapacity)/length(dataCapacity);
end

% The algorithm starts here
vectorBScapacity = zeros(nPicos,2);
for p = 1:(nPicos)
    cell = ['P' num2str(p)];
    nActiveUEs = sum(scheme.(['S' num2str(sec)]).(cell).UEpos(:,1) ~= -1);
    if nActiveUEs == 0
        X(p) = 0; % Turn off BS
    else
        BScapacity = scheme.(['S' num2str(sec)]).(cell).BScapacity;
        vectorBScapacity(p,:) = [BScapacity p];
    end
end
vectorBScapacity = vectorBScapacity(vectorBScapacity(:,1) ~= 0,:);
indexesPossiblePico = vectorBScapacity(:,2);
%currentActivePicos = length(indexesPossiblePico);
%indexesCurrentActivePicos = indexesPossiblePico;

while ~isempty(vectorBScapacity)
    [~,minIndex] = min(vectorBScapacity(:,1));
    cell = ['P' num2str(vectorBScapacity(minIndex,2))];
    nActiveUEs = sum(scheme.(['S' num2str(sec)]).(cell).UEpos(:,1) ~= -1);
    auxGains = scheme.(['S' num2str(sec)]).(cell).UEgains;
    %auxGains = [auxGains(:,1:nInterferingMacros) auxGains(:,nInterferingMacros+indexesCurrentActivePicos)];
    indexesPossiblePico(minIndex) = [];
    
    reAssociatedUEs = 0;
    for ue = 1:nActiveUEs
        
        uePos = scheme.(['S' num2str(sec)]).(cell).UEpos(ue,:);
        
        % Obtaining capacities
        capacityVector = zeros(length(indexesPossiblePico)+1,1);
        for bs = 1:(length(indexesPossiblePico)+1)
            if bs > length(indexesPossiblePico)
                signalPower = macroPower*auxGains(ue,scheme.(secStr).macro.autoIndex);
                indexBS = scheme.(secStr).macro.autoIndex;
                macroBScapacity = scheme.(['S' num2str(sec)]).macro.BScapacity;
            else
                signalPower = picoPower*auxGains(ue,nInterferingMacros+indexesPossiblePico(bs));
                indexBS = nInterferingMacros+indexesPossiblePico(bs);
            end
            
            auxGains2 = auxGains;
            auxGains2(ue,indexBS) = 0;
            noisePowerWorstCase = [repmat(macroPower,1,nInterferingMacros) repmat(picoPower,1,nActivePicos)] .* auxGains2(ue,:);
            noisePowerBestCase = [zeros(1,nInterferingMacros) repmat(picoPower,1,nActivePicos)] .* auxGains2(ue,:);
            SINRWorst = signalPower / (n0 + sum(noisePowerWorstCase,2)); % We calculate the SINR of each UE
            SINRBest = signalPower / (n0 + sum(noisePowerBestCase,2));
            capacityWorst = W * log2(1 + SINRWorst);
            capacityBest = W * log2(1 + SINRBest);
            capacityVector(bs) = (capacityBest + capacityWorst)/2;
        end
        
        w = capacityVector/W; % Spectral efficiency
        b = Robjective./w; % UE bandwith. Divided by 180 KHz -> necessary nRBs
        auxVectorBScapacity = vectorBScapacity(:,1);
        auxVectorBScapacity(minIndex) = []; %delete actual pico BS
        auxVectorBScapacity = [auxVectorBScapacity; macroBScapacity]; %% add macro BS capacity
        condition = (W*auxVectorBScapacity + b)< W;
        
        if ~(sum(condition) == 0) % UE can be re-associated
            [~,selectedBS] = max(condition.*w); %puede que no sea ninguna
            
            if selectedBS > length(indexesPossiblePico)
                winner = 'macro';
            else
                winner = ['P' num2str(indexesPossiblePico(selectedBS))];
            end
            
            index = find(scheme.(secStr).(winner).UEpos(:,1) == -1,1);
            if isempty(index)
                error('The struct is full of UEs. Increase the ''nMaxUE'' variable.');
            end
            
            % Once the station is selected, the vectors are stored in the appropriate part of the struct.
            scheme.(secStr).(winner).UEpos(index, :) = uePos;
            scheme.(secStr).(winner).UEgains(index,:) = scheme.(['S' num2str(sec)]).(cell).UEgains(ue,:);
            % Update BS capacity
            newBScapacity = scheme.(secStr).(winner).BScapacity + b(selectedBS)/(bandwithRB*nRBs);
            scheme.(secStr).(winner).BScapacity = newBScapacity;
            scheme.(secStr).(cell).UEpos(ue, :) = -1;
            reAssociatedUEs = reAssociatedUEs + 1;
        end
    end
    % Update state vector
    if reAssociatedUEs == nActiveUEs
        X(vectorBScapacity(minIndex,2)) = 0;
        scheme.(secStr).(cell).BScapacity = 0;
    end
    % Delete BS from the list
    vectorBScapacity(minIndex,:) = [];
end
optimumX = X;
time = toc;
end
