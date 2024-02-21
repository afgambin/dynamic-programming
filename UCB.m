function [optimalABS, optimalCRE] = UCB(X,sec, rho, RupperObjective)
% Código del UCB que me pasaste, lo único que cambia es la función de
% cálculo de la recompensa

%debug
% sec = 1;
% RupperObjective = 250000;
% rho = 35;

nIterations = 200;
ABSvalues = 1:7;
CREvalues = [0 6 9 12 18];
armsValues = [];

for i=1:length(ABSvalues)
    aux = [repmat(ABSvalues(i),length(CREvalues),1) CREvalues'];
    armsValues = [armsValues; aux];
end

nArms = length(armsValues);

counts = zeros(nArms,1);
values = zeros(nArms,1);
ucb_values = zeros(1, nArms);

selectedArmSlot = zeros(1, nIterations);
RewardSlot = zeros(1, nIterations);
valuesSlot = zeros(nArms, nIterations);
UCBvaluesSlot = zeros(nArms, nIterations);

for i = 1:nIterations
    %fprintf('Iteration %d\n', i);
    %select Arm 
    selectedArm = 0;
    for arm = 1:nArms
        if counts(arm) == 0
            selectedArm = arm;
            break
        end
    end
    if selectedArm == 0
        ucb_values = zeros(1, nArms);
        total_counts = sum(counts);
        for arm = 1:nArms
            bonus = sqrt((2*log(total_counts))/(counts(arm)));
            ucb_values(arm) = values(arm) +  bonus;
        end 
        [~, selectedArm] =  max(ucb_values);
    end    
    %fprintf('selectedArm =  %d\n', selectedArm);
    
    %Update
    nABSsubframes = armsValues(selectedArm,1);
    CRE = armsValues(selectedArm,2);
    reward = calculateReward(X,sec, nABSsubframes, CRE, rho, RupperObjective);
   
    counts(selectedArm) = counts(selectedArm) + 1;
    n = counts(selectedArm);
    value = values(selectedArm);
    values(selectedArm) = ((n-1)/n)*value + (1/n)*reward;

    selectedArmSlot(i) = selectedArm;
    RewardSlot(i) =  reward;
    valuesSlot(:,i) = values;
    UCBvaluesSlot(:,i) = ucb_values;
    
%     data.selectedArmSlot = selectedArmSlot;
%     data.RewardSlot = RewardSlot;
%     data.valuesSlot = valuesSlot;
%     data.UCBvaluesSlot = UCBvaluesSlot;
%     
%     save('UCB_data', 'data');

end

% Determine Optimal Arm
optimalArm = mode(selectedArmSlot);
optimalABS = armsValues(optimalArm,1);
optimalCRE = armsValues(optimalArm,2);

end