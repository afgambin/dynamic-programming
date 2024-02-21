function [] = storeDataFrentePareto(dailyConsumption, dailyBlockingProb, blockingThreshold)

pathFolder = 'C:\Universidad\TFM\código simulador\resultados\frentePareto';

% name = ['cDP' num2str(blockingThreshold*100)];
% nameBlock = ['probDP' num2str(blockingThreshold*100)];
% 
% pareto.consumption.(name) = dailyConsumption;
% pareto.block.(nameBlock) = dailyBlockingProb;

cDP = dailyConsumption;
probDP = dailyBlockingProb;

mat = ['\pareto' num2str(blockingThreshold*1000)];
save([pathFolder mat],'cDP','probDP');

display('Stored Pareto Data')
end