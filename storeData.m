function [] = storeData(nPicos, type, dailyConsumption, dailyBlockingProb)

pathFolder = 'C:\Universidad\TFM\código simulador\v8\simulaciones\differentTraffic';

if strcmp(type,'NS')
    cNS = dailyConsumption;
    probNS = dailyBlockingProb;
    mat = ['\gNS' num2str(nPicos)];
    save([pathFolder mat],'cNS','probNS');
    
elseif strcmp(type,'Bench')
    cBench = dailyConsumption;
    probBench = dailyBlockingProb;
    mat = ['\gBench' num2str(nPicos)];
    save([pathFolder mat],'cBench','probBench');
    
else
    cDP = dailyConsumption;
    probDP = dailyBlockingProb;
    mat = ['\gDP' num2str(nPicos)];
    save([pathFolder mat],'cDP','probDP');
end
display('Stored Data')
end