function [] = storeDataTiming(nPicos, type, computingTime, firstTime)

pathFolder = 'C:\Universidad\TFM\código simulador\resultados\computingTime';

if strcmp(type,'Bench')
    timeBench = computingTime;
    firstBench = firstTime;
    mat = ['\gBench' num2str(nPicos)];
    save([pathFolder mat],'timeBench','firstBench');
    
else
    timeDP = computingTime;
    firstDP = firstTime;
    mat = ['\gDP' num2str(nPicos)];
    save([pathFolder mat],'timeDP','firstDP');
end
display('Stored Timing Data')
end