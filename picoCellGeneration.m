function [s] = picoCellGeneration(apothem, nPicos, simulatedSectorCenters, macroPositions, s, plotting)

pico2macroMinDistance = 75;
pico2picoMinDistance = 40;

picoPosition = zeros(nPicos * size(simulatedSectorCenters,1),2);

for i = 1:size(simulatedSectorCenters,1) % for each sector we generate nPicos picocells
    
    for j = 1:nPicos
        while 1 % this loop assures the minimum distance constraint
            distance = rand * apothem;
            angle = rand * 2*pi;
            picoPosAux = simulatedSectorCenters(i,:) + ...
                [distance*cos(angle) distance*sin(angle)];
            pico2picoDistance = sqrt(sum((repmat(picoPosAux, nPicos,1) - s.(['S' num2str(i)]).macro.picosPos) .^ 2,2));
            pico2macroDistance = sqrt(sum((repmat(picoPosAux, size(macroPositions,1),1) - macroPositions) .^ 2,2));
            
            if min(pico2macroDistance)> pico2macroMinDistance && min(pico2picoDistance) > pico2picoMinDistance
                break;
            end
        end

        picoPosition((i-1)*nPicos + j,:) = picoPosAux;
        s.(['S' num2str(i)]).macro.picosPos(j,:) = picoPosAux;
        s.(['S' num2str(i)]).(['P' num2str(j)]).picoPos = picoPosAux;
        s.(['S' num2str(i)]).(['P' num2str(j)]).autoIndex = j;
    end
    
end

if plotting
    plot(picoPosition(:,1), picoPosition(:,2),'MarkerFaceColor','g','LineStyle','o', 'MarkerSize', 6)
end
    

end