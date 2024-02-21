function [s, macroPositions, simulatedSectorCenters] = inicialization(ISD, nPicos, nMaxUE, plotting)

%TODO: irregular cell shape
VISD = 2 * ISD * sind(60); % Vertical Inter-Site Distance
nMacro = 19;    % # of macrocells
macroPositions = zeros(nMacro, 2);  % geographical positions of all macros

%simulatedMacros = [5 7 8 10 12 13 15]; % Index of simulated macros in which 
                                        % picos and UEs will be generated 
simulatedMacros = [10];                 % now, we are focused on simulate only one macro
nSectors = 3*length(simulatedMacros);

for i = 1:nSectors  % Inicialization of struct 's'
    s.(['S' num2str(i)]).macro.macroPos = zeros(1,2);
    s.(['S' num2str(i)]).macro.sectCenter = zeros(1,2);
    s.(['S' num2str(i)]).macro.picosPos = zeros(nPicos,2);
    s.(['S' num2str(i)]).macro.UEpos = ones(nMaxUE,2)*-1;
    %s.(['S' num2str(i)]).macro.UEgains = zeros(nMaxUE);
    s.(['S' num2str(i)]).macro.InterfMacroIndex = zeros(4,1);
    s.(['S' num2str(i)]).macro.autoIndex = 0;  
    s.(['S' num2str(i)]).macro.BScapacity = 0;
    s.(['S' num2str(i)]).macro.nUEsBlocking = 0;
    for j = 1:nPicos
        s.(['S' num2str(i)]).(['P' num2str(j)]).picoPos = zeros(1,2);
        s.(['S' num2str(i)]).(['P' num2str(j)]).UEpos = ones(nMaxUE,2)*-1;
        %s.(['S' num2str(i)]).(['P' num2str(j)]).UEgains = zeros(nMaxUE);
        s.(['S' num2str(i)]).(['P' num2str(j)]).BScapacity = 0;
        s.(['S' num2str(i)]).(['P' num2str(j)]).nUEsBlocking = 0;
    end
    s.(['S' num2str(i)]).X = ones(1,nPicos);
    s.(['S' num2str(i)]).X_previous = ones(1,nPicos);
end

% The position of each macro is calculated row by row 
currentX = 0;
currentY = 0;
final_index = 0;

% Row 1
macroPositions(final_index+1, :) = [currentX, VISD];
final_index =  final_index + 1;

% Row 2
currentX = currentX + ISD/2;
currentY = currentY + VISD/2;
for i=1:2
    macroPositions(final_index+i, :) = [currentX, currentY+(i-1)*VISD];
end
final_index =  final_index + 2;

% Row 3
currentX = currentX + ISD/2;
currentY = 0;
for i=1:3
    macroPositions(final_index+i, :) = [currentX, currentY+(i-1)*VISD];
end
final_index =  final_index + 3;

% Row 4
currentX = currentX + ISD/2;
currentY = VISD/2;
for i=1:2
    macroPositions(final_index+i, :) = [currentX, currentY+(i-1)*VISD];
end
final_index =  final_index + 2;

% Row 5
currentX = currentX + ISD/2;
currentY = 0;
for i=1:3
    macroPositions(final_index+i, :) = [currentX, currentY+(i-1)*VISD];
end
final_index =  final_index + 3;

% Row 6
currentX = currentX + ISD/2;
currentY = VISD/2;
for i=1:2
    macroPositions(final_index+i, :) = [currentX, currentY+(i-1)*VISD];
end
final_index =  final_index + 2;

% Row 7
currentX = currentX + ISD/2;
currentY = 0;
for i=1:3
    macroPositions(final_index+i, :) = [currentX, currentY+(i-1)*VISD];
end
final_index =  final_index + 3;

% Row 8
currentX = currentX + ISD/2;
currentY = VISD/2;
for i=1:2
    macroPositions(final_index+i, :) = [currentX, currentY+(i-1)*VISD];
end
final_index =  final_index + 2;

% Row 9
currentX = currentX + ISD/2;
currentY = VISD;
macroPositions(final_index+1, :) = [currentX, currentY];

simulatedMacroPositions = macroPositions(simulatedMacros,:);

if plotting
    figure,plot(macroPositions(:,1), macroPositions(:,2),'MarkerFaceColor','b','LineStyle','o', 'MarkerSize', 10), hold on
    title('Escenario de Simulación')
    xlabel('Distancia (m)')
    ylabel('Distancia (m)')
end

% the geographical center of each simulated sector is calculated (to be used later)
% and saved in the struct 's'.
simulatedSectorCenters =  zeros(3*length(simulatedMacros),2);
r = ISD/3;
for i = 1:length(simulatedMacros)
    simulatedSectorCenters((i-1)*3 + 1,:) = simulatedMacroPositions(i,:) + [r 0];
    s.(['S' num2str((i-1)*3+1)]).macro.macroPos = simulatedMacroPositions(i,:);
    s.(['S' num2str((i-1)*3+1)]).macro.sectCenter = simulatedSectorCenters((i-1)*3 + 1,:);
    
    simulatedSectorCenters((i-1)*3 + 2,:) = simulatedMacroPositions(i,:) + [-r*cosd(60) r*sind(60)];
    s.(['S' num2str((i-1)*3+2)]).macro.macroPos = simulatedMacroPositions(i,:);
    s.(['S' num2str((i-1)*3+2)]).macro.sectCenter = simulatedSectorCenters((i-1)*3 + 2,:);
    
    simulatedSectorCenters((i-1)*3 + 3,:) = simulatedMacroPositions(i,:) + [-r*cosd(60) -r*sind(60)];
    s.(['S' num2str((i-1)*3+3)]).macro.macroPos = simulatedMacroPositions(i,:);
    s.(['S' num2str((i-1)*3+3)]).macro.sectCenter = simulatedSectorCenters((i-1)*3 + 3,:);   
end

% The interfering macros for each sector are calculated in this loop
for i = 1:nSectors
    pos = s.(['S' num2str(i)]).macro.sectCenter;
    distances = sqrt(sum((repmat(pos,size(macroPositions,1),1) - macroPositions).^2,2));
        
    [~, interferingMacros] = sort(distances, 'ascend'); 
    ownMacro = interferingMacros(1);
    s.(['S' num2str(i)]).macro.InterfMacroIndex = sort(interferingMacros(1:4), 'ascend');
    s.(['S' num2str(i)]).macro.autoIndex = find(s.(['S' num2str(i)]).macro.InterfMacroIndex == ownMacro);
end


if plotting

    sectorCenters =  zeros(3*length(macroPositions),2);
    for i = 1:length(macroPositions)
        sectorCenters((i-1)*3 + 1,:) = macroPositions(i,:) + [r 0];
        sectorCenters((i-1)*3 + 2,:) = macroPositions(i,:) + [-r*cosd(60) r*sind(60)];
        sectorCenters((i-1)*3 + 3,:) = macroPositions(i,:) + [-r*cosd(60) -r*sind(60)];
    end

    figure,plot(simulatedSectorCenters(:,1), simulatedSectorCenters(:,2),'MarkerFaceColor','r','LineStyle','o', 'MarkerSize', 8), hold on
    plot(macroPositions(:,1), macroPositions(:,2),'MarkerFaceColor','b','LineStyle','o', 'MarkerSize', 11), hold on
    %plot(simulatedMacroPositions(:,1), simulatedMacroPositions(:,2),'MarkerFaceColor','b','LineStyle','o', 'MarkerSize', 11), hold on
    voronoi(sectorCenters(:,1), sectorCenters(:,2))
    %voronoi(simulatedSectorCenters(:,1), simulatedSectorCenters(:,2))    
    title('Escenario de Simulación')
    xlabel('Distancia (m)')
    ylabel('Distancia (m)')
end

end