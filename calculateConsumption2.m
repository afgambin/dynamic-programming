function [C] = calculateConsumption2(X, t, thetaBlockade, muBlockade, sigmaBlockade, thetaPico, muPico, sigmaPico, thetaMacro, muMacro, sigmaMacro,  blockingThreshold, nABSsubframes, CRE)
% Función que calcula el consumo total del sistema. Se hace uso de la
% estructura de datos obtenido en regresión lineal múltiple.

%debug
% clc, clear all, close all
% X= [1 0 1 0];
% t = 35;
% sec = 1;
% blockingThreshold = 0.2;

picoPowerSleep = 39; %[W]
nPicos = length(X);
nActivePicos = sum(X == 1); % Number of active Pico BS

% Blocking Probability
Xnorm = [1 (nActivePicos-muBlockade(1))/sigmaBlockade(1) (t-muBlockade(2))/sigmaBlockade(2) (nABSsubframes-muBlockade(3))/sigmaBlockade(3) (CRE-muBlockade(4))/sigmaBlockade(4)];
blockingProb = Xnorm*thetaBlockade;
%blockingProb = 1./(1+ exp(-(Xnorm*thetaBlockade)));

if (blockingProb > blockingThreshold)   % Si la prob de bloqueo obtenida supera el umbral fijado (parámetro de diseño), se considera un estado prohibido (indicado como -1) 
    C = -1;
else
    
    % Consumption off Picos
    offConsumption = (nPicos - nActivePicos)*picoPowerSleep;
    
    % Consumption Macro
    XnormMacro  = [1 (nActivePicos-muMacro(1))/sigmaMacro(1) (t-muMacro(2))/sigmaMacro(2) (nABSsubframes-muMacro(3))/sigmaMacro(3) (CRE-muMacro(4))/sigmaMacro(4)];
    macroConsumption = XnormMacro*thetaMacro;
    
    %Consumption active Picos
    XnormPico = [1 (nActivePicos-muPico(1))/sigmaPico(1) (t-muPico(2))/sigmaPico(2) (nABSsubframes-muPico(3))/sigmaPico(3) (CRE-muPico(4))/sigmaPico(4)];
    picoConsumption = nActivePicos*(XnormPico*thetaPico);

    C = offConsumption + macroConsumption + picoConsumption;
    
end
%fprintf('Eficiente: %f \n',toc);
end

