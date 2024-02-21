function [activationArray] = calculateActivationCostMatrix(s,nPicos)
% Devuelve un array indicando el número de estaciones base pico que se
% activan para cada uno de los estados posibles en función del estado actual s (siendo s el número de picos activas).
activationArray = zeros(1,nPicos+1);

for i=1:nPicos+1
    if ((i-1) - s) > 0
        activationArray(i) = i-1-s;
    end
end

end