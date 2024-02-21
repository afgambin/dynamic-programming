function [activationArray] = calculateActivationCostMatrix(s,nPicos)
% Devuelve un array indicando el n�mero de estaciones base pico que se
% activan para cada uno de los estados posibles en funci�n del estado actual s (siendo s el n�mero de picos activas).
activationArray = zeros(1,nPicos+1);

for i=1:nPicos+1
    if ((i-1) - s) > 0
        activationArray(i) = i-1-s;
    end
end

end