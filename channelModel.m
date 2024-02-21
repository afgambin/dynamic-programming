function [gain] = channelModel(pos1, pos2, type)
distance = sqrt(sum((pos1 - pos2).^2,2)) / 1000; %in Km

    if strcmp(type, 'm2ue') % Macro to UE
        L = 128.1+37.6 * log10(distance);
    elseif strcmp(type, 'p2ue') % Pico to UE
        L = 140.71+36.7 * log10(distance);
    else
        error('Invalid argument.')
    end
    gain =  10.^(-L/10);
end