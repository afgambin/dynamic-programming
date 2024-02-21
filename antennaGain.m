function [gains] =  antennaGain(macroIndex, macroPositions, UEPos)

    gains_dB = zeros(length(macroIndex), 1);
    phi_3db = 70;
    A_m = 25;
    referenceVectors = [166.6667 0; -83.3333 144.3376; -83.3333 -144.3376]; % Directions of the directional antennas.
    
    for i = 1:length(macroIndex) % for each intefering macro
        evaluationVector = UEPos - macroPositions(macroIndex(i),:); 
        
        angles = zeros(1,3);
        for j = 1:3
            % http://en.wikipedia.org/wiki/Cosine_similarity
           theta = acosd( (referenceVectors(j,:) * evaluationVector') /  (norm(evaluationVector) * norm(referenceVectors(j,:))) );
           angles(j) = min(theta, 360 - theta);
        end
        phi =  min(angles); % The angle of the interfering BS is the minimum
        gains_dB(i) = - min(12 * (phi/phi_3db).^2, A_m); 
    end
    gains = 10.^(gains_dB/10);
end