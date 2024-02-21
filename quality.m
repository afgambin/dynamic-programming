function oneMinusAlpha = quality(relativeTol, num, sum, sumcuadrado)
% Calcula el intervalo de confianza asociado a las muestras y
% Devuelve la probabilidad asociada al intervalo de confianza 

%calcula media muestral y cuasi-varianza            
media_muestral = sum / num;
Cuasi_Var = (sumcuadrado - ((sum.^2)/num))/(num-1);
indexes = Cuasi_Var < 0;
Cuasi_Var(indexes) = 0;

%tolerancia total
    tol = relativeTol*media_muestral;

%Obtener rango del intervalo
    z_a = tol./(sqrt(Cuasi_Var/num));

%calcular calidad 1-alpha

    alfa = (1 - normcdf(z_a))*2;
    oneMinusAlpha = 1 - alfa;
end