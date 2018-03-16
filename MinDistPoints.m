function d = MinDistPoints( Points, PlusMinusPoints )
% Calcula la distancia al punto m�s cercano dentro de Points en vertical con la funci�n lineal a trozos
% func
    m = length(Points(:,1));
    d = abs(Points(1)-PlusMinusPoints);
    for i = 2:m
       d = min(d , abs(Points(i)-PlusMinusPoints));
    end
end

