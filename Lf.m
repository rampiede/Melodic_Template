function [ Points ] = Lf( Matrix, events, M )
%Lf te devuelve las intersecciones de las funciones +-epsilon con la sweep line en el punto x
% 
%   Pregunta como hacer que solo tenga en cuenta las funciones del
%   intervalo
    x = events(:,1);
    y = events(:,2);
    z = events(:,3);
    func1 = events(:,4);
    func2 = events(:,5);
    eps1 =  events(:,6);
    eps2 =  events(:,7);
    coincident =  events(:,8);
    sweep_line = [x(1) -M x(1) M];
    Points = uniquetol(Cruces2(sweep_line, Matrix,x,y,z, func1, func2, eps1, eps2, coincident),'ByRows',true);%los de intervalo 0 vienen señalados por 0    
end

function Points  = Cruces2( XY1,F,x,y,z, func1, func2, eps1, eps2, coincident)
% intersecciones de una funcion de la forma XY1 matriz
% Nx4 donde las filas son [x1 y1 x2 y2] con (x1,y1) vertice inicial (x2,y2)
% vertice final
%[x,y] es el evento cruce
    XY2 = F(:,1:4);
    out = lineSegmentIntersect(XY1,XY2);
    P =[out.intMatrixX(:),out.intMatrixY(:)];
    Points = [];
    m = length(P);
    tol = 10^-10;
    c = length(func1);
    for j = 1:c
        for i = 1:m
            if all(abs(P(i,2) - y(j))< tol) & z(j)==0 %estamos en un cruce?
                if func1(j) == F(i,6)
                    if coincident
                        if isempty(Points)
                            Points = [Points ; P(i,:) 0 F(i,6) func2(j) eps1(j) eps2(j)];
                        else
                            t =find( Points(:,4)== func2(j)& Points(:,6)== -eps2(j));
                            t =[t find(Points(:,5)== func2(j)& Points(:,7)== -eps2(j))];
                            if ~isempty(t)&Points(t,2)<P(i,2)
                                Points = [Points ; P(i,:) 0 func2(j) F(i,6) eps2(j) eps1(j)];
                            else
                                Points = [Points ; P(i,:) 0 F(i,6) func2(j) eps1(j) eps2(j)];
                            end
                        end
                    else
                        k1 = find(F(:,1)<=x(j)+tol & F(:,6) == func1(j)&F(:,5) == eps1(j) ,1,'last');
                        k2 = find(F(:,1)<=x(j)+tol & F(:,6) == func2(j)&F(:,5) == eps2(j),1,'last');
                        if F(k1,4)< F(k2,4)
                            Points = [Points ; P(i,:) 0  func2(j) F(i,6)  eps2(j) eps1(j)];
                        else
                            Points = [Points ; P(i,:) 0 F(i,6) func2(j) eps1(j) eps2(j)];
                        end
                    end
                end
            elseif out.intAdjacencyMatrix(i) %corta con la sweep line pero no es un cruce
                diff =[];
                for k = 1:c
                    if k ~=j
                       diff = [diff abs(P(i,2) - y(k))< tol];%si no es un cruce futuro
                    end
                end                  
                if ~any(diff)|| z(j)==1
                    Points = [Points ; P(i,:) 1  F(i,6) F(i,6) F(i,5) F(i,5) ];
                end
            end
        end
    end
    Points = sortrows(Points,2);
end






