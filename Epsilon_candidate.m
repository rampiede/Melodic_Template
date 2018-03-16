function [epsilon_cand,F,perturbacion] = Epsilon_candidate(F, p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    epsilon_cand = [];
    [events, F, perturbacion] = Le_bin(F);
    m = length(events(:,1));
    Matrix = [];
    tol = 10^-10;
    n = length(F);
    for i = 1:n
        Matrix = [Matrix; F{1,i}];
    end
    L = Matrix(:,4);
    M = max([max(Matrix(:,2)),max(L)]) + 1; %NÚMERO SUFICIENTEMENTE GRANDE
    i = 1;
    epsilon0_plus = perturbacion;
    epsilon1_plus = perturbacion;
    epsilon0_minus = perturbacion;
    epsilon1_minus = perturbacion;
    for k = 1:m
        t = find(abs(events(:,1) - events(i,1)) < tol);
        e = events(i:t(end),:);
        i = t(end);
        %j2 = 1; %por si se repiten varios
        intersection = Lf_bin( Matrix, e, M );
        t = find(abs(events(i,2) - intersection(:,2)) < tol);
        len = length(intersection(:,1));
        if t + p-1 <= len
            functionp_1 = intersection(t+p-1,2);
            if t+p <= len
                functionp = intersection(t+p,2);
                epsilon0_plus = abs(events(i,2) - functionp) + perturbacion;
            end
            epsilon1_plus = abs(events(i,2) - functionp_1) + perturbacion;
        end
        if t-p+1 >= 1
            functionp_1 = intersection(t-p+1,2);
            if t-p >= 1
                functionp = intersection(t-p,2);
                epsilon0_minus = abs(events(i,2) - functionp) + perturbacion;
            end
            epsilon1_minus = abs(events(i,2) - functionp_1) + perturbacion;
        end
        epsilon_cand = [epsilon_cand epsilon0_plus epsilon1_plus epsilon0_minus epsilon1_minus];
        if i == m
            break
        end
        i = i+1;
    end
    epsilon_cand = sort(uniquetol(epsilon_cand,tol))/2;
end


function [x, F, perturbacion] = Le_bin( F )
%Le ordena los eventos que son vertices y cruces de las funciones F
    [x,F, perturbacion] = CrucesFunciones_bin(F);%con las mas menos epsilon
    x = Vertices_bin(F,x);   
    x = unique(sortrows(x), 'rows');
end

function y = Vertices_bin(F,x)
    %Encuentra los vertices de las funciones f_i de F
    % los cruces los asignamos con un indicador = 1
    m = length(F);
    y = x;
    tol = 10^-10;
    if isempty(x)||~any(abs(y(:,1) - F{1,1}(1,1)) < tol)
        y = [y; F{1,1}(1,1), F{1,1}(1,2), 1,1 1 0];%Si no pasa nada al principio ponemos uno cualquiera
    end
    if isempty(x)||~any(abs(y(:,1) - F{1,1}(end,3)) < tol)
        y = [y; F{1,1}(end,3), F{1,1}(end,4), 1,length(F) length(F) 0];%Si no pasa nada al final ponemos uno cualquiera
    end
    for i = 1:m
       n = length(F{i}(:,1));
       for j = 1:n-1
          if  ~any(abs(F{1,i}(j,3) - y(:,1)) < tol)
              if (F{i}(j,2) - F{i}(j,4) ~= 0)||(F{i}(j+1,2) - F{i}(j+1,4) ~= 0)
                  y = [y; F{1,i}(j,3), F{1,i}(j,4), 1,i i 0];
              end
          end
       end  
    end

end
function [Points]  = Cruces_bin( F,XY2 )
% intersecciones de una funcion de la forma XY1 matriz
% Nx4 donde las filas son [x1 y1 x2 y2] con (x1,y1) vertice inicial (x2,y2)
% vertice final
% los cruces los asignamos con un indicador = 0
    XY1 = F(:,1:4);
    XY2_1  = XY2(:,1:4);
    out = lineSegmentIntersect(XY1,XY2_1);
    P = [out.intMatrixX(:),out.intMatrixY(:)];
    Points = [];
    m = length(P(:,1));
    c = 0;
    for i = 1:m
        if ~isempty(find(out.coincAdjacencyMatrix,1))
           [I,list] =  find(out.coincAdjacencyMatrix);
           len = length(list);
           for k = 1:len
              J = list(k); 
              if I == length(F(:,1))
                  Points = [Points ; [XY2(J,3) XY2(J,4) 0 F(I,5) XY2(J,5) 1]];
              else
                  Points = [Points ; [XY2(J,1) XY2(J,2) 0 F(I,5) XY2(J,5) 1]];
              end
           end
        end
        if any(P(i,:) ~= [0 0])
            if ~isnan(P(i,:))
                c = c + 1;
                [listx,list] = find(out.intAdjacencyMatrix);
                J = list(c);
                I = listx (c);
                Points = [Points ; P(i,:), 0, F(I,5) XY2(J,5)  0];
            end
        end
    end
end
function [ Points,X, perturbacion ] = CrucesFunciones_bin( X )
% Calcula todos los cruces de todas las funciones
% X es una lista de m matrices Nix4 donde Ni es la cantidad de vértices que
% tiene la función contando el punto inicial.
    m = length(X);
    Points = [];
    X2 = [];
    perturbacion = 0;
    tol = 10^-10;
    for i = 2:m
       X1 = [X{i} i*ones(length(X{:,i}(:,1)),1)];
       X2 = [X2; X{i-1} (i-1)*ones(length(X{i-1}(:,1)),1)]; 
       Cr = Cruces_bin(X1, X2);
       Cr = uniquetol(Cr,tol,'ByRows',true);
       if ~isempty(Cr) 
           if ~isempty(Points)
               P = sortrows(Points);
               t = find(P(:,end)==1);
               len = length(P(1,:))-1;
               if ~isempty(t)
                   xs = P(t,1); %x inicial
                   ys = P(t,2); %y inicial
                   funcs = P(t,4:len); %funciones coincidentes
                   j = length(funcs(:,1));
                   for j2 = 1:j
                      ti = find(ismember(P(:,4:len),funcs(j2,:),'rows'));
                      ti(ti==t(j2))=[]; %elimina las dos funciones coincidentes, solo queremos cuando corta a una coincidente
                      if ~isempty(ti) %si hay alguna funcion coincidente (f1 y f2)
                          x2 =  P(ti(1),1); %x final
                          f1 = funcs(j2,1);
                          %f2 = funcs(j2,2);
                          Cr1 = Cr(:,5);%funciones que cortan con la nueva
                          tp1 = find(ismember(Cr1,f1,'rows')); %si se corta con f1
                          %tp2 = find(ismember(Cr1,f2,'rows')); %si se corta con f2
                          if ~isempty(tp1) && any(xs(j2) <= Cr(tp1,1) + tol & Cr(tp1,1) <= x2 + tol,1)% si corta con f1 en la coincidencia de f1 y f2
                              rng('default');
                              perturbacion = rand(1, 1)*0.00001+tol;
                              P1 = Points(:,5);%
                              P2 = Points(:,4);% funciones en los puntos
                              tq1 = find(ismember(P1, f1,'rows')); %posicion f1 en P
                              tq2 = find(ismember(P2, f1,'rows'));
                              tq = sort([tq1',tq2']); 
                              Points(tq,:) = [];%borro los cruces con f1 de la lista 
                              j = funcs(j2,1); %el numero de funcion f1
                              %modificar, sabiendo que los dos cruces
                              %de dos funciones coincidentes son
                              %vértices
                              t3 = find(ismember(X{j}(:,1:2),[xs(j2),ys(j2)],'rows'),1);%MIRAR
                              X{j}(t3,4) = X{j}(t3,4)+perturbacion;
                              if t3 ~= length(FMinus{j}(:,1))%si no es el ultimo tambien tengo que cambiar
                                  X{j}(t3+1,2) = X{j}(t3,4);
                              end
                              X2 = [cat(1,X{:,1:i-1}) X2(:,end)];
                              X1 = [X{i} i*ones(length(X{:,i}(:,1)),1)];
                              %calcular de nuevo los cruces con la
                              %func1
                              X1_new = X{:,j};
                              X2_new = X{:,1:j-1};
                              Cr_new = Cruces_bin(X1_new, X2_new);
                              Cr_new = uniquetol(Cr_new,tol,'ByRows',true);
                              if ~isempty(Cr_new)
                                   Points =[Points; Cr_new ];
                              end
                              Cr = Cruces_bin(X1, X2);
                              Cr = uniquetol(Cr,tol,'ByRows',true);
                          end
                      end
                   end
               end
           end
           N = length(Cr(:,1));
           for k = 1:N
               if ~isempty(Points)                         
                   if all(ismember(Cr(k,1:2),Points(:,1:2),'rows') & ~any(all(ismember( Points,Cr(k,:)),2))) % si se repite pero no es exactamente el mismo           
                       y = find(abs(Points(:,1) - Cr(k,1)) < tol & abs(Points(:,2) - Cr(k,2)) > tol );
                       rng('default');
                       perturbacion = rand(1, 1)*0.00001 + tol;
                       if ~isempty(y)
                           d = MinDistPoints(Points(y,2), Cr(k,2) );
                           while d/2 < perturbacion
                               perturbacion = perturbacion*0.1;
                           end 
                       end
                       rowmax = max([find( X{i}(:,3) <= Cr(k,1))+1;1]);
                       fmax = X{i}(rowmax,:);
                       for j = 1:m
                           t = find(all(X{j} == fmax,2));
                           if ~isempty(t)
                               X{j}(t,4) = X{j}(t,4) + perturbacion;
                               if t ~= length(X{j}(:,1))  
                                   X{j}(t+1,2) = X{j}(t,4);
                               end
                           end
                       end
                       tq1 = find(Points(:,4) == i);
                       tq2 = find(Points(:,5) == i);
                       tq = sort([tq1',tq2']);
                       Points(tq,:)=[];
                       X1 = [X{i} i*ones(length(X{:,i}(:,1)),1)];
                       X2 = [cat(1,X{:,1:i-1}) X2(:,end)];
                       Cr = Cruces_bin(X1, X2);
                       Points = cat(1,Points,Cr); 
                       break
                   end
               end
               Points =[Points; Cr(k,:) ];
           end
       end
    end
end
function [ Points ] = Lf_bin( Matrix, events, M )
%Lf te devuelve las intersecciones de las funciones +-epsilon con la sweep line en el punto x
    x = events(:,1);
    y = events(:,2);
    z = events(:,3);
    tol = 10^-10;
    sweep_line = [x(1) -M x(1) M];
    Points = uniquetol(Cruces2_bin(sweep_line, Matrix,y,z),tol,'ByRows',true);%los de intervalo 0 vienen señalados por 0    
end

function Points  = Cruces2_bin( XY1,XY2,y, z)
% intersecciones de una funcion de la forma XY1 matriz
% Nx4 donde las filas son [x1 y1 x2 y2] con (x1,y1) vertice inicial (x2,y2)
% vertice final
%[x,y] es el evento cruce
    out = lineSegmentIntersect(XY1,XY2);
    P =[out.intMatrixX(:),out.intMatrixY(:)];
    Points = [];
    m = length(P);
    tol = 10^-10;
    c = length(y);
    for j = 1:c
        for i = 1:m
            if all(abs(P(i,2) - y(j))< tol) && z(j) == 0 && out.intAdjacencyMatrix(i)
                Points = [Points ; P(i,:) 0];
            elseif out.intAdjacencyMatrix(i) 
                diff = [];
                for k = 1:c
                    if k ~= j
                       diff = [diff abs(P(i,2) - y(k))< tol];
                    end
                end                  
                if ~any(diff)|| z(j) == 1
                    Points = [Points ; P(i,:) 1];
                end
            end
        end
    end
    Points = sortrows(Points,2);
end
