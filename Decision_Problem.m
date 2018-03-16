function [ Answer,FPlusMinusEpsilonLabel,FPlus,FMinus,Events, boolean ] = Decision_Problem(F, epsilon, p)
% Decision_Problem resuelve el problema de decisión del template, dado m1 funciones lineales a trozos decidir si puedo obtener un tubo de amplitud 2epsilon
% tal que en todo momento contenga a p funciones (o %p), esas p funciones pueden ir cambiando no tienen que ser siempre las mismas.
%% preprocesado
    m1 = length(F);% # of functions
    %p = ceil(p/m); porcentaje
    FPlusMinusEpsilon = F;
    FPlusMinusEpsilonLabel = F;
    FPlus  = F;
    FMinus = F;
    for i = 1:m1 
        FPlusMinusEpsilonLabel{1,i} = [F{1,i}(:,1) F{1,i}(:,2)+epsilon F{1,i}(:,3) F{1,i}(:,4)+epsilon ones(length(F{1,i}(:,1)),1) i*ones(length(F{1,i}(:,1)),1);
                                      F{1,i}(:,1) F{1,i}(:,2)-epsilon F{1,i}(:,3) F{1,i}(:,4)-epsilon -ones(length(F{1,i}(:,1)),1) i*ones(length(F{1,i}(:,1)),1)];
        FPlusMinusEpsilon{1,i}      = [F{1,i}(:,1) F{1,i}(:,2)+epsilon F{1,i}(:,3) F{1,i}(:,4)+epsilon;
                                      F{1,i}(:,1) F{1,i}(:,2)-epsilon F{1,i}(:,3) F{1,i}(:,4)-epsilon];
        FPlus{1,i}  = [F{1,i}(:,1) F{1,i}(:,2)+epsilon F{1,i}(:,3) F{1,i}(:,4)+epsilon];
        FMinus{1,i} = [F{1,i}(:,1) F{1,i}(:,2)-epsilon F{1,i}(:,3) F{1,i}(:,4)-epsilon];
    end
    
    [Events,FPlusMinusEpsilonLabel,FPlus, FMinus] = Le(F,FPlusMinusEpsilonLabel, FPlusMinusEpsilon,FPlus,FMinus);
    m = length(Events(:,1));
    Matrix = [];
    for i = 1:length(FPlusMinusEpsilonLabel)
        Matrix = [Matrix; FPlusMinusEpsilonLabel{1,i}];
    end
    L  = Matrix(:,4);
    M  = max([max(Matrix(:,2)),max(L)]) + 1; %a big M
    l2 = m1*2 -1;% # of intervals
    boolean = zeros(m, l2);
    Answer  = 0;
    tol = 10^-12;
%     c = distinguishable_colors(m1+1);%colors
%     for i = 1:m1 % Si quieres ver visualemente las funciones
%         col1 = rand;
%         col2 = rand;
%         col3 = rand;
% 
%         plot([F{1,i}(:,1),F{1,i}(:,3)],[F{1,i}(:,2),F{1,i}(:,4)],'x-','Color', c(i,:),'LineWidth',3);
%         hold on; 
%         plot([F{1,i}(:,1),F{1,i}(:,3)],[F{1,i}(:,2)+epsilon,F{1,i}(:,4)+epsilon],'x--','Color',c(i,:));
%         hold on; 
%         plot([F{1,i}(:,1),F{1,i}(:,3)],[F{1,i}(:,2)-epsilon,F{1,i}(:,4)-epsilon],'x--','Color',c(i,:));
%         hold on; 
%     end
      
    %% caso i = 1
    t = find(abs(Events(:,1)- Events(1,1))<tol);
    e = Events(1:t(end),:);
    intersection  = Lf(Matrix,  e,M);
    intersection0 = intersection;
    position = [];
    n  = zeros(l2,1);   
    j2 = 1; %j2 cuenta los vertices repetidos (los de longitud 0)
    aux  = 0;
    aux2 = 1;
    r2 = [];
    q2 =[];
    for j = 1:l2
        if intersection(j2,3) ~= 0 || aux == 1 %si este punto no es el cruce    
            int = [intersection(j2,2), intersection(j2+1,2)]; 
            if aux
                position = [position; position(j-1,2) intersection(j2+1,4) position(j-1,4) intersection(j2+1,6)];
            else
                position = [position; intersection(j2,5) intersection(j2+1,4) intersection(j2,7) intersection(j2+1,6)];
            end
            j2 = j2+1;
            aux = 0;
        else %si en este evento j es un cruce
            int = [intersection(j2,2), intersection(j2,2)];
            k1  = intersection(j2,4);
            k2  = intersection(j2,5);
            t1  = find(FPlusMinusEpsilonLabel{k1}(:,5) == intersection(j2,6),1);
            t2  = find(FPlusMinusEpsilonLabel{k2}(:,5) == intersection(j2,7),1);
            position = [position; intersection(j2,4:7) ]; %MAL LA POSICION
            intersection0 =[intersection0(1:j2+(aux2-1),:); intersection0(j2+(aux2-1):end,:)];
            if FPlusMinusEpsilonLabel{k1}(t1,4) >= FPlusMinusEpsilonLabel{k2}(t2,4) || abs(FPlusMinusEpsilonLabel{k1}(t1,4) - FPlusMinusEpsilonLabel{k2}(t2,4)) < tol
                if j ~= 1
                    position(j-1,:) = [position(j-1,1) intersection(j2,5) position(j-1,3) intersection(j2,7)];
                end
                position(j,:) =  [ intersection(j2,5) intersection(j2,4) intersection(j2,7) intersection(j2,6)] ;
                intersection0(j2+1+(aux2-1),:) = [intersection0(j2+1+(aux2-1),1:2) 1 intersection0(j2+1+(aux2-1),4) intersection0(j2+1+(aux2-1),4) intersection0(j2+1+(aux2-1),6) intersection0(j2+1+(aux2-1),6)];
                intersection0(j2+(aux2-1),:)   = [intersection0(j2+(aux2-1),1:2) 1 intersection0(j2+(aux2-1),5) intersection0(j2+(aux2-1),5) intersection0(j2+(aux2-1),7) intersection0(j2+(aux2-1),7)];

            else
                intersection0(j2+(aux2-1),:)   = [intersection0(j2+(aux2-1),1:2) 1 intersection0(j2+(aux2-1),4) intersection0(j2+(aux2-1),4) intersection0(j2+(aux2-1),6) intersection0(j2+(aux2-1),6)];
                intersection0(j2+1+(aux2-1),:) = [intersection0(j2+1+(aux2-1),1:2) 1 intersection0(j2+(aux2-1)+1,5) intersection0(j2+(aux2-1)+1,5) intersection0(j2+(aux2-1)+1,7) intersection0(j2+(aux2-1)+1,7)];

            end
            aux = 1;
            r2  = [r2 intersection(j2,6)];
            q2  = [q2 intersection(j2,7)];
            
            aux2 = aux2 + 1;
            
        end
        n(j) = Count2(int, FPlus, FMinus,intersection, j2, M ); 
        if n(j) >= p
            boolean(1,j) = 1;
        end 
    end
    intersection2 = intersection;
    n2 = n;
    
    %% casos k = 2:m
    i = length(e(:,1)) + 1;
    for k = 2:m
        d = e;
        t = find(abs(Events(:,1) - Events(i,1))<tol);
        e = Events(i:t(end),:);
        i = t(end); 
        aux  = 0;
        aux2 = 0;
        h = 1;
        if all(e(:,3) == 1) && all(d(:,3) == 1) || all(e(:,3) == 1) && all(d(:,end) == 1)% el evento es un vertice y el anterior tambien|| es coincidente y hay un vertice
            boolean(k,:) = boolean(k-1,:);
            e = d;
        else %el evento es un cruce || es un vertice pero el anterior es un cruce
            [intersection,r,q, position] = Lf2(FPlusMinusEpsilonLabel,e, intersection0, position, d);
            j2 = 1; %j2 cuenta los vertices repetidos (los de longitud 0)
            j3 = 1;
            for j = 1:l2
                if intersection(j2,3) ~= 0 || aux == 1%si este punto no es el cruce
                    j2  = j2 + 1;
                    aux = 0;
                    if intersection2(j3,3) == 0 && aux2 == 0%si en el evento pasado el punto j era cruce  
                        aux2 = 1;
                        if r2(h) == 1 && q2(h) == 1
                            h    = h+1;
                            n(j) = n2(j) - 1;
                        elseif r2(h) == 1 && q2(h) == -1
                            h    = h+1;
                            n(j) = n2(j);
                        elseif r2(h) == -1 && q2(h) == -1
                            h    = h+1;
                            n(j) = n2(j) - 1;
                        else
                            h    = h+1;
                            n(j) = n2(j) - 2;
                        end
                    else %si en el evento pasado el punto j no era cruce
                        aux2 = 0;
                        j3   = j3+1;
                        n(j) = n2(j);
                    end
                else %si en este evento j es un cruce
                    
                    aux = 1;
                    if intersection2(j3,3) == 0 && aux2 == 0 %el evento pasado cruce
                        n(j) = n2(j);
                        aux2 = 1;
                    else
                        j3   = j3+1;
                        int  = [intersection(j2,2), intersection(j2,2)];
                        n(j) = Count2(int, FPlus, FMinus,intersection, j2, M);
                    end
                end 

                if n(j) >= p
                    if j == 1
                        boolean(k,j) = boolean(k-1,j)|boolean(k-1,j+1);
                    elseif j == l2
                        boolean(k,j) = boolean(k-1,j)|boolean(k-1,j-1);
                    else
                        boolean(k,j) = boolean(k-1,j)|boolean(k-1,j-1)|boolean(k-1,j+1); 
                    end
                end 
            end
            intersection2 = intersection;
            n2 = n;
            r2 = r;
            q2 = q;
        end
        if i == m && any(boolean(k,:))
            Answer = 1;
            break;
        end
        if boolean(k,:) == zeros(1,l2) % si en hasta el evento i no se puede llegar no hace falta seguir
            break;
        end
        i = i+1;
    end
   
end

function [x,FPlusMinusEpsilonLabel,FPlus, FMinus] = Le( F,FPlusMinusEpsilonLabel, FPlusMinusEpsilon,FPlus, FMinus )
%Le ordena los eventos que son vertices y cruces de las funciones F
    [x,FPlusMinusEpsilonLabel,FPlus, FMinus] = CrucesFunciones(F,FPlusMinusEpsilonLabel,FPlusMinusEpsilon,FPlus, FMinus);%con las mas menos epsilon
    x = Vertices(FPlusMinusEpsilonLabel,x);   
    x = unique(sortrows(x), 'rows');
end

function y = Vertices(F,x)
    %Encuentra los vertices de las funciones f_i de F
    % los cruces los asignamos con un indicador = 1
    m = length(F);
    y = x;
    tol = 10^-12;
    if isempty(x)||~any(abs(y(:,1) - F{1,1}(1,1))<tol)
        y = [y; F{1,1}(1,1), F{1,1}(1,2), 1,F{1,1}(1,6), F{1,1}(1,6), F{1,1}(1,5), F{1,1}(1,5) 0];%Si no pasa nada al principio ponemos uno cualquiera
    end
    if isempty(x)||~any(abs(y(:,1) - F{1,1}(end,3))<tol)
        y = [y; F{1,1}(end,3), F{1,1}(end,4), 1,F{1,1}(end,6), F{1,1}(end,6), F{1,1}(end,5), F{1,1}(end,5) 0];%Si no pasa nada al final ponemos uno cualquiera
    end
    for i = 1:m
       n = length(F{i}(:,1));
       for j = 1:n-1
          if  ~any(abs(F{1,i}(j,3) - y(:,1))<tol)%para no repetir vertices y cruces
              if (F{i}(j,2) - F{i}(j,4)~=0)||(F{i}(j+1,2) - F{i}(j+1,4)~=0)
                  y = [y; F{1,i}(j,3), F{1,i}(j,4), 1,F{1,i}(j,6), F{1,i}(j,6), F{1,i}(j,5), F{1,i}(j,5) 0];
              end
          end
       end  
    end

end
function Points  = Cruces( F,XY2, X2_label )
% intersecciones de una funcion de la forma XY1 matriz
% Nx4 donde las filas son [x1 y1 x2 y2] con (x1,y1) vertice inicial (x2,y2)
% vertice final
% los cruces los asignamos con un indicador = 0
    XY1 = F(:,1:4);
    out = lineSegmentIntersect(XY1,XY2);
    P = [out.intMatrixX(:),out.intMatrixY(:)];
    Points = [];
    m = length(P(:,1));
    tol = 10^-12;
    c = 0;
    for i = 1:m
        if ~isempty(find(out.coincAdjacencyMatrix,1))
           [I,list] =  find(out.coincAdjacencyMatrix);
           len = length(list);
           for k = 1:len
              J = list(k); 
               
              if I == length(F(:,1))
                  Points = [Points ; [XY2(J,3) XY2(J,4) 0 F(1,6) X2_label(J,6) F(1,5) X2_label(J,5) 1]];
              else
                  t1 = find(F(:,2)>= XY2(J,2),1);
                  t2 = find(X2_label(:,2)>= XY2(J,2),1);
                  if F(t1,2) > X2_label(t2,2)+tol
                    Points = [Points ; [XY2(J,1) XY2(J,2) 0 F(1,6) X2_label(J,6) F(1,5) X2_label(J,5) 1]];
                  elseif abs(F(t1,2) - X2_label(t2,2))<tol
                      if F(t1,4) > X2_label(t2,4)
                          Points = [Points ; [XY2(J,1) XY2(J,2) 0 F(1,6) X2_label(J,6) F(1,5) X2_label(J,5) 1]];
                      else
                          Points = [Points ; [XY2(J,1) XY2(J,2) 0  X2_label(J,6) F(1,6)  X2_label(J,5) F(1,5) 1]];
                      end
                  else
                      Points = [Points ; [XY2(J,1) XY2(J,2) 0  X2_label(J,6) F(1,6)  X2_label(J,5) F(1,5) 1]];
                  end
              end
           end
        end
        if any(P(i,:) ~= [0 0])
            if ~isnan(P(i,:))
                c = c + 1;
                [listx,list] = find(out.intAdjacencyMatrix);
                J = list(c);
                I = listx (c);
                Points = [Points ; P(i,:), 0, F(I,6) X2_label(J,6)  F(I,5) X2_label(J,5) 0];
            end
        end
    end
end
function [ Points,FPlusMinusEpsilonLabel,FPlus, FMinus ] = CrucesFunciones( X,FPlusMinusEpsilonLabel,FPlusMinusEpsilon,FPlus, FMinus )
% Calcula todos los cruces de todas las funciones
% X es una lista de m matrices Nix4 donde Ni es la cantidad de vértices que
% tiene la función contando el punto inicial.
    m = length(X);
    Points = [];
    X2 = [];
    X2_label = [];
    tol = 10^-12;
    for i = 2:m
       aux = 0;
       L = FPlusMinusEpsilonLabel{:,i};
       X1Plus = L(1:length(L(:,1))/2,:);
       X1Minus = L(length(L(:,1))/2+1:end,:);
       X2 = [X2; FPlusMinusEpsilon{:,i-1}]; 
       X2_label = [X2_label; FPlusMinusEpsilonLabel{:,i-1}];
       %X2 = X{:,i+1:m};
       CPlus = Cruces(X1Plus, X2, X2_label);
       CPlus = uniquetol(CPlus,'ByRows',true);
       CMinus = Cruces(X1Minus, X2, X2_label);
       CMinus = uniquetol(CMinus,'ByRows',true);
       if ~isempty(CPlus) 
           if ~isempty(Points)
               P = sortrows(Points);
               t = find(P(:,end) == 1);
               if ~isempty(t)
                   xs = P(t,1); %x inicial
                   ys = P(t,2); %y inicial
                   funcs = P(t,4:7); %funciones coincidentes
                   j = length(funcs(:,1));
                   for j2 = 1:j
                      ti = find(ismember(P(:,4:7),funcs(j2,:),'rows'));
                      ti = [ti find(ismember(P(:,4:7),[funcs(j2,2) funcs(j2,1) funcs(j2,4) funcs(j2,3)],'rows'))];
                      ti(ti == t(j2)) = []; %elimina las dos funciones coincidentes, solo queremos cuando corta a una coincidente
                      if ~isempty(ti) %si hay alguna funcion coincidente (f1 y f2)
                          x2 =  P(ti(1),1); %x final
                          %f1 = [funcs(j2,1) funcs(j2,3)];
                          f2 = [funcs(j2,2) funcs(j2,4)];
                          CPlus1 = [CPlus(:,5) CPlus(:,7)];%funciones que cortan con la nueva
                          %tp1 = find(ismember(CPlus1,f1,'rows')); %si se corta con f1
                          tp2 = find(ismember(CPlus1,f2,'rows')); %si se corta con f2
                          if all(~isempty(tp2) & any(xs(j2)<=CPlus(tp2,1)& CPlus(tp2,1)<=x2,1))% si corta con f1 en la coincidencia de f1 y f2
                              rng('default');
                              perturbacion = rand(1, 1)*0.09;
                              P1 = Points(:,5);%
                              P2 = Points(:,4);% funciones en los puntos
                              tq1 = find(ismember(P1, f2(1),'rows')); %posicion f1 en P
                              tq2 = find(ismember(P2, f2(1),'rows'));
                              tq = sort([tq1',tq2']); 
                              Points(tq,:) = [];%borro los cruces con f1 de la lista 
                              j = funcs(j2,2); %el numero de funcion f1
                              %modificar, sabiendo que los dos cruces
                              %de dos funciones coincidentes son
                              %vértices
                              eps1 = funcs(j2,4);
                              %modificar, sabiendo que los dos cruces
                              %de dos funciones coincidentes son
                              %vértices
                              t3 = find(ismember(FPlusMinusEpsilon{j}(:,1:2),[xs(j2),ys(j2)],'rows'),1);
                              if eps1 == -1
                                  FPlusMinusEpsilonLabel{j}(t3,4) = FPlusMinusEpsilonLabel{j}(t3,4)+perturbacion;
                                  FPlusMinusEpsilonLabel{j}(floor(t3/2),4) = FPlusMinusEpsilonLabel{j}(floor(t3/2),4)+perturbacion ;
                                  FPlusMinusEpsilon{j}(t3,4) = FPlusMinusEpsilon{j}(t3,4)+perturbacion;
                                  FPlusMinusEpsilon{j}(floor(t3/2),4) = FPlusMinusEpsilon{j}(floor(t3/2),4)+perturbacion ;
                                  FPlus{j}(floor(t3/2),4) = FPlus{j}(floor(t3/2),4)+perturbacion;
                                  FMinus{j}(floor(t3/2),4) = FMinus{j}(floor(t3/2),4)+perturbacion;
                                  X{j}(floor(t3/2),4) = X{j}(floor(t3/2),4)+perturbacion;
                                  if t3 ~= length(FMinus{j}(:,1))%si no es el ultimo tambien tengo que cambiar
                                      FPlusMinusEpsilonLabel{j}(t3+1,2) = FPlusMinusEpsilonLabel{j}(t3,4);
                                      FPlusMinusEpsilonLabel{j}(floor(t3/2)+1,2) = FPlusMinusEpsilonLabel{j}(floor(t3/2),4) ;
                                      FPlusMinusEpsilon{j}(t3+1,2) = FPlusMinusEpsilon{j}(t3,4);
                                      FPlusMinusEpsilon{j}(floor(t3/2)+1,2) = FPlusMinusEpsilon{j}(floor(t3/2),4) ;
                                      FPlus{j}(floor(t3/2)+1,2) = FPlus{j}(floor(t3/2),4);
                                      FMinus{j}(floor(t3/2)+1,2) = FMinus{j}(floor(t3/2),4);
                                      X{j}(floor(t3/2)+1,2) = X{j}(floor(t3/2),4);
                                  end
                              else
                                  L = length(FPlusMinusEpsilonLabel{j}(:,1))/2;
                                  FPlusMinusEpsilonLabel{j}(t3,4) = FPlusMinusEpsilonLabel{j}(t3,4)+perturbacion;
                                  FPlusMinusEpsilonLabel{j}(t3+L,4) = FPlusMinusEpsilonLabel{j}(t3+L,4)+perturbacion ;
                                  FPlusMinusEpsilon{j}(t3+L,4) = FPlusMinusEpsilon{j}(t3+L,4)+perturbacion;
                                  FPlusMinusEpsilon{j}(t3,4) = FPlusMinusEpsilon{j}(t3,4)+perturbacion ;
                                  FPlus{j}(t3,4) = FPlus{j}(t3,4)+perturbacion;
                                  FMinus{j}(t3,4) = FMinus{j}(t3,4)+perturbacion;
                                  X{j}(t3,4) = X{j}(t3,4)+perturbacion;
                                  if t3 ~= length(FMinus{j}(:,1))%si no es el ultimo tambien tengo que cambiar
                                      FPlusMinusEpsilonLabel{j}(t3+L+1,2) = FPlusMinusEpsilonLabel{j}(t3+L,4);
                                      FPlusMinusEpsilonLabel{j}(t3+1,2) = FPlusMinusEpsilonLabel{j}(t3,4) ;
                                      FPlusMinusEpsilon{j}(t3+L+1,2) = FPlusMinusEpsilon{j}(t3+L,4);
                                      FPlusMinusEpsilon{j}(t3+1,2) = FPlusMinusEpsilon{j}(t3,4) ;
                                      FPlus{j}(t3+1,2) = FPlus{j}(t3,4);
                                      FMinus{j}(t3+1,2) = FMinus{j}(t3,4);
                                      X{j}(t3+1,2) = X{j}(t3,4);
                                  end
                              end
                              X2 = cat(1,FPlusMinusEpsilon{:,1:i-1});
                              X2_label =  cat(1,FPlusMinusEpsilonLabel{:,1:i-1});
                              %calcular de nuevo los cruces con la
                              %func1
                              X1_new = FPlusMinusEpsilonLabel{:,j};
                              X2_new = FPlusMinusEpsilon{:,1:j-1};
                              X2_new_label = FPlusMinusEpsilonLabel{:,1:j-1};
                              X1Plus_new = X1_new(1:length(X1_new(:,1))/2,:);
                              X1Minus_new = X1_new(length(X1_new(:,1))/2+1:end,:);
                              CPlus_new = Cruces(X1Plus_new, X2_new, X2_new_label);
                              CPlus_new = uniquetol(CPlus_new,'ByRows',true);
                              CMinus_new = Cruces(X1Minus_new, X2_new, X2_new_label);
                              CMinus_new = uniquetol(CMinus_new,'ByRows',true);
                              if ~isempty(CPlus_new)
                                   Points =[Points; CPlus_new ];
                              end
                               if ~isempty(CMinus_new)
                                   Points = [Points; CMinus_new ];
                              end
                              CPlus = Cruces(X1Plus, X2, X2_label);
                              CPlus = uniquetol(CPlus,'ByRows',true);
                              CMinus = Cruces(X1Minus, X2, X2_label);
                              CMinus = uniquetol(CMinus,'ByRows',true);
                          end
                      end
                   end
               end
           end
           N = length(CPlus(:,1));
           for k = 1:N
               if ~isempty(Points)                         
                   if all(any(all(abs(Points(:,1:2)-CPlus(k,1:2))<tol,2))& ~any(all(ismember( Points,CPlus(k,:)),2))) % si se repite pero no es exactamente el mismo           
                       y = find(abs(Points(:,1)-CPlus(k,1))<tol&abs(Points(:,2)-CPlus(k,2))>tol );
                       rng('default');
                       perturbacion = rand(1, 1)*0.09;
                       if ~isempty(y)
                           d = MinDistPoints(Points(y,2), CPlus(k,2) );
                           while d/2 < perturbacion
                               perturbacion = perturbacion*0.1;
                           end 
                       end
                       rowmax = max([find(FPlusMinusEpsilonLabel{i}(:,5) == CPlus(k,6) & FPlusMinusEpsilonLabel{i}(:,3) < CPlus(k,1))+1;1]);
                       fmax = FPlusMinusEpsilonLabel{i}(rowmax,:);
%                        t2 = find(all(abs(X2_label- fmax) < tol,2));
%                        if  ~isempty(t2)
%                            X2_label(t2,4) = X2_label(t2,4) + perturbacion;
%                        end
                       for j = 1:m
                           t = find(all(abs(FPlusMinusEpsilonLabel{j} - fmax) < tol,2));
                           if ~isempty(t)
                               FPlusMinusEpsilonLabel{j}(t,4) = fmax(4)+perturbacion;
                               FPlusMinusEpsilonLabel{j}(2*t,4) = FPlusMinusEpsilonLabel{j}(2*t,4)+perturbacion ;
                               FPlusMinusEpsilon{j}(t,4) = fmax(4)+perturbacion;
                               FPlusMinusEpsilon{j}(2*t,4) = FPlusMinusEpsilon{j}(2*t,4)+perturbacion ;
                               FPlus{j}(t,4) = fmax(4)+perturbacion;
                               FMinus{j}(t,4) = FMinus{j}(t,4)+perturbacion;
                               X{j}(t,4) = X{j}(t,4)+perturbacion;
                               if t ~= length(X{j}(:,1))
                                   FPlusMinusEpsilonLabel{j}(t+1,2)   = FPlusMinusEpsilonLabel{j}(t,4);
                                   FPlusMinusEpsilonLabel{j}(2*t+1,2) = FPlusMinusEpsilonLabel{j}(2*t,4) ;
                                   FPlusMinusEpsilon{j}(t+1,2)   = FPlusMinusEpsilon{j}(t,4);
                                   FPlusMinusEpsilon{j}(2*t+1,2) = FPlusMinusEpsilon{j}(2*t,4) ;
                                   FPlus{j}(t+1,2)  = FPlus{j}(t,4);
                                   FMinus{j}(t+1,2) = FMinus{j}(t,4);
                                   X{j}(t+1,2) = X{j}(t,4);
                               end
                           end
                       end
                       tq1 = find(Points(:,4) == i);
                       tq2 = find(Points(:,5) == i);
                       tq = sort([tq1',tq2']);
                       Points(tq,:) = [];
                       L = FPlusMinusEpsilonLabel{:,i};
                       X1Plus = L(1:length(L(:,1))/2,:);
                       X1Minus = L(length(L(:,1))/2+1:end,:);
                       X2 = cat(1,FPlusMinusEpsilon{:,1:i-1});
                       X2_label = cat(1,FPlusMinusEpsilonLabel{:,1:i-1});
                       CPlus = Cruces(X1Plus, X2, X2_label);
                       CMinus = Cruces(X1Minus, X2, X2_label);
                       Points = cat(1,Points,CPlus); 
                       break
                   end
               end
               Points =[Points; CPlus(k,:) ];
           end
       end
       if ~isempty(CMinus)
           if ~isempty(Points)
               P = sortrows(Points);
               t = find(P(:,end) == 1);
               if ~isempty(t)
                   xs = P(t,1); %x inicial
                   ys = P(t,2); %y inicial
                   funcs = P(t,4:7); %funciones coincidentes
                   j = length(funcs(:,1));
                   for j2 = 1:j
                      ti = find(ismember(P(:,4:7),funcs(j2,:),'rows'));
                      ti = [ti find(ismember(P(:,4:7),[funcs(j2,2) funcs(j2,1) funcs(j2,4) funcs(j2,3)],'rows'))];
                      ti(ti == t(j2)) = []; %elimina las dos funciones coincidentes, solo queremos cuando corta a una coincidente
                      if ~isempty(ti) %si hay alguna funcion coincidente (f1 y f2)
                          x2 =  P(ti(1),1); %x final
                          %f1 = [funcs(j2,1) funcs(j2,3)];
                          f2 = [funcs(j2,2) funcs(j2,4)];
                          CMinus1 = [CMinus(:,5) CMinus(:,7)];%funciones que cortan con la nueva
                          %tp1 = find(ismember(CMinus1,f1,'rows')); %si se corta con f1
                          tp2 = find(ismember(CMinus1,f2,'rows')); %si se corta con f2
                          if all(~isempty(tp2) & any(xs(j2) <= CMinus(tp2,1)& CMinus(tp2,1) <= x2,1))% si corta con f1 en la coincidencia de f1 y f2
                              rng('default');
                              perturbacion = rand(1, 1)*0.09;
                              P1 = Points(:,5);%
                              P2 = Points(:,4);% funciones en los puntos
                              tq1 = find(ismember(P1, f2(1),'rows')); %posicion f2 en P
                              tq2 = find(ismember(P2, f2(1),'rows'));
                              tq = sort([tq1',tq2']); 
                              Points(tq,:)=[];%borro los cruces con f1 de la lista 
                              j = funcs(j2,2); %el numero de funcion f2
                              eps1 = funcs(j2,4);
                              %modificar, sabiendo que los dos cruces
                              %de dos funciones coincidentes son
                              %vértices
                              t3 = find(ismember(FPlusMinusEpsilon{j}(:,1:2),[xs(j2),ys(j2)],'rows'),1);
                              if eps1 == -1
                                  FPlusMinusEpsilonLabel{j}(t3,4) = FPlusMinusEpsilonLabel{j}(t3,4)+perturbacion;
                                  FPlusMinusEpsilonLabel{j}(floor(t3/2),4) = FPlusMinusEpsilonLabel{j}(floor(t3/2),4)+perturbacion ;
                                  FPlusMinusEpsilon{j}(t3,4) = FPlusMinusEpsilon{j}(t3,4)+perturbacion;
                                  FPlusMinusEpsilon{j}(floor(t3/2),4) = FPlusMinusEpsilon{j}(floor(t3/2),4)+perturbacion ;
                                  FPlus{j}(floor(t3/2),4) = FPlus{j}(floor(t3/2),4)+perturbacion;
                                  FMinus{j}(floor(t3/2),4) = FMinus{j}(floor(t3/2),4)+perturbacion;
                                  X{j}(floor(t3/2),4) = X{j}(floor(t3/2),4)+perturbacion;
                                  if t3 ~= length(FMinus{j}(:,1))%si no es el ultimo tambien tengo que cambiar
                                      FPlusMinusEpsilonLabel{j}(t3+1,2) = FPlusMinusEpsilonLabel{j}(t3,4);
                                      FPlusMinusEpsilonLabel{j}(floor(t3/2)+1,2) = FPlusMinusEpsilonLabel{j}(floor(t3/2),4) ;
                                      FPlusMinusEpsilon{j}(t3+1,2) = FPlusMinusEpsilon{j}(t3,4);
                                      FPlusMinusEpsilon{j}(floor(t3/2)+1,2) = FPlusMinusEpsilon{j}(floor(t3/2),4) ;
                                      FPlus{j}(floor(t3/2)+1,2) = FPlus{j}(floor(t3/2),4);
                                      FMinus{j}(floor(t3/2)+1,2) = FMinus{j}(floor(t3/2),4);
                                      X{j}(floor(t3/2)+1,2) = X{j}(floor(t3/2),4);
                                  end
                              else
                                  L = length(FPlusMinusEpsilonLabel{j}(:,1))/2;
                                  FPlusMinusEpsilonLabel{j}(t3,4) = FPlusMinusEpsilonLabel{j}(t3,4)+perturbacion;
                                  FPlusMinusEpsilonLabel{j}(t3+L,4) = FPlusMinusEpsilonLabel{j}(t3+L,4)+perturbacion ;
                                  FPlusMinusEpsilon{j}(t3+L,4) = FPlusMinusEpsilon{j}(t3+L,4)+perturbacion;
                                  FPlusMinusEpsilon{j}(t3,4) = FPlusMinusEpsilon{j}(t3,4)+perturbacion ;
                                  FPlus{j}(t3,4) = FPlus{j}(t3,4)+perturbacion;
                                  FMinus{j}(t3,4) = FMinus{j}(t3,4)+perturbacion;
                                  X{j}(t3,4) = X{j}(t3,4)+perturbacion;
                                  if t3 ~= length(FMinus{j}(:,1))%si no es el ultimo tambien tengo que cambiar
                                      FPlusMinusEpsilonLabel{j}(t3+L+1,2) = FPlusMinusEpsilonLabel{j}(t3+L,4);
                                      FPlusMinusEpsilonLabel{j}(t3+1,2) = FPlusMinusEpsilonLabel{j}(t3,4) ;
                                      FPlusMinusEpsilon{j}(t3+L+1,2) = FPlusMinusEpsilon{j}(t3+L,4);
                                      FPlusMinusEpsilon{j}(t3+1,2) = FPlusMinusEpsilon{j}(t3,4) ;
                                      FPlus{j}(t3+1,2) = FPlus{j}(t3,4);
                                      FMinus{j}(t3+1,2) = FMinus{j}(t3,4);
                                      X{j}(t3+1,2) = X{j}(t3,4);
                                  end
                              end
                              X2 = cat(1,FPlusMinusEpsilon{:,1:i-1});
                              X2_label =  cat(1,FPlusMinusEpsilonLabel{:,1:i-1});
                              %calcular de nuevo los cruces con la
                              %func1
                              X1_new = FPlusMinusEpsilonLabel{:,j};
                              X2_new = FPlusMinusEpsilon{:,1:j-1};
                              X2_new_label = FPlusMinusEpsilonLabel{:,1:j-1};
                              X1Plus_new = X1_new(1:length(X1_new(:,1))/2,:);
                              X1Minus_new = X1_new(length(X1_new(:,1))/2+1:end,:);
                              CPlus_new = Cruces(X1Plus_new, X2_new, X2_new_label);
                              CPlus_new = uniquetol(CPlus_new,'ByRows',true);
                              CMinus_new = Cruces(X1Minus_new, X2_new, X2_new_label);
                              CMinus_new = uniquetol(CMinus_new,'ByRows',true);
                              if ~isempty(CPlus_new)
                                   Points = [Points; CPlus_new ];
                              end
                               if ~isempty(CMinus_new)
                                   Points = [Points; CMinus_new ];
                              end
                              CPlus = Cruces(X1Plus, X2, X2_label);
                              CPlus = uniquetol(CPlus,'ByRows',true);
                              CMinus = Cruces(X1Minus, X2, X2_label);
                              CMinus = uniquetol(CMinus,'ByRows',true);
                              aux = 1;
                          end
                      end
                   end
               end
           end
           N = length(CMinus(:,1));
           for k = 1:N
               if ~isempty(Points)                   
                   if all(any(all(abs(Points(:,1:2) - CMinus(k,1:2)) < tol,2)) & ~any(all(abs(Points- CMinus(k,:))<tol,2)) )
                       y = find(abs(Points(:,1) - CMinus(k,1)) < tol & abs(Points(:,2)-CMinus(k,2)) > tol ); % asegurarnos que la perturbacion es pequeña
                       rng('default');
                       perturbacion = rand(1, 1)*0.09;
                       if ~isempty(y)
                           d = MinDistPoints(Points(y,2), CMinus(k,2) );
                           while d/2<perturbacion
                               perturbacion = perturbacion*0.1;
                           end 
                       end
                       s = length(FPlusMinusEpsilonLabel{i}(:,1))/2;
                       rowmax = max([find(FPlusMinusEpsilonLabel{i}(:,5) == CMinus(k,6) & FPlusMinusEpsilonLabel{i}(:,3) < CMinus(k,1))+1;s+1]);%calcula la maxima fila de F que tiene el mismo signo y es menor que x de CMINUS 
                       fmax = FPlusMinusEpsilonLabel{i}(rowmax,:);%la fila en cuestion
%                        t2 = find(all(abs(X2_label - fmax) < tol,2));
%                        if  ~isempty(t2)
%                            X2_label(t2,4) = X2_label(t2,4) + perturbacion;
%                        end
                       for j = 1:m
                           t = find(all(abs(FPlusMinusEpsilonLabel{j} - fmax) < tol,2));
                           
                           if ~isempty(t)
                               FPlusMinusEpsilonLabel{j}(t,4) = fmax(4) + perturbacion;
                               FPlusMinusEpsilonLabel{j}(floor(t/2),4) = FPlusMinusEpsilonLabel{j}(floor(t/2),4) + perturbacion ;
                               FPlusMinusEpsilon{j}(t,4) = fmax(4) + perturbacion;
                               FPlusMinusEpsilon{j}(floor(t/2),4) = FPlusMinusEpsilon{j}(floor(t/2),4) + perturbacion ;
                               FPlus{j}(floor(t/2),4) = fmax(4) + perturbacion;
                               FMinus{j}(floor(t/2),4)=FMinus{j}(floor(t/2),4) + perturbacion;
                               X{j}(floor(t/2),4) = X{j}(floor(t/2),4) + perturbacion;
                               if floor(t/2) ~= length(X{j}(:,1))
                                   FPlusMinusEpsilonLabel{j}(floor(t/2)+1,2) = FPlusMinusEpsilonLabel{j}(floor(t/2),4);
                                   FPlusMinusEpsilonLabel{j}(t+1,2) = FPlusMinusEpsilonLabel{j}(t,4) ;
                                   FPlusMinusEpsilon{j}(floor(t/2)+1,2) = FPlusMinusEpsilon{j}(floor(t/2),4);
                                   FPlusMinusEpsilon{j}(t+1,2) = FPlusMinusEpsilon{j}(t,4) ;
                                   FPlus{j}(floor(t/2)+1,2) = FPlus{j}(floor(t/2),4);
                                   FMinus{j}(floor(t/2)+1,2) = FMinus{j}(floor(t/2),4);
                                   X{j}(floor(t/2)+1,2) = X{j}(floor(t/2),4);
                               end
                           end
                       end
                       tq1 = find(Points(:,4) == i);
                       tq2 = find(Points(:,5) == i);
                       tq = sort([tq1',tq2']);
                       Points(tq,:) = [];
                       L = FPlusMinusEpsilonLabel{:,i};
                       X1Plus = L(1:length(L(:,1))/2,:);
                       X1Minus = L(length(L(:,1))/2+1:end,:);
                       X2 = cat(1,FPlusMinusEpsilon{:,1:i-1}); % ver como meter
                       X2_label = cat(1,FPlusMinusEpsilonLabel{:,1:i-1});
                       %X2 = X{:,i+1:m}
                       CPlus = Cruces(X1Plus, X2, X2_label);
                       CMinus = Cruces(X1Minus, X2, X2_label);
                       Points = cat(1,Points,CPlus);
                       Points = cat(1,Points,CMinus);
                       break
                   end
               end
               Points = [Points; CMinus(k,:) ];
               
           end
           if all(~isempty(CPlus) & aux == 1)
                Points = [Points; CPlus ];
           end
       end
    end
    %Points = sort(Points);
end
function [ n ] = Count2( int, FPlus, FMinus,event,j, M ) %ESTA MAL
%Count te devuelve el número de intersecciones en el intervalo int
    a = int(1);
    b = int(2);
    m = length(FPlus);
    n = 0;
    tol = 10^-12; %tolerancia pues uno es aprox
    x = event(j,1);
    if a ~= b
        L1 = event(:,6);
        L2 = event(:,7);
        index = find(abs(event(:,2) - a) < tol);
        nfmin1 = find(L1(1:index) == -1);
        nfmin2 = find(L2(1:index) == -1);
        nfplus1 = find(L1(1:index) == 1);
        nfplus2 = find(L2(1:index) == 1);           
        nmin1 = length(nfmin1);
        nmin2 = length(nfmin2);
        n1 = max(nmin1,nmin2);
        nmax1 = length(nfplus1);
        nmax2= length(nfplus2);
        n2 = max(nmax1,nmax2);
        n = n + (n1-n2);
        index2 = find(event(1:index,3) == 0);
        for i =1:length(index2) 
            if all(event(index2,6) == 1 & event(index2,7) == 1)
                n = n-1;
            elseif all(event(index2,6) == -1 & event(index2,7) == -1)
                n = n+1;
            end
        end
    else
        h = 0;
        for i = 1:m 
            long = length(FMinus{i}(:,1));
            for j = 1:long
                if all(FMinus{i}(j,1) <= x && FPlus{i}(j,3) >= x && h ~= i)
                    yFmin = nCruces(FMinus{i}(j,:),[x,-M,x,M]);
                    yFmax = nCruces(FPlus{i}(j,:),[x,-M,x,M]);
                    n2 = n;
                    n = n + sum((abs(a - yFmin(2)) < tol || a >= yFmin(2)) && (b <= yFmax(2)||(abs(b - yFmax(2)) < tol )));
                    if n2 ~= n
                        h = i;
                    end
                    
                end
            end
        end
   end
end
function [ Points, r ,q, position ] = Lf2(  F, events, points0, position, d )
%Lf2 te devuelve las intersecciones de las funciones +-epsilon con la sweep line en el punto x
    Points = points0;
    events = sortrows(events,2);
    k = find(ismember(position, events(:,4:7),'rows'));
    j = find(ismember(position, [fliplr(events(:,4:5)),fliplr(events(:,6:7))],'rows'));
    j = sort([k',j']); 
    l = length(j);
    len = length(events(:,1));
    r =[];
    q =[];
    tol = 10^-12;
    if l ~= 0
        k = j(1);
    end
    for m = 1:l
        i = j(m);
        if events(m,end) == 0  
            j1 = events(m,4);
            j2 = events(m,5);
            k1 = find(all([F{j1}(:,1) <= events(m,1) + tol F{j1}(:,5) == events(m,6)],2),1,'last');
            k2 = find(all([F{j2}(:,1) <= events(m,1) + tol F{j2}(:,5) == events(m,7)],2),1,'last');
            x = min(F{j1}(k1,3),F{j2}(k2,3));
            M = max(abs([F{j1}(k1,2),F{j2}(k2,2),F{j1}(k1,4),F{j2}(k2,4)]))+1;
            if x == F{j1}(k1,3)
                y1 = F{j1}(k1,4);
                y2 = Cruces([x -M x M 0 0 0 0],F{j2}(k2,1:4),F{j2}(k2,:));
                y2 = y2(2);
                if y1 > y2 
                    r = [r events(m,6)];
                    q = [q events(m,7)];
                    position(i,:) = [j2 j1 events(m,7) events(m,6)]; 
                    if i ~= 1
                        position(i-1,:) = [position(i-1,1) position(i,1) position(i-1,3) position(i,3)];
                    end
                    if i~= length(position(:,1))
                         position(i+1,:) = [position(i,2) position(i+1,2) position(i,4) position(i+1,4)];
                    end
                else
                    r = [r events(m,7)];
                    q = [q events(m,6)];
                    position(i,:) = [j1 j2 events(m,6) events(m,7)]; 
                    if i ~= 1
                        position(i-1,:) = [position(i-1,1) position(i,1) position(i-1,3) position(i,3)];
                    end
                    if i~= length(position(:,1))
                         position(i+1,:) = [position(i,2) position(i+1,2) position(i,4) position(i+1,4)];
                    end
                end
            else
                y1 = F{j2}(k2,4);
                y2 = Cruces([x -M x M 0 0 0 0],F{j1}(k1,1:4),F{j1}(k1,:));
                y2 = y2(2);
                if y1 > y2 
                    r = [r events(m,7)];
                    q = [q events(m,6)];
                    position(i,:) = [j1 j2 events(m,6) events(m,7)]; 
                    if i ~= 1
                        position(i-1,:) = [position(i-1,1) position(i,1) position(i-1,3) position(i,3)];
                    end
                    if i~= length(position(:,1))
                         position(i+1,:) = [position(i,2) position(i+1,2) position(i,4) position(i+1,4)];
                    end
                else
                    r = [r events(m,6)];
                    q = [q events(m,7)];
                    position(i,:) = [j2 j1 events(m,7) events(m,6)]; 
                    if i ~= 1
                        position(i-1,:) = [position(i-1,1) position(i,1) position(i-1,3) position(i,3)];
                    end
                    if i~= length(position(:,1))
                         position(i+1,:) = [position(i,2) position(i+1,2) position(i,4) position(i+1,4)];
                    end
                end
            end   
        else
            r = [r events(m,6)];
            q = [q events(m,7)];
        end
        if Points(k,3) == 0 
            Points(k,:) = [events(m,1:3) position(i,:)];
            if m ~= l
                k = j(m+1);
            end
        else
            Points(k,:) = [events(m,1:3) position(i,:)];
            Points(k+1,:)=[];
            if m ~= l
                k = j(m+1)-1;
            end
        end
    end
end
function y = nCruces(F, int)
    out = lineSegmentIntersect(F,int);
    P =[out.intMatrixX(:),out.intMatrixY(:)];
    y = [];
    m = length(P(:,1));
    for i = 1:m
        if ~isempty(find(out.coincAdjacencyMatrix,1))
           [J,list] = find(out.coincAdjacencyMatrix);
           len = length(list);
           for k = 1:len
               y = [y ; [int(J,1) int(J,2)]];
           end
        end
        if out.intAdjacencyMatrix(i) 
            y = [y ; P(i,:)];
        end
    end
end


    
