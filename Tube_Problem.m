function  Tube_Problem(F, epsilon, p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here    
    [ Answer,FPlusMinusEpsilonLabel,~,~,Events, boolean ] = Decision_Problem(F, epsilon, p);
    if Answer
        m1 = length(F);
        %p = ceil(p/m); si queremos que p sea porcentaje
        m = length(Events(:,1));
        l2 = m1*2 -1;%#intervalos
        intersections = cell(1,m);
        Matrix  = [];
        tol = 10^-12;
        for i = 1:length(FPlusMinusEpsilonLabel)
            Matrix = [Matrix; FPlusMinusEpsilonLabel{1,i}];
        end
        L = Matrix(:,4);
        M = max([max(Matrix(:,2)),max(L)]) + 1; %NÚMERO SUFICIENTEMENTE GRANDE
        Y = cell(m, l2);%donde estan las intersecciones en cada intervalo al que se puede llegar en continuo
        c = distinguishable_colors(m1+1);%colors
        e1 = [];
        for i = 1:m1 % Si quieres ver visualemente las funciones
            
            plot([F{1,i}(:,1),F{1,i}(:,3)],[F{1,i}(:,2),F{1,i}(:,4)],'x-','Color', c(i,:),'LineWidth',3);
            hold on; 
        end
        i = 1;
        for k = 1:m
            t = find(abs(Events(:,1) - Events(i,1)) < tol);
            e = Events(i:t(end),:);
            e1 = [e1, Events(i,1)];
            i = t(end);
            j2 = 1; %por si se repiten varios
            intersections{k} = Lf(Matrix, e, M);
            aux = 0;   
            for j = 1:l2
                if intersections{k}(j2,3) ~= 0 || aux == 1 %si este punto no es el cruce             
                    int = [intersections{k}(j2,2), intersections{k}(j2+1,2)]; 
                    j2 = j2+1;
                    aux = 0;
                else %si en este evento j es un cruce
                    int = [intersections{k}(j2,2), intersections{k}(j2,2)];
                    aux = 1;
                end
                %n(j) = Count2(int, FPlus, FMinus,intersections{k},j2,M ); 
                if boolean(k,j)
                        Y{k,j} =  int; 
                end
            end
            if i == m 
                Y2 = [];
                for x = k:-1:1
                    z = find(boolean(x,:));                
                    if x == k
                        Y2 = [Y2 (Y{x,z(1)}(1)+Y{x,z(1)}(2))/2 ];
                        z2 = z(1);
                    else
                        ind = 1;
                        while ~(z2 == z(ind) || z2-1 == z(ind) || z2+1 == z(ind))
                            ind = ind + 1;
                        end
                        Y2 = [Y2 (Y{x,z(ind)}(1)+Y{x,z(ind)}(2))/2];
                        z2 = z(ind);
                    end
                end
                Y2 = flip(Y2);
                %e1 = uniquetol(Events(:,1), tol);
                xq = [e1(1):0.1:e1(end)];
                hold on
                vq1 = interp1(e1,Y2,xq);
                plot(e1,Y2,'ro',xq,vq1,':.','LineWidth',0.5);
                hold on
                plot(xq, vq1, 'x--','Color', c(end,:),'LineWidth',0.5);
                hold on
                patch([e1, fliplr(e1)],[Y2-epsilon,fliplr(Y2+epsilon)],c(end,:));
                alpha(0.25);
                hold off
                break;
            end
            i = i+1;
        end
    else
        error('it is not possible to construct the tube');       
    end
end
