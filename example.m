F = {[0 1 5 1],[0 0.5 2 0.5; 2 0.5 5 0],[0 1.5 2 1.5; 2 1.5 5 2],[0 2 5 0],[0 0 5 2]};

% F = {matriz Nix4 de manera que cada fila es de la forma [x1 y1 x2 y2]

F2 = {[0 1 3 1; 3 1 5 3; 5 3 7 3],[0 4 7 4],[0 5 4.5 5; 4.5 5 6 4.5; 6 4.5 7 4.5]};

%ejemplo de Nadine

F_Nadine = {[0 1 5 1],[0 2.2 1 2.2 ; 1 2.2 5 1.2], [0 3.5 5 3.5],[0 4.2 3 4.2 ; 3 4.2 5 6]};

% aleatorio 
segment1 = [0 1 3 1];
segment2 = [0 1 3 8];

n = 10;
F_ran3 = cell(1,n);
s = 25;
x = 0;
for i = 1:n
    yz = [0 rand(1,1)*50];
    m = randi(s,1);
    for j = 1:m     
       xy = [0.1+rand(1,1)*5+yz(1) rand(1,1)*50];
       F_ran3{i}(j,:)=[yz xy];
       yz = xy;
       if j==m
            x2=F_ran3{i}(j,3);
       end
    end
   if x<x2
       for k = 1:i-1
          F_ran3{k}(end,3) = x2;
       end
       x = x2;
   else
       F_ran3{i}(m,3) = x;
   end
end

% F3 

F3 = {[0 1 3 1; 3 1 3.5 2; 3.5 2 4 1; 4 1 5 6],[0 2 2 0 ; 2 0 3 6; 3 6 5 2], [0 2.5 1.5 1;1.5 1 5 4],[0 3 2.5 3 ; 2.5 3 3.5 0; 3.5 0 5 1]};





