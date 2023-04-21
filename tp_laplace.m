% Définition du nombre de mailles
    % Nombre de noeuds en hauteur
imax=30;
    % Nombre de noeuds en longueur
jmax=120;

% Initialisation de la matrice T
T = zeros(imax,jmax);

hauteur=3;
largueur=12;

dx=largueur/jmax;
% Attention, dy est une fonction de i, donc ce trouve dans la boucle

lamda=150;
h=400;
Text=900;

% Paramètre pour la précision et son test
precision=0.001;
testpr=1;
testT=1;

% Faire le tour de la matrice autant de fois de besoin pour avoir la
% précision demandé
iter=0;
while (testpr>precision)
  testpr=0;
    for i = 1:imax
        % Initialisation du dy avant (dy1) et après (dy2) cellule 
        if (i+1) < (1+imax/2)
            dy1=(hauteur*4/(imax^2))*(i-1);
            dy2=(hauteur*4/(imax^2))*(i);
        elseif i <= (1+imax/2) && (i+1) >= (1+imax/2)
            dy1=(hauteur*4/(imax^2))*(i-1);
            dy2=(hauteur*4/(imax^2))*(imax-i);
        elseif i > (1+imax/2)
            dy1=(hauteur*4/(imax^2))*(imax-i+1);
            dy2=(hauteur*4/(imax^2))*(imax-i);
        end 
        

        for j = 1:jmax
            testT = T(i, j);

            if i > 1 && j > 1 && i < imax && j < jmax
                % Si la cellule n'est pas sur un bord ou un coin
                T(i,j) = (T(i-1,j) + T(i+1,j) + T(i,j-1) + T(i,j+1)) / 4;
            elseif i == 1 && j > 1 && j < jmax
                % Si la cellule est sur le bord supérieur
                T(i,j) = 100;
            elseif i == imax && j > 1 && j < jmax
                % Si la cellule est sur le bord inférieur
                T(i,j) = 200;
            elseif j == 1 && i > 1 && i < imax
                % Si la cellule est sur le bord gauche
                T(i,j) = (T(i-1,j) + T(i+1,j) + T(i,j+1)) / 3;
            elseif j == jmax && i > 1 && i < imax
                % Si la cellule est sur le bord droit
                T(i,j) = (T(i-1,j) + T(i+1,j) + T(i,j-1)) / 3;
            elseif i == 1 && j == 1
                % Si la cellule est dans le coin supérieur gauche
                T(i,j) = (T(i+1,j) + T(i,j+1)) / 2;
            elseif i == 1 && j == jmax
                % Si la cellule est dans le coin supérieur droit
                T(i,j) = (T(i+1,j) + T(i,j-1)) / 2;
            elseif i == imax && j == 1
                % Si la cellule est dans le coin inférieur gauche
                T(i,j) = (T(i-1,j) + T(i,j+1)) / 2;
            elseif i == imax && j == jmax
                % Si la cellule est dans le coin inférieur droit
                T(i,j) = (T(i-1,j) + T(i,j-1)) / 2;
            end

            testpr=max(testpr,abs(T(i,j)-testT));
        end
    end
    iter=iter+1;
end

T= reshape(T, imax, jmax);
contourf(T);
colorbar;
title(sprintf('Itération = %d',iter)); 