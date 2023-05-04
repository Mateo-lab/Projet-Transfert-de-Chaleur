format longG
% Définition du maillage
    % Nombre de noeuds en hauteur
imax=30;
    % Nombre de noeuds en longueur
jmax=120;

% Initialisation de la matrice T
T = zeros(imax,jmax);
% Température sur la surface extérieur de la vitre
Tpe = zeros(1,jmax);

% Dimension de la surface étudié en mètre
hauteur=0.005; %épaisseur de la vitre
largueur=0.014; %moitié de la distane entre deux source

dx=largueur/jmax;
% Attention, dy est une fonction de i, donc ce trouve dans la boucle

% Caractéristique du problème
    % Coefficient de conductivité du matériau
k=0.8;
    % Coéfficient de convection extérieur
he=11;
    % Température extérieur
Te=258;
    % Coéfficient de convection extérieur
hi=8;
    % Température extérieur
Ti=258;
    % Température minimum 
T0=278;
    % Source
q=15;

% Coefficient intermédiaire
    % Nombre de Biot supérieur
Biote=he*dx/k;
    % Nombre de Biot inférieur
Bioti=hi*dx/k;

% Paramètre pour la précision et son test
precision=0.000001;
testpr=1;
testT=1;

% Faire le tour de la matrice autant de fois de besoin pour avoir la
% précision demandé
iter=0;
while (testpr>precision)
  testpr=0;
    for i = 1:imax
        % Initialisation du dy avant (dy1) et après (dy2) maille 
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
        % dy1 représente la distance entre la maille et celle juste au-dessus (la précédante), pour la première maille, celle distance vaut 0 mais il n'y a pas de précédante, ce n'est donc pas un prblème
        % dy2 de la même manière est la distance entre le maille et celle au-dessous (la suivante), dy2 pour la dernière maille vaut 0, plus de maille suivante
        % L'évolution des dy est linéaire croissante jusqu'à la moitié puis décroissante de pente opposé après la moitié 
        

        for j = 1:jmax
            testT = T(i, j);

            if i > 1 && j > 1 && i < imax && j < jmax
                % Si la cellule n'est pas sur un bord ou un coin
                T(i,j) = (T(i-1,j)*dx^2*dy2 + T(i+1,j)*dx^2*dy1 + (T(i,j-1)+T(i,j+1))*((dy1+dy2)/2)*dy1*dy2) / (dx^2*dy1 + dx^2*dy2 + dy1^2*dy2 + dy2^2*dy1);
                %T(i,j) = (T(i-1,j) + T(i+1,j) + T(i,j-1) + T(i,j+1)) / 4;
            elseif i == 1 && j > 1 && j < jmax
                % Si la cellule est sur le bord supérieur
                %T(i,j) = (Biote*Te + T(i+1,j)*dx/dy2 + (T(i,j-1)+T(i,j+1))*dy2/(2*dx)) / (Biote + 2);

                T(i,j)=(((dy2/(2*dx))*T(i,j-1))+((dy2/(2*dx))*T(i,j+1))+((dx/dy2)*T(i+1,j))+(((he*dx)/k)*Te)) /(((he*dx)/k)+(dy2/dx)+(dx/dy2));
         
                %T(i,j) = (Ti + T(i+1,j) + T(i,j-1) + T(i,j+1)) / 4;
            elseif i == imax && j > 1 && j < jmax
                % Si la cellule est sur le bord inférieur
                %T(i,j) = (Bioti*Ti + T(i-1,j)*dx/dy1 + (T(i,j-1)+T(i,j+1))*dy1/(2*dx)) / (Bioti + 2);

                T(i,j)=(((dy1/(2*dx))*T(i,j-1))+((dy1/(2*dx))*T(i,j+1))+((dx/dy1)*T(i-1,j))+(((hi*dx)/k)*Ti)) /(((hi*dx)/k)+(dy1/dx)+(dx/dy1));

                %T(i,j) = (T(i-1,j) + Te + T(i,j-1) + T(i,j+1)) / 4;
            elseif j == 1 && i > 1 && i < imax
                % Si la cellule est sur le bord gauche
                T(i,j) = (T(i-1,j)*dx^2*dy2 + T(i+1,j)*dx^2*dy1 + T(i,j+1)*(dy1+dy2)*dy1*dy2) / (dx^2*dy1 + dx^2*dy2 + dy1^2*dy2 + dy2^2*dy1);
                %T(i,j) = (T(i-1,j) + T(i+1,j) + 2*T(i,j+1)) / 4;
            elseif j == jmax && i > 1 && i < imax
                % Si la cellule est sur le bord droit
                T(i,j) = (T(i-1,j)*dx^2*dy2 + T(i+1,j)*dx^2*dy1 + T(i,j-1)*(dy1+dy2)*dy1*dy2) / (dx^2*dy1 + dx^2*dy2 + dy1^2*dy2 + dy2^2*dy1);
                %T(i,j) = (T(i-1,j) + T(i+1,j) + 2*T(i,j-1)) / 4;
            elseif i == 1 && j == 1
                % Si la cellule est dans le coin supérieur gauche
                %T(i,j) = (Biote*Te + T(i+1,j)*dx/dy2 + T(i,j+1)*dy2/dx) / (Biote + 2);

                T(i,j)=(((dy2/dx)*T(i,j+1))+((dx/dy2)*T(i+1,j))+(((he*dx)/k)*Te)) /(((he*dx)/k)+(dy2/dx)+(dx/dy2));

                %T(i,j) = (Ti + T(i+1,j) + 2*T(i,j+1)) / 4;
            elseif i == 1 && j == jmax
                % Si la cellule est dans le coin supérieur droit
                %T(i,j) = (Biote*Te + T(i+1,j)*dx/dy2 + T(i,j-1)*dy2/dx) / (Biote + 2);

                T(i,j)=(((dy2/dx)*T(i,j-1))+((dx/dy2)*T(i+1,j))+(((he*dx)/k)*Te)) /(((he*dx)/k)+(dy2/dx)+(dx/dy2));

                %T(i,j) = (Ti + T(i+1,j) + 2*T(i,j-1)) / 4;
            elseif i == imax && j == 1
                % Si la cellule est dans le coin inférieur gauche
                %T(i,j) = (Bioti*Ti + T(i-1,j)*dx/dy1 + T(i,j+1)*dy1/dx + q/k) / (Bioti + 2);

                T(i,j)=(((dy1/dx)*T(i,j+1))+((dx/dy1)*T(i-1,j))+(((hi*dx)/k)*Ti+(q/k))) /(((hi*dx)/k)+(dy1/dx)+(dx/dy1));

                %T(i,j) = (T(i-1,j) + Te + 2*T(i,j+1) + q) / 4;
            elseif i == imax && j == jmax
                % Si la cellule est dans le coin inférieur droit
                %T(i,j) = (Bioti*Ti + T(i-1,j)*dx/dy1 + T(i,j-1)*dy1/dx) / (Bioti + 2);

                T(i,j)=(((dy1/dx)*T(i,j-1))+((dx/dy1)*T(i-1,j))+(((hi*dx)/k)*Ti)) /(((hi*dx)/k)+(dy1/dx)+(dx/dy1));

                %T(i,j) = (T(i-1,j) + Te + 2*T(i,j-1)) / 4;
            end

            if i == 1
                % On regarde les températures sur la surface éxtérieur de
                % la vitre
                Tpe(1,j)=T(i,j);
            end

            testpr=max(testpr,abs(T(i,j)-testT));
        end
    end
    iter=iter+1;
end

% Affichage de la température sur la vitre extérieur
plot(Tpe);

% Affichage des température en tout point

colormap(jet); %choix de la palette de couleur : "jet"
imagesc(T); %affichage sans la grille
%contourf(T); %affichage avec les courbes de températures
colorbar;
%axis equal; %mettre la même echelle pour les deux axes

% Vérification de la température minimum
if T(imax,jmax) >= T0
    title(sprintf('Itération = %d, Température minimum respecté',iter)); 
else
    title(sprintf('Itération = %d, Température minimum non respecté',iter)); 
end
% affichage de la matrice avec les temperatures exactes

%disp(T)
