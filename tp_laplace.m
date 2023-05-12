format longG
% ON POSE LE PROBLEME (DONNEES, CONSTANTES, TAILLE DU MAILLAGE)
    
        % Définition du maillage
imax=30;  % Nombre de noeuds en hauteur
jmax=120; % Nombre de noeuds en longueur


        % Initialisation des matrices utilisées
T = zeros(imax,jmax);  % Matrice des températures
Tpe = zeros(1,jmax);   % Vecteur des températures sur la surface extérieur de la vitre

        %Initialisation vecteurs positions
X = dx:dx:largueur; %contient les valeurs pour
Y = 1:imax;
        % Dimension de la surface étudié en mètre
hauteur=0.005; % Epaisseur de la vitre
largueur=0.014; % Moitié de la distane entre deux source
dx=largueur/jmax; % Ecart entre chaque maille sur l'axe X
% dy est variable donc ce trouve dans la boucle

        % Données du problème
k=0.8;    % Coefficient de conductivité du matériau
he=11;    % Coéfficient de convection extérieur
Te=258;    % Température extérieur
hi=8;    % Coéfficient de convection extérieur
Ti=258;    % Température extérieur
T0=278;    % Température minimum 
q=410/33;     % Source q = Puissance de la source / la longueur de fils


        % Coefficient intermédiaire
Biote=he*dx/k;    % Nombre de Biot supérieur
Bioti=hi*dx/k;    % Nombre de Biot inférieur


        % Paramètre pour la précision et son test
precision=0.000001;
testpr=1; %Initialisation de la variable de test de la précision 
testT=1;  %Initialisation de la variable de test de la Temparute

%PARTIE ALGORITHME DE CALCUL
        
iter=0; % Initialisation de la variable qui permet de donner le nombre d'itération que la boucle a éffectué
    
while (testpr>precision)
    % La boucle 'While' va permettre de faire le tour de la matrice autant de fois que besoin pour avoir la
    % précision demandé
  testpr=0;
  LargeurTampon = 0.005; %permet de créer la matrice Y
    for i = 1:imax
    
        %Partie du pas variable
        % dy1 représente la distance entre la maille et celle juste au-dessus (la précédante), pour la première maille, celle distance vaut 0 mais il n'y a pas de précédante, ce n'est donc pas un prblème
        % dy2 de la même manière est la distance entre le maille et celle au-dessous (la suivante), dy2 pour la dernière maille vaut 0, plus de maille suivante
        % L'évolution des dy est linéaire croissante jusqu'à la moitié puis décroissante de pente opposé après la moitié 
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
        LargeurTampon = LargeurTampon - dy1; % Va permettre de contruire le vecteur position Y
        Y(i) = LargeurTampon;  % On remarque que le vecteur Y est contruit à l'envers, il commence à 0.005 et non à 0
                               % Il sera inversé lors de son utilisation.
        if i>=imax
            Y(i)=0;
        end
        
        %Partie Calcul
        
        for j = 1:jmax
            testT = T(i, j);

            if i > 1 && j > 1 && i < imax && j < jmax       % Si la cellule n'est pas sur un bord ou un coin
                
                T(i,j) = (T(i-1,j)*dx^2*dy2 + T(i+1,j)*dx^2*dy1 + (T(i,j-1)+T(i,j+1))*((dy1+dy2)/2)*dy1*dy2) / (dx^2*dy1 + dx^2*dy2 + dy1^2*dy2 + dy2^2*dy1);

            elseif i == 1 && j > 1 && j < jmax      % Si la cellule est sur le bord supérieur
            
                T(i,j)=(((dy2/(2*dx))*T(i,j-1))+((dy2/(2*dx))*T(i,j+1))+((dx/dy2)*T(i+1,j))+(((he*dx)/k)*Te)) /(((he*dx)/k)+(dy2/dx)+(dx/dy2));
         
            elseif i == imax && j > 1 && j < jmax    % Si la cellule est sur le bord inférieur
              

                T(i,j)=(((dy1/(2*dx))*T(i,j-1))+((dy1/(2*dx))*T(i,j+1))+((dx/dy1)*T(i-1,j))+(((hi*dx)/k)*Ti)) /(((hi*dx)/k)+(dy1/dx)+(dx/dy1));

            elseif j == 1 && i > 1 && i < imax      % Si la cellule est sur le bord gauche
               
                T(i,j) = (T(i-1,j)*dx^2*dy2 + T(i+1,j)*dx^2*dy1 + T(i,j+1)*(dy1+dy2)*dy1*dy2) / (dx^2*dy1 + dx^2*dy2 + dy1^2*dy2 + dy2^2*dy1);

            elseif j == jmax && i > 1 && i < imax       % Si la cellule est sur le bord droit
                
                T(i,j) = (T(i-1,j)*dx^2*dy2 + T(i+1,j)*dx^2*dy1 + T(i,j-1)*(dy1+dy2)*dy1*dy2) / (dx^2*dy1 + dx^2*dy2 + dy1^2*dy2 + dy2^2*dy1);

            elseif i == 1 && j == 1     % Si la cellule est dans le coin supérieur gauche
                
                T(i,j)=(((dy2/dx)*T(i,j+1))+((dx/dy2)*T(i+1,j))+(((he*dx)/k)*Te)) /(((he*dx)/k)+(dy2/dx)+(dx/dy2));

            elseif i == 1 && j == jmax      % Si la cellule est dans le coin supérieur droit
               

                T(i,j)=(((dy2/dx)*T(i,j-1))+((dx/dy2)*T(i+1,j))+(((he*dx)/k)*Te)) /(((he*dx)/k)+(dy2/dx)+(dx/dy2));

            elseif i == imax && j == 1      % Si la cellule est dans le coin inférieur gauche

                T(i,j)=(((dy1/dx)*T(i,j+1))+((dx/dy1)*T(i-1,j))+(((hi*dx)/k)*Ti+(q/k))) /(((hi*dx)/k)+(dy1/dx)+(dx/dy1));

            elseif i == imax && j == jmax       % Si la cellule est dans le coin inférieur droit
                
                T(i,j)=(((dy1/dx)*T(i,j-1))+((dx/dy1)*T(i-1,j))+(((hi*dx)/k)*Ti)) /(((hi*dx)/k)+(dy1/dx)+(dx/dy1));

            end

            if i == 1       % On regarde les températures sur la surface éxtérieur de la vitre
                
                Tpe(1,j)=T(i,j);

            end

            testpr=max(testpr,abs(T(i,j)-testT));

        end
    end
    
    iter=iter+1;

end

      % Calcul des gradients de la température en x et y
[qx, qy] = gradient(-k * flipud(T));

      % Vérification de la température minimum sur l'extérieur de la vitre
if T(1,jmax) >= T0 %T(1,jmax) represente la T° en haut à droite, car 1ère ligne et jème colonne
    title(sprintf('Itération = %d, Température minimum respecté',iter)); 
else
    title(sprintf('Itération = %d, Température minimum non respecté',iter)); 
end

%PARTIE AFFICHAGE 

      % Affichage de la température sur la vitre extérieur
plot(Tpe);

        % Affichage des température en tout point
subplot(2, 1, 1); % Mettre plusieurs sous graphiques dans la même fenêtre, ici le 1
colormap(jet); % Choix de la palette de couleur : "jet"
%imagesc(X, flip(Y), T); % Affichage sans la grille, 'flip()' permet d'inverser un vecteur
pcolor(X, Y, T); % Permet d'afficher me maillage 

        % Affichage flux de chaleur
subplot(2, 1, 2); % Mettre plusieurs sous graphiques dans la même fenêtre, ici le 2
contour(X, Y, T,10); % Tracer les courbes de niveaux de la matrice T, le '8' correspond au nombre de courbes desirées, 'flip()' permet d'inverser un vecteur
hold on;
quiver(X, flip(Y), qx, qy, 2); % Tracer les vecteurs de flux de chaleur représentés par les composantes qx et qy
title('Flux de chaleur et isotherme');

