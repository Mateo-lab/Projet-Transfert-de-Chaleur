%test banche thomas
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
    % Source q = Puissance de la source / la longueur de fils
q=410/33; 

% Coefficient intermédiaire
y=zeros(jmax);
b=zeros(jmax);

% Paramètre pour la précision et son test
precision=0.0001;
testpr=1;
testT=1;

% Faire le tour de la matrice autant de fois de besoin pour avoir la
% précision demandé
iter=0;
while (testpr>precision)
  testpr=0;

  % Première ligne, on fixe d'abord pour le premier point puis une boucle
  % pour les autres
  i=1;
  j=1;
  
  dy2=(hauteur*4/(imax^2))*(i);

  % Les expression ne sont pas exactement les memes pour les bords
  y(j+1)=-(dy2/dx)/(((he*dx)/k)+(dy2/dx)+(dx/dy2));
  b(j+1)=(((dx/dy2)*T(i+1,j))+(((he*dx)/k)*Te))/(((he*dx)/k)+(dy2/dx)+(dx/dy2));
  for j = 2:jmax-1
      y(j+1)=-(dy2/2*dx)/((dy2/2*dx)*y(j)+((he*dx)/k)+(dy2/dx)+(dx/dy2));
      b(j+1)=(((dx/dy2)*T(i+1,j))+(((he*dx)/k)*Te)-(dy2/2*dx)*b(j))/((dy2/2*dx)*y(j)+((he*dx)/k)+(dy2/dx)+(dx/dy2));
  end
  testT=T(i,jmax);
  T(i,jmax)=(((dx/dy2)*T(i+1,jmax))+(((he*dx)/k)*Te)-(dy2/dx)*b(jmax))/((dy2/dx)*y(jmax)+((he*dx)/k)+(dy2/dx)+(dx/dy2));
  testpr=max(testpr,abs(T(i,jmax)-testT));
  for j = jmax:-1:2
      testT=T(i,j-1);
      T(i,j-1)=y(j)*T(i,j)+b(j);
      testpr=max(testpr,abs(T(i,j-1)-testT));
  end


   
  % Toutes les ligne entre la première et la dernière
  for i = 2:imax-1
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
      
      j=1;
      y(j+1)=-((dy1+dy2)/dx)/((dx/dy1) + (dx/dy2) + ((dy1+dy2)/dx));      
      b(j+1)=(T(i-1,j)*(dx/dy1) + T(i+1,j)*(dx/dy2))/((dx/dy1) + (dx/dy2) + ((dy1+dy2)/dx));
      for j = 2:jmax-1      
          y(j+1)=-(((dy1+dy2)/2*dx))/(((dy1+dy2)/2*dx)*y(j) + (dx/dy1) + (dx/dy2) + ((dy1+dy2)/dx));      
          b(j+1)=(T(i-1,j)*(dx/dy1) + T(i+1,j)*(dx/dy2)-((dy1+dy2)/2*dx)*b(j))/(((dy1+dy2)/2*dx)*y(j) + (dx/dy1) + (dx/dy2) + ((dy1+dy2)/dx));
      end
      testT=T(i,jmax);
      T(i,jmax) = (T(i-1,jmax)*(dx/dy1) + T(i+1,jmax)*(dx/dy2)-((dy1+dy2)/dx)*b(jmax))/(((dy1+dy2)/dx)*y(jmax) + (dx/dy1) + (dx/dy2) + ((dy1+dy2)/dx));
      testpr=max(testpr,abs(T(i,jmax)-testT));
      for j = jmax:-1:2
          testT=T(i,j-1);
          T(i,j-1)=y(j)*T(i,j)+b(j);
          testpr=max(testpr,abs(T(i,j-1)-testT));
      end
  end



  % La dernière ligne
  i=imax;
  j=1;
  
  dy1=(hauteur*4/(imax^2))*(imax-i+1);  
  
  y(j+1)=-(dy1/dx)/(((hi*dx)/k)+(dy1/dx)+(dx/dy1));
  b(j+1)=(((dx/dy1)*T(i-1,j))+(((hi*dx)/k)*Ti+(q/k)))/(((hi*dx)/k)+(dy1/dx)+(dx/dy1));
  for j = 2:jmax-1
      y(j+1)=-(dy1/2*dx)/((dy1/2*dx)*y(j)+((hi*dx)/k)+(dy1/dx)+(dx/dy1));
      b(j+1)=(((dx/dy1)*T(i-1,j))+(((hi*dx)/k)*Ti)-(dy1/2*dx)*b(j))/((dy1/2*dx)*y(j)+((hi*dx)/k)+(dy1/dx)+(dx/dy1));
  end
  testT=T(i,jmax);
  T(i,jmax)=(((dx/dy1)*T(i-1,jmax))+(((hi*dx)/k)*Ti)-(dy1/dx)*b(jmax))/((dy1/dx)*y(jmax)+((hi*dx)/k)+(dy1/dx)+(dx/dy1));
  testpr=max(testpr,abs(T(i,jmax)-testT));
  for j = jmax:-1:2
      testT=T(i,j-1);
      T(i,j-1)=y(j)*T(i,j)+b(j);
      testpr=max(testpr,abs(T(i,j-1)-testT));
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
if T(1,jmax) >= T0 %T(1,jmax) represente la T° en haut à droite, car 1ère ligne et jème colonne
    title(sprintf('Itération = %d, Température minimum respecté',iter)); 
else
    title(sprintf('Itération = %d, Température minimum non respecté',iter)); 
end
% affichage de la matrice avec les temperatures exactes

%disp(T)