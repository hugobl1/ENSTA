%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Projet OPT201-OPT202                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pilotage du solveur

% Variables pilotant le solveur

indic=1;

% indic = 1 : chs fait un tracé de la chaîne (fonction plot);
% indic = 2 : chs calcule e et c;
% indic = 4 : chs calcule e, c, g et a;
% indic = 5 : chs calcule hl.

%% Définition des variables du problème
warning('off','all')
warning
%Coordonnée du point de fixation extrème
global A; A = 1;

global B; B=-1;

%Longueur des barres et taille du vecteur associé

global L; 

%Cas test 2:
L=[0.7;0.5;0.3;0.2;0.5];

%Quelles actions souhaite on effectuer:
global test; test="2d";
deriv=false;

%Test cas simple:
%L=[0.5;1;0.5];
%Test cas simple non surjectif:
%L=[0.5;0.5];



% Coordonnées des noeuds non fixés
if strcmp(test,"2a")
  %Cas test 2a:
  xy=[0.2;0.4;0.6;0.8;-1.0;-1.5;-1.5;-1.3];
end

if strcmp(test,"2b")
  %Cas test 2b:
  xy=[0.2;0.4;0.6;0.8;1.0;1.5;1.5;1.3];
end

if strcmp(test,"2c")
  %Cas test 2c:
  xy=[0.2;0.4;0.6;0.8;-1.0;-1.5;1.5;-1.3];
end

if strcmp(test,"2d")
  %Cas test 2d:
  xy=[0.2;0.4;0.6;0.8;1.0;-1.2;1.5;-1.3];
end

if strcmp(test,"3a")
  %Cas test 3a:
  A=1;
  B=0;
  L=[0.6;0.6];
  xy=[0.5;0.4];
end

if strcmp(test,"3b")
  %Cas test 3b:
  A=1;
  B=0;
  L=[2;1];
  xy=[0.5;0.3];
end

if strcmp(test,"3c")
  %Cas test 3c:
  A=0;
  B=-1;
  L=[2;1];
  xy=[0.3;0.3];
end

if strcmp(test,"non_surjectif")
  xy=[0.5;0];
end 

%Test cas simple (TP1):
%xy=[0.5;0.5;0;-1];
global Nb; Nb=length(L);
global N; N=length(xy);


################################################################################
%%Appel au solveur et ? l'optimiseur

%On affiche la cha?ne de d?part:
[~,~,~,~,~,~]=chs(1,xy);


%On calcule a et g afin de pouvoir d?terminer le multiplicateur lagrangien
lm=ones(Nb,1);
[~,~,g,a,~,~]=chs(4,xy,lm);

lm=-a'\g;

%g+a'*lm

%On d?finit la tol?rance et le nombre maximale d'it?rations

options.tol=[0.001 0.001];
options.maxit=1000;
options.rl=0;
options.verb=2;
%On lance l'optimiseur:

[x,lm,info] = sqp(@chs,xy,lm,options);

%On trace la cha?ne obtenue:

if strcmp(test,"3a")
  %Cas test 3a:
  xy=[0.5;sqrt(0.11)];%%first
  [~,~,~,~,~,~]=chs(1,xy);
  [e,c,~,a,~,~]=chs(4,xy,lm);
   c
   e
   xy=[0.5;-sqrt(0.11)];%%second
  [~,~,~,~,~,~]=chs(1,xy);
  [e,c,~,a,~,~]=chs(4,xy,lm);
  c
  e

elseif strcmp(test,"3b")
  %Cas test 3b:
  xy=[2;0];
  [~,~,~,~,~,~]=chs(1,xy);
  [e,c,~,a,~,~]=chs(4,xy,lm);
  c
  e
elseif strcmp(test,"3c")
  %Cas test 3b:
  xy=[0;-2];
  [~,~,~,~,~,~]=chs(1,xy);
  [e,c,~,a,~,~]=chs(4,xy,lm);
  c
  e
end 

[~,~,~,~,~,~]=chs(1,x);
[e,c,~,a,~,~]=chs(4,x,lm);
[~,~,~,~,L,~]=chs(5,x,lm);


K=null(a);
[vecteurs_propres,valeurs_propres]=eig((K')*L*K);
valeurs_propres



################################################################################

%Exactitude des dérivées:

if(deriv)


deriv_e=g;%dérivée calculée de l'énergie

deriv_c=a;%dérivée calculée des contraintes

%1/Exactitude dérivée de e
t_i=(sqrt(eps))*max(1,abs(xy));


DF=[];
for i=1:N
  df=0;
  XY_plus=xy;
  XY_plus(i)+=t_i(i);
  [e,c,g,a,hl,indic]=chs(2,XY_plus,lm);
  df+=e;
  XY_moins=xy;
  XY_moins(i)-=t_i(i);
  [e,c,g,a,hl,indic]=chs(2,XY_moins,lm);
  df-=e;
  df=df/(2*t_i(i));
  DF=[DF; df];
  
end
relativ=(deriv_e-DF)./deriv_e;
fprintf('\n %g\n',(deriv_e-DF)./deriv_e);

%2/Exactitude dérivée de c

t_i=(sqrt(eps))*max(1,abs(xy));


DF_bis=zeros(Nb,N);
for i=1: Nb
for j=1:N
  df=0;
  XY_plus=xy;
  XY_plus(j)+=t_i(j);
  [e,c,g,a,hl,indic]=chs(2,XY_plus,lm);
  df+=c(i);
  XY_moins=xy;
  XY_moins(j)-=t_i(j);
  [e,c,g,a,hl,indic]=chs(2,XY_moins,lm);
  df-=c(i);
  df=df/(2*t_i(j));
  DF_bis(i,j)=df;
  end
end

relativ_bis=(deriv_c-DF_bis)./deriv_c;
fprintf('\n %g\n',relativ_bis);

[e,c,g,a,hl,indic]=chs(4,xy,lm);
Z=a';

end

