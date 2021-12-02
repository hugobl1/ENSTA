function [e,c,g,a,hl,indic]=chs(indic,xy,lm)
  global A;
  global B;
  global L; 
  global Nb;
  global N;
  % On liste les abcisses et les coordonnées des noeuds 
  % dans deux vecteurs distincts ce qui sera utile pour 
  % afficher la chaîne
  %On initialise chacune des variables à renvoyer afin d'éviter 
  % tout problème de non existence à la compilation
  e=0;
  c= zeros(Nb,1);
  g=zeros(2*(Nb-1),1);
  a=zeros(Nb,2*(Nb-1));
  hl=zeros(2*(Nb-1),2*(Nb-1));
  %On vérifie le bon dimensionnement des données en entrée, 
  %dans le cas où il n'est pas bon on renvoie indic=1
  if (mod(N,2)==1)
    indic=1;
  else
    x=[0; xy(1:N/2); A];
    y=[0; xy((N/2+1):N);B];
    
  if indic==1
    plot(x,y);
    indic=0;
    
  elseif indic==2
    y_i=y(2:Nb+1); 
    y_j=y(1:Nb);
    x_i=x(2:Nb+1);
    x_j=x(1:Nb);
    e=(L')*(y_i+y_j)/2;
    c=(x_i-x_j).^2+(y_i-y_j).^2-L.^2;
    indic=0;
    
  elseif indic==4
    y_i=y(2:Nb+1);
    y_j=y(1:Nb);
    x_i=x(2:Nb+1);
    x_j=x(1:Nb);
    L_i=L(1:(Nb-1));
    L_j=L(2:Nb);
    e=(L')*(y_i+y_j)/2;
    c=(x_i-x_j).^2+(y_i-y_j).^2-L.^2;
    g=[zeros(Nb-1,1) ; (L_i+L_j)/2];
    a= [spdiags([2*(x_j-x_i)(2:Nb) 2*(x_i-x_j)(1:Nb-1)], -1:0, Nb,Nb-1),...
    spdiags([2*(y_j-y_i)(2:Nb) 2*(y_i-y_j)(1:Nb-1)], -1:0, Nb,Nb-1)]; 
    indic=0;
    
  elseif indic==5
    zero=zeros(Nb-1,Nb-1);
    aux=zero;
    aux(1,1)+=lm(1)*2;
    hl =[aux, zero;
         zero, aux;];
    for i=2:(Nb-1)
      aux=zero;
      aux(i,i)=2;
      aux(i-1,i-1)=2;
      aux(i-1,i-1)=2;
      aux(i-1,i)=-2;
      aux(i,i-1)=-2;
      hl+=lm(i)*[aux, zero;
         zero, aux;];
      end
    aux=zero;
    aux(Nb-1,Nb-1)=lm(Nb)*2;
    hl +=[aux, zero;
         zero, aux;];
  end
 end
 end
 
 