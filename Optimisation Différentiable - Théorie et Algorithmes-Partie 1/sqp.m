function [x,lm,info] = sqp (simul,x,lm,options)
  %% Méthode sans recherche linéaire
  if (options.rl==1)
  ##On récupère le hession associé au point de départ
  
  [e,c,g,a,hl,indic]=simul(5,x,lm);
  
  ##On définit la variable Lk pour stocké la valeur courante du hessien
  
  Lk=hl;
  
  ##On récupère le gradient des contraintes
  
  [e,c,g,a,hl,indic]=simul(4,x,lm);
  
  ##On définit la variable Ak pour stocké la valeur courante du gradient des contraintes
  
  Ak=a;
  
  info.niter=0;
  
  k=size(x);
  
  n=size(lm);
  
  %%Si l'on veut calculer la vitesse de convergence:
  norm=[];
  #On vérifie la bonne taille des variables
  
  if (k!=(2*n-2))
    info.status=1;
  else
    
    ##On vérifie que les conditions d'entrée de boucle sont vérifiées
    
    while (((max(abs(g+a'*lm)) > options.tol(1))||(max(abs(c)) > options.tol(2)))&&(info.niter<options.maxit))
      info.niter+=1;
      X=[Lk Ak'; Ak zeros(size(Ak)(1))]\[-g;-c];
      x=x+X(1:k); %%On met à jour la valeur de x
      lm=X((k+1):(k+n));  %%On met à jour la valeur du coefficient de lagrange
      
      ##On récupère le hessien associé au point courant##
      
      [e,c,g,a,hl,indic]=simul(5,x,lm);
      Lk=hl;
      
      ##On récupère le gradient des contraintes associé au point courant##
      
      [e,c,g,a,hl,indic]=simul(4,x,lm);
      Ak=a;
      
      %[~,~,~,~,~,~]=chs(1,x);
      
      %pause(1);
      
      %%Si l'on veut calculer la vitesse de convergence:
      norm=[norm max(max(abs(g+a'*lm)),max(abs(c)))];
      
    end
    
    %%On affiche la vitesse de convergence: 
    figure
    plot(log(norm));
    title("Etude la converge de la suite")
    xlabel("Nombre d'itérations")
    ylabel("Norme de Fk (échelle logarithmique)")
    
    if ((max(abs(g+a'*lm)) <= options.tol(1))&&(max(abs(c)) <= options.tol(2)))
      info.status=0;
    elseif (info.niter>=options.maxit)
      info.status=2;
    end
  end
  elseif (options.rl==0)  ##Méthode avec recherche linéaire
##On récupère le hession associé au point de départ
  
  [e,c,g,a,hl,indic]=simul(5,x,lm);
  
  ##On définit la variable Lk pour stocké la valeur courante du hessien
  
  Lk=hl;
  
  ##On récupère le gradient des contraintes
  
  [e,c,g,a,hl,indic]=simul(4,x,lm);
  
  ##On définit la variable Ak pour stocké la valeur courante du gradient des contraintes
  
  Ak=a;
  
  info.niter=0;
  
  k=size(x);
  
  n=size(lm);
  
  %%Si l'on veut calculer la vitesse de convergence:
  norme=[];
  norme2=[];
  
  #On vérifie la bonne taille des variables
  
  if (k!=(2*n-2))
    info.status=1;
  else
    
    ##On vérifie que les conditions d'entrée de boucle sont vérifiées:
    %%Si l'on veut afficher des informations(2 modes possibles):
    if(options.verb==1)
fprintf("------------------------------------------------------------------------------------------\n");
      fprintf("\n");
      fprintf("iter       |gl|           |ce|         |x|       |lm|        alpha         phi\n");
    end
    if(options.verb==2)
      appelsimul=2;%%On a déjà fait 2 appels ici
    end
    
    %%On initialise les variables intermédiaires du problème 
    
    omega=10^(-4);
   
    Fk=[+g+Ak'*lm;+c];
    
    X=[Lk Ak'; Ak zeros(size(Ak)(1))]\(-Fk);
    
    phi=(0.5)*((norm(Fk))^2);
    
    ##On met à jour les données tant que la tolérance désirée n'est pas atteinte
    while (((max(abs(g+a'*lm)) > options.tol(1))||(max(abs(c)) > options.tol(2)))&&(info.niter<options.maxit))
      info.niter+=1;
      if(options.verb==2)
        fprintf("------------------------------------------------------------------------------------------\n");
        fprintf("\n");
        fprintf("iter  %d, simul  %d,  phi  %.5e,  pente %.5e",info.niter,appelsimul,phi,-2*phi);
        fprintf("\n");
        fprintf("\n");
        fprintf("recherche lineaire d’Armijo: |d| = %.2e",norm(X(1:k)));
        fprintf("\n");    
        fprintf("    alpha     phip-phi     DF(phi)\n");  
      end
      ##Déterminons alpha_k
      ##On initialise r à 0
      r=0;
      ##On prend comme premier pas alpha_k=1
      alpha_k=2^(-r);
      
      [eplus,cplus,gplus,aplus,~,indicplus]=simul(4,x+alpha_k*X(1:k),lm+alpha_k*X((k+1):(k+n)));
      
      Fkplus=[+gplus+aplus'*(lm+alpha_k*X((k+1):(k+n)));+cplus];
      if(options.verb==2)
        appelsimul+=1;
      end
      ##On calcule phi(z_k+alpha_k*p_k)
      phiplus=(0.5)*((norm(Fkplus))^2);
      
      ##On incrémente tant que la règle d'Armijo n'est pas vérifié
        if(options.verb==2)
          fprintf(" %.4e  %.5e  %.5e\n",alpha_k,phiplus-phi,(phiplus-phi)/alpha_k);  
        end
      while(phiplus>phi*(1-2*omega*alpha_k))
        r=r+1;

        alpha_k=2^(-r);
        
        [eplus,cplus,gplus,aplus,~,indicplus]=simul(4,x+alpha_k*X(1:k),lm+alpha_k*X((k+1):(k+n)));

        Fkplus=[+gplus+aplus'*(lm+alpha_k*X((k+1):(k+n)));+cplus];
        
        phiplus=(0.5)*((norm(Fkplus))^2);
        if(options.verb==2)
          appelsimul+=1;
          fprintf(" %.4e  %.5e  %.5e\n",alpha_k,phiplus-phi,(phiplus-phi)/alpha_k);  
        end
      end

      x=x+alpha_k*X(1:k); %%On met à jour la valeur de x
      lm=lm+alpha_k*X((k+1):(k+n));  %%On met à jour la valeur du coefficient de lagrange
      
      ##On récupère le hessien associé au point courant##
      
      [~,~,~,~,hl,indic]=simul(5,x,lm);
      if(options.verb==2)
        appelsimul+=1;
      end
      Lk=hl;
      
      ##On met à jour les différentes variables intermédiaires du problème     
      e=eplus;
      
      c=cplus;
      
      g=gplus;

      a=aplus;
      Ak=aplus;
      
      Fk=[+g+Ak'*lm;+c];
    
      X=[Lk Ak'; Ak zeros(size(Ak)(1))]\(-Fk);
    
      phi=(0.5)*((norm(Fk))^2);
      
      %[~,~,~,~,~,~]=chs(1,x);
      
      %pause(1);
      
      %%Si l'on veut calculer la vitesse de convergence:
      norme=[norme max(max(abs(g+a'*lm)),max(abs(c)))];
      global test;
      if strcmp(test,"3a")
        
      norme2=[norme2 norm(x-[0.5;sqrt(0.11)])];
    elseif strcmp(test,"3b")
      norme2=[norme2 norm(x-[2;0])];
    elseif strcmp(test,"3c")
      norme2=[norme2 norm(x-[0;-2])];
      end
      
      %%Si l'on veut afficher les informations de chaque itération
      if(options.verb==1)
        fprintf('    %d    %.4e    %.4e    %.1e    %.1e    %.3e    %.5e\n',info.niter, norm(g+Ak'*lm,Inf), norm(c,Inf),norm(x,Inf),norm(lm,Inf),alpha_k,phi);
      end
      if(options.verb==2)
        fprintf("\n");
        fprintf("  |gl| = %.3e  |ce| = %.3e \n", norm(g+Ak'*lm,Inf), norm(c,Inf));  
      end
 
      
    end
    if((options.verb==1)||(options.verb==2))
      fprintf("\n");
      fprintf("------------------------------------------------------------------------------------------\n");
      fprintf("\n");
    end

    %%On affiche la vitesse de convergence: 
    if (strcmp(test,"3a")||strcmp(test,"3b")||strcmp(test,"3c"))
          figure
          plot(norme2);
          title("Etude la converge de la suite")
          xlabel("Nombre d'itérations")
          ylabel("||x_k-x*||")
    else
          figure
          plot(log(norme));
          title("Etude la converge de la suite")
          xlabel("Nombre d'itérations")
          ylabel("Norme de Fk (échelle logarithmique)")
    end
    
    if ((max(abs(g+a'*lm)) <= options.tol(1))&&(max(abs(c)) <= options.tol(2)))
      info.status=0;
    elseif (info.niter>=options.maxit)
      info.status=2;
    end

  end  
end
end

  




