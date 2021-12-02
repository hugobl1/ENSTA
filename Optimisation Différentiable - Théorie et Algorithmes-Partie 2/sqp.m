function [x,lme,lmi,info] = sqp(simul,x,lme,lmi,options)
    info = struct;
    info.niter = 0;     %nombre d'itérations effectuées
    n = length(x);
    nn = n/2;
    nb = nn + 1;
    epsilon = 10^-3;
    global norme;
    %Si lm est vide à l'initialisation, on choisit lm tel que le
    %solveur n'itère pas si x est stationnaire
    if(isempty([lme;lmi]))
        [~,ce,ci,g,ae,ai,~,~] = feval(simul,4,x,lme,lmi);
        me=length(ce);
        mi=length(ci);
        a=[ae;ai];
        A=[zeros(me,me) zeros(me,mi) ; zeros(mi,me) eye(mi,mi)];
        Aeq=[zeros(me,me) zeros(me,mi) ; zeros(mi,me) diag(ci)];
        size(transpose([zeros(+me,1) ;ci]))
        lm = quadprog(2*a*transpose(a),2*a*g,A,zeros(mi+me,1),Aeq,zeros(mi+me,1));
        lme = lm(1:nb);
        if(~isempty(ai))
            lmi = lm(nb+1:end);
        end
        
    %Sinon on vérifie que lm a les bonnes dimensions
    elseif(length([lme;lmi])~=nb)
        info.status = 1;
        return
    end
    
    %Puis on fixe la valeur finale de m et p
    m = length([lme;lmi]);
    p = length(lmi)/nn;
    % Calcul des valeurs utiles
    [~,ce,ci,g,ae,ai,hl] = feval(simul,5,x,lme,lmi);
    
    if(options.verb==1)
        fprintf("------------------------------------------------------------------------------------------\n");
        fprintf("\n");
        fprintf("iter       |gl|           |ce|         (ci,lmi)         |x|     |lm|        Powell       cond(M)\n");
    end
    % Boucle principale
    while true
        % Test d'arrêt du maximum d'itérations
        if(info.niter >= options.maxit)
            info.status = 2;
            break;
        end
        

        if ((options.deriv == 1) &&(info.niter==0))
            M = speye(n);
        end
        if (options.deriv == 2 )

            [l,d] = cholmod(hl,1.e-5,1.e+5);
            M = l*diag(d)*l';
        end
        
        
        % Test d'arrêt de la convergence
        if( (norm(g+[ae;ai]'*[lme;lmi],inf) <= options.tol(1)) && (norm(ce,inf) <= options.tol(2)) && (norm(min(lmi,-ci),inf)) <= options.tol(3))
            info.status = 0;
            break;
        end
        
        % Calcul des valeurs suivantes de x et lme et lmi
        lb=[];
        ub=[];
        x0=[];
        optionsquad = optimset('Display', 'off');
        [d, ~ , ~, ~, lambda] = quadprog(M,g,ai,-ci,ae,-ce,lb,ub,x0,optionsquad);
        
        if ((options.deriv == 1)||(options.deriv == 2))
            x_k = x;
            difflm=[lambda.eqlin;lambda.ineqlin]-lm;
        end
        lm = [lambda.eqlin;lambda.ineqlin];
        x = x + d;
        lme = lambda.eqlin;
        if(~isempty(ai))
            lmi = lambda.ineqlin;
        end
        % On actualise les valeurs utiles
        [~,ce,ci,g,ae,ai,hl] = feval(simul,5,x,lme,lmi);
        %on actualise M dans le cas de la méthode de Quasi-Newton
        if (options.deriv == 1)
        %on calcule le vecteur variation du gradient du lagrangien
            [~,~,~,g_k,ae_k,ai_k,~,indic] = feval(simul,4,x_k,lme,lmi);
            if indic == 1
                error('Problème avec la fonction chs');
            end
            if (~isempty(lmi)) 
                gamma_l = g+transpose(ae)*lme+transpose(ai)*lmi-g_k-transpose(ae_k)*lme-transpose(ai_k)*lmi;
            else
                gamma_l = g+transpose(ae)*lme-g_k-transpose(ae_k)*lme;
            end
            %on actualise theta en fonction de la correction de Powell
            if (transpose(gamma_l)*d < 0.2*transpose(d)*M*d)
               theta = 0.8*(transpose(d)*M*d)/(transpose(d)*M*d-transpose(gamma_l)*d);
            else
               theta = 1;
            end
            %on calcule gamma
            gamma = (1-theta)*M*d+theta*gamma_l;

            if (info.niter == 0)
                eta = (transpose(gamma)*gamma)/(transpose(gamma)*d);
                %on actualise M
                M = eta*M;
            end
            M = M -((M*d*transpose(d)*M)/(transpose(d)*M*d))+(gamma*transpose(gamma)/(transpose(gamma)*d));
            
        end
        info.niter = info.niter + 1;
        %%Si l'on veut calculer la vitesse de convergence:
        norme=[norme norm([d ;difflm] , inf)];
             
        %%Si l'on veut afficher les informations de chaque itération
        if(options.verb==1)
             fprintf('    %d    %.2e    %.2e      %.2e        %.1e    %.2e    %.1e    %.1e\n',info.niter, norm(g+[ae' ai']*[lme;lmi]), norm(ce,Inf), norm(-ci.*lmi),norm(x,Inf),norm(lm,Inf),theta,cond(M));
        end
        
    end
        if((options.verb==1))
            fprintf("\n");
            fprintf("------------------------------------------------------------------------------------------\n");
            fprintf("\n");
        end
    
end