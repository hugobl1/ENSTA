function [e,ce,ci,g,ae,ai,hl,indic] = chs(indic,xy,lme,lmi)
    global L;
    global A;
    global B;
    global R;
    global S;
    n = length(xy);
    nn = n/2;
    nb = nn + 1;
    p = length(S);
    x = [0;xy(1:nn);A];     % Coordonnées de la totalité de la chaîne
    y = [0;xy(nn+1:n);B];
    hl = sparse(2*nn,2*nn);
    if(indic == 1)
        % Tracé de la chaîne
        p1 = plot(x,y);
        title('Chaîne de travail');
        xlabel('x');
        ylabel('y');
        xlim([min(x) max(x)]);
        ylim([min(y) max(y)]);
        hold on;
        
        % Tracé du plancher
        if(p>0)
            x_plancher_list = linspace(-20,20,1000);   % Liste des abcisses du plancher
            y_plancher_list = [];                   % Liste des ordonnées du plancher
            for x = x_plancher_list
               y_plancher_list = [y_plancher_list,max(x*S+R)];
            end
            p2 = area(x_plancher_list,y_plancher_list,'basevalue',-100,'Facecolor',[0.7 0.7 0.7]);
        end
        hold off;
        uistack(p1,'top');
    end
    
    if(indic >= 2)
        % Calcul de e
        Mp = spdiags([ones(nb,1) ones(nb,1)],0:1,nb,nb+1);
        e = 0.5*sum(L.*(Mp*y));
        
        % Calcul de ce
        Mn = spdiags([-ones(nb,1) ones(nb,1)],0:1,nb,nb+1);
        ce = (Mn*x).^2 + (Mn*y).^2 - L.^2;
        
        % Calcul de ci
        if(p==0)
            ci = [];
        else
            Kn = kron(eye(p),ones(nn,1));
            Rn = Kn*R;
            Sn = Kn*S;
            Xn = repmat(xy(1:nn),p,1);
            Yn = repmat(xy(nn+1:n),p,1);
            ci = Rn + Sn.*Xn - Yn;
        end
    end
        
    if(indic >= 4)
        % Calcul de g
        L1 = L(1:nn);
        L2 = L(2:nn+1);
        g = [zeros(nn,1);0.5*(L1+L2)];
        
        % Calcul de ae
        aex = Mn.*repmat(Mn*x,1,nb+1);
        aex = aex(1:nb,2:nb);
        aey = Mn.*repmat(Mn*y,1,nb+1);
        aey = aey(1:nb,2:nb);
        ae = 2*[aex,aey];
        
        % Calcul de ai
        if(p==0)
            ai=[];
        else
            I1 = repmat(eye(nn),p,1);
            aix = bsxfun(@times,I1,Sn);
            aiy = repmat(-eye(nn),p,1);
            ai = [aix,aiy];
        end
    end
        
    if(indic >= 5)
        % Calcul de hl
        u = [1 -1;-1 1];
        b = sparse(nn,nn);
        hl = sparse(2*nn,2*nn);
        for i=1:nb
            if(i==1)
                b(1,1) = 2;
            elseif(i==nb)
                b(nn,nn) = 2;
            else
                b(i-1:i,i-1:i)=2*u;
            end
            hci = [[b;zeros(nn)],[zeros(nn);b]];
            hl = hl + lme(i)*hci;
            b = sparse(nn,nn);
        end
    end
    indic=0;
end

