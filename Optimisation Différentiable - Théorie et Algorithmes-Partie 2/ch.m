% Initialisation des donnees
global A;
global B;
global L;
global R;
global S;
global norme; norme=[];
% Tests de l'optimiseur
test_ref = '5a';
options = struct;
options.tol = [10^-6;10^-6;10^-6];
options.maxit = 1000;
options.deriv=1;
options.verb=1;

warning('off','all')
warning

%Test de l'optimiseur pour le cas test 4a
if strcmp(test_ref,'3c')
    A = 0;
    B = -1;
    L = [2;1];
    R = [];
    S = [];
    xy = [0.3;0.3];
            chs(1,xy);
    [x,lme,lmi,info] = sqp('chs',xy,[],[],options);
    chs(1,x,lme,lmi);
    info.status;
    info.niter;
    [e,ce,ci,g,ae,ai,hl,indic] = chs(5,x,lme,lmi);
    cprime=full([ae; ai]);
    size(cprime)
    rank(cprime)
    %[e,c,g,a,hl] = chs(5,x,lme,lmi);
    %norm(F(c,g,a,lm),inf);
    %norm(c,inf);
    %eig(hl);
end
%Test de l'optimiseur pour le cas test 4a
if strcmp(test_ref,'4a')
    A = 1;
    B = -1;
    L = [0.7;0.5;0.3;0.2;0.5];
    R = [];
    S = [];
    xy = [0.2;0.4;0.6;0.8;1;1.5;1.5;1.3];
            chs(1,xy);
    [x,lme,lmi,info] = sqp('chs',xy,[],[],options);
    chs(1,x,lme,lmi);
    info.status;
    info.niter;
    [e,ce,ci,g,ae,ai,hl,indic] = chs(5,x,lme,lmi);
    cprime=full([ae; ai]);
    size(cprime)
    rank(cprime)
    
    norme1=norme(1:(end-1));
    norme2=norme(2:end);
    Quot_list=norme2./norme1;
    figure
    plot(Quot_list);
    title("Etude la converge de la suite")
    xlabel("Nombre d'itérations")
    ylabel("||s_{k+1}||/||s_{k}|| ")
    %[e,c,g,a,hl] = chs(5,x,lme,lmi);
    %norm(F(c,g,a,lm),inf);
    %norm(c,inf);
    %eig(hl);
end

%Test de l'optimiseur pour le cas test 4b
if strcmp(test_ref,'4b')
    A = 1;
    B = 0;
    L = [0.2;0.2;0.2;0.3;0.3;0.5;0.2;0.2;0.3;0.1];
    R = [-0.25];
    S = [-0.5];
    xy = [0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;-0.5;-0.9;-1.2;-1.4;-1.5;-1.4;-1.2;-0.9;-0.5];

    [x,lme,lmi,info] = sqp('chs',xy,[],[],options);
    chs(1,x,lme,lmi);
    info.status;
    info.niter;
    
    norme1=norme(1:(end-1));
    norme2=norme(2:end);
    Quot_list=norme2./norme1;
    figure
    plot(Quot_list);
    title("Etude la converge de la suite")
    xlabel("Nombre d'itérations")
    ylabel("||s_{k+1}||/||s_{k}|| ")
    %[e,c,g,a,hl] = chs(5,x,lme,lmi);
    %norm(F(c,g,a,lm),inf);
    %norm(c,inf);
    %eig(hl)
end

%Test de l'optimiseur pour le cas test 4c
if strcmp(test_ref,'4c')
    A = 1;
    B = 0;
    L = [0.2;0.2;0.2;0.3;0.3;0.5;0.2;0.2;0.3;0.1];
    R = [-0.25;-0.5];
    S = [-0.5;0];
    xy = [0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;-0.5;-0.9;-1.2;-1.4;-1.5;-1.4;-1.2;-0.9;-0.5];

    [x,lme,lmi,info] = sqp('chs',xy,[],[],options);
    chs(1,x,lme,lmi);
    info.status;
    info.niter;
    
    norme1=norme(1:(end-1));
    norme2=norme(2:end);
    Quot_list=norme2./norme1;
    figure
    plot(Quot_list);
    title("Etude la converge de la suite")
    xlabel("Nombre d'itérations")
    ylabel("||s_{k+1}||/||s_{k}|| ")

    %norm(F(c,g,a,lm),inf);
    %norm(c,inf);
    %eig(hl)
end
%Test de l'optimiseur pour le cas test 5a
if strcmp(test_ref,'5a')
    A = 0;
    B = 0;
    L = [0.5;0.3;0.4;1.2;0.3;0.3];
    R = [-1];
    S = [-0.1];
    xy = [0.2;0.5;0.8;1.0;1.2;-0.4;-0.6;-0.4;-0.2;0];
    chs(1,xy);
    [x,lme,lmi,info] = sqp('chs',xy,[],[],options);
    chs(1,x,lme,lmi);
    info.status;
    info.niter
    [e,ce,ci,g,ae,ai,hl,indic] = chs(5,x,lme,lmi);
    cprime=full([ae; ai]);
    size(cprime);
    rank(cprime);
    eig(hl)
    min(eig(hl))
    norme1=norme(1:(end-1));
    norme2=norme(2:end);
    Quot_list=norme2./norme1;
    figure
    plot(Quot_list);
    title("Etude la converge de la suite")
    xlabel("Nombre d'itérations")
    ylabel("||s_{k+1}||/||s_{k}|| ")
    
    %[e,c,g,a,hl] = chs(5,x,lme,lmi);
    %norm(F(c,g,a,lm),inf);
    %norm(c,inf);
    %eig(hl);
end
%Test de l'optimiseur pour le cas test 5b
if strcmp(test_ref,'5b')
    A = 0;
    B = -4;
    L = [3;2.5;2.5];
    R = [-6;-10];
    S = [-2;100];
    xy = [-2;0;1;-2];
    chs(1,xy);
    [x,lme,lmi,info] = sqp('chs',xy,[],[],options);
    chs(1,x,lme,lmi);
    info.status;
    info.niter
    [e,ce,ci,g,ae,ai,hl,indic] = chs(5,x,lme,lmi);
    cprime=full([ae; ai]);

    
    norme1=norme(1:(end-1));
    norme2=norme(2:end);
    Quot_list=norme2./norme1;
    figure
    plot(Quot_list);
    title("Etude la converge de la suite")
    xlabel("Nombre d'itérations")
    ylabel("||s_{k+1}||/||s_{k}||")

    %[e,c,g,a,hl] = chs(5,x,lme,lmi);
    %norm(F(c,g,a,lm),inf);
    %norm(c,inf);
    %eig(hl);
end
%Test de l'optimiseur pour le cas test 5c
if strcmp(test_ref,'5c')
    A = 0;
    B = 0;
    L = [0.1;0.2;0.3;0.4;0.5;0.4;0.3;0.1];
    R = [-1;-0.2;-1.0];
    S = [-7;0.0;7];
    xy = [0.1; 0.2; 0.1; -0.1; -0.3;-0.2;-0.2; -0.3;-0.1;0.2; 0.4; 0.3;0.2;0];%A trouver
    chs(1,xy);
    [x,lme,lmi,info] = sqp('chs',xy,[],[],options);
    chs(1,x,lme,lmi);
    info.status;
    info.niter;
    [e,ce,ci,g,ae,ai,hl,indic] = chs(5,x,lme,lmi);
    cprime=full([ae; ai]);
    size(cprime);
    rank(cprime);
    [e,c,g,a,hl] = chs(5,x,lme,lmi);
    
        
    norme1=norme(1:(end-1));
    norme2=norme(2:end);
    Quot_list=norme2./norme1;
    figure
    plot(Quot_list);
    title("Etude la converge de la suite")
    xlabel("Nombre d'itérations")
    ylabel("||s_{k+1}||/||s_{k}|| ")
    %norm(F(c,g,a,lm),inf);
    %norm(c,inf);
    %eig(hl);
end
%Test de l'optimiseur pour le cas test 5d
if strcmp(test_ref,'5d')
    A = 4;
    B = 14;
    L = [5;4;3;2;2;3;4;5];
    R = [-30;-12;-3;0;0;-3;-12;-30];
    S = [-10;-6;-3;-1;1;3;6;10];
    xy = [1; 2; 1; 1; 3;2;2;3; 1;2; 4; 3;2;0];
    chs(1,xy);
    [x,lme,lmi,info] = sqp('chs',xy,[],[],options);
    chs(1,x,lme,lmi);
    info.status;
    info.niter;
    [e,ce,ci,g,ae,ai,hl,indic] = chs(5,x,lme,lmi);
    cprime=full([ae; ai]);
    size(cprime);
    rank(cprime);
    [e,c,g,a,hl] = chs(5,x,lme,lmi);
    
    
    norme1=norme(1:(end-1));
    norme2=norme(2:end);
    Quot_list=norme2./norme1;
    figure
    plot(Quot_list);
    title("Etude la converge de la suite")
    xlabel("Nombre d'itérations")
    ylabel("||s_{k+1}||/||s_{k}|| ")

    %norm(F(c,g,a,lm),inf);
    %norm(c,inf);
    %eig(hl);
end


