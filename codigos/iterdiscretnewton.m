function [p,step,errorelativo,flowrate,flowresult]=iterdiscretnewton(M_old,RHS_old,tolnewton,kmap,...
    parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,nflagno,benchmark,...
    M_old1,RHS_old1,p_old1,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded)
% inicializando dados para iteração Picard
Rk= M_old1*p_old1-RHS_old1;
R0= M_old*p_old-RHS_old;
er=1;
step=0;

% calculo do jacobiano discreto
[J]=aproxmjacobian(R0,p_old1,p_old,nflagface,w,s,metodoP,parameter,...
    kmap,nflagno,benchmark,fonte,auxflag,gamma,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
%
while tolnewton<er
    step=step+1;
    % calculo inversa da Jaconiano
    
    pr=-J\Rk;
    % calculo da pressão
    p_new=p_old1+pr;
    %% plotagem no visit
     %S=ones(size(p_new,1),1);
     %postprocessor(p_new,S,step)
    % Interpolação das pressões na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);

    % Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflagface,nflagno...
        ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
  
    % calculo do residuo
    Rkk= M_new*p_new - RHS_new;
    
    % calculo do Jacobiano discreto
    [J]=aproxmjacobian(Rk,p_new,p_old1,nflagface,w,s,metodoP,parameter,kmap,...
        nflagno,benchmark,fonte,auxflag,gamma,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
    
    % calculo do erro
    A=logical(norm(R0) ~= 0.0);
    B=logical(norm(R0) == 0.0);
    er=A*abs(norm(Rkk)/norm(R0))+B*0
    errorelativo(step)=er;
    
    % atualizar
    p_old1=p_new;
    Rk=Rkk;
end
p=p_old1;
pinterp=pressureinterp(p,nflag,w,s,auxflag,metodoP,parameter,weightDMP);
if strcmp(metodoP,'nlfvDMPSY')
    
    % implementação do fluxo NLFV-DMP
    [flowrate,flowresult]=flowrateNLFVDMP(p, pinterp, parameter,nflag,kmap,gamma,weightDMP);
else
    
    % implementação do fluxo NLFV
    [flowrate,flowresult]=flowrateNLFV(p, pinterp, parameter);
end