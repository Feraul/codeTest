function [p,step,erro,flowrate,flowresult]=picardRRE(M_old,RHS_old,nitpicard,tolpicard,kmap,...
    parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,nflagno,benchmark,...
    weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface)

%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);

%% inicializando dados para iteração Picard
step=0;
er=1;
while (tolpicard<er || tolpicard==er) && (step<nitpicard)
    
    %% atualiza iterações
    step=step+1
    %########################################################
    p_extra(:,step)=p_old;
    if 1e-6<er
        p_new=M_old\RHS_old;  % inversão sem pivotamento
    else
        
        p_new = extrapRRE(p_extra);
        r=p_new(:)<0;
        x=find(r==1);
        if max(x)>0
            p_new(x)=0;
        end
        
    end
    
    %% plotagem no visit
    %S=ones(size(p_new,1),1);
    %postprocessor(p_new,S,step)
    p_max=max(p_new);
    p_min=min(p_new);
    %% Interpolação das pressões na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflagface,nflagno,w,s,auxflag,...
        metodoP,parameter,weightDMP);
    
    %% Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflagface,nflagno...
        ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,...
        wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
    
    %% calculo do erro
    R = norm(M_new*p_new - RHS_new);
    
    if (R0 ~= 0.0)
        er = abs(R/R0)
        %er = R/norm(RHS_new)
    else
        er = 0.0; %exact
    end
    erro(step)=er;
    
    %% atualizar
    M_old=full(M_new);
    RHS_old=RHS_new;
    p_old=p_new;  % inversão sem pivotamento
    
end
p=M_old\RHS_old;
pinterp=pressureinterp(p,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);

if strcmp(metodoP,'nlfvDMPSY')
    % implementação do fluxo NLFV-DMP
    [flowrate,flowresult]=flowrateNLFVDMP(p, pinterp, parameter,...
        nflagface,kmap,gamma,weightDMP,mobility);
else
    % implementação do fluxo NLFV
    [flowrate,flowresult]=flowrateNLFV(p, pinterp, parameter,mobility);
end
residuo=er;
niteracoes=step;

name = metodoP;
X = sprintf('Calculo o campo de pressão pelo método: %s ',name);
disp(X)

x=['Erro:',num2str(residuo)];
disp(x);
y=['Número de iterações:',num2str(niteracoes)];
disp(y);

end