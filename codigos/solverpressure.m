function [pressure,errorelativo,flowrate,flowresult]=solverpressure(kmap,nflagface,nflagno,fonte,...
    tol, nit,p_old,mobility,gamma,wells,S_old,nw,no,parameter,metodoP,...
    auxflag,interptype,Hesq, Kde, Kn, Kt, Ded,weightDMP,auxface,benchmark,...
    iteration,nflag,simu,calnormface)

errorelativo=0;
w=0;
s=0;
if interptype==1
    % calculo dos pesos que correspondem ao LPEW1
    [w,s]=Pre_LPEW_1(kmap,mobility,V,S_old,nw,no,N);
elseif interptype==2
    % calculo dos pesos que correspondem ao LPEW2
    [w,s]=Pre_LPEW_2(kmap);
end
W=w';
% incializando variaveis

flowrate=0;
flowresult=0;


% interpola��o nos n�s ou faces
[pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);

% calculo da matriz globlal inicial
[M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
    parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
    mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);

switch metodoP
    
    case {'nlfvLPEW', 'nlfvDMPSY','nlfvDMPV1'}
        
        if strcmp(iteration,'AA')
            
            [pressure,step,errorelativo,flowrate,flowresult]=picardAA(M_old,RHS_old,nit,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface);
        elseif strcmp(iteration,'RRE')
            
            [pressure,step,errorelativo,flowrate,flowresult]=picardRRE(M_old,RHS_old,nit,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface);
        elseif strcmp(iteration,'MPE')
            [pressure,step,errorelativo,flowrate,flowresult]=picardMPE(M_old,RHS_old,nit,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface);
        elseif strcmp(iteration,'fullpicard')
            [pressure,step,errorelativo,flowrate,flowresult]=fullpicard(M_old,RHS_old,nit,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,accelerator,calnormface);
            
        elseif strcmp(iteration,'iterbroyden')
            
            p_old1=M_old\RHS_old;
            
            % interpola��o nas faces
            [pinterp1]=pressureinterp(p_old1,nflagface,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
            % residuo inicial
            R0_old=M_old1*p_old1-RHS_old1;
            
            % solver de press�o pelo m�todo Broyden
            [pressure,step,errorelativo,flowrate,flowresult]=iterbroyden(M_old,RHS_old,p_old,tol,kmap,parameter,...
                metodoP,auxflag,w,s,nflagface,fonte,gamma,nflagno,benchmark,R0_old,p_old1,...
                weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
            
        elseif strcmp(iteration, 'iterdiscretnewton')
            
            p_old1=M_old\RHS_old;
            % interpola��o nos n�s ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
            % resolvedor de press�o pelo m�todo de Newton-Discreto
            [pressure,step,errorelativo,flowrate,flowresult]=iterdiscretnewton(M_old,RHS_old,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,M_old1,RHS_old1,p_old1,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
            
        elseif strcmp(iteration, 'iterhybrid')
            
            p_old1=M_old\RHS_old;
            
            % interpola��o nos n�s ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
            % solver pressure pelo m�todo hybrido
            [pressure,step,errorelativo,flowrate,flowresult]=iterhybrid(M_old1,RHS_old1,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,p_old1,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
        elseif strcmp(iteration, 'JFNK')
            
            
            p_old1=M_old\RHS_old;
            % calculo do residuo
            R0=norm(M_old*p_old-RHS_old);
            
            % interpola��o nos n�s ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,...
                auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
            % calculo da press�o
            [pressure,step,errorelativo,flowrate,flowresult]= JFNK1(tol,kmap,parameter,metodoP,auxflag,w,s,nflagface,fonte,gamma,...
                nflagno,benchmark,M_old1,RHS_old1,p_old1,R0,weightDMP,...
                auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
        end
        
    case {'lfvHP','lfvLPEW','mpfad','tpfa'}
        
        pressure=M_old\RHS_old;
        %%%%%% Indica de a matriz � diagonal dominante ou n�o
        %flag1 = domdiag( full(M_old), 'strict' );
        %flag1 = domdiag( full(M_old));
        % flag1---> 1 a matriz � estrictamente diagonal dominante
        % flag1---> 0 a matriz n�o diagonal dominante
        %%%%%%%
        if strcmp(metodoP, 'lfvHP')
            % interpola��o nos n�s ou faces
            [pinterp]=pressureinterp(pressure,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            % calculo das vaz�es
            [flowrate,flowresult]=flowratelfvHP(parameter,weightDMP,mobility,pinterp,pressure);
        elseif strcmp(metodoP, 'lfvLPEW')
            [pinterp]=pressureinterp(pressure,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            % calculo das vaz�es
            [flowrate,flowresult]=flowratelfvLPEW(parameter,weightDMP,mobility,pinterp,pressure);
        elseif  strcmp(metodoP, 'tpfa')
            [flowrate, flowresult]=flowrateTPFA(pressure,Kde,Kn,Hesq,nflag,mobility);
        else
            % calculo das vaz�es
            [flowrate,flowresult]=calflowrateMPFAD(pressure,w,s,Kde,Ded,Kn,Kt,Hesq,nflag,auxflag,mobility);
        end
        residuo=0;
        niteracoes=0;
        
        name = metodoP;
        X = sprintf('Calculo da press�o pelo m�todo: %s ',name);
        disp(X)
        
        x=['Residuo:',num2str(residuo)];
        disp(x);
        y=['N�mero de itera��es:',num2str(niteracoes)];
        disp(y);
        
end

end