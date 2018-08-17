function [J]=aproxmjacobian(Fk,p_new,p_old1,nflagface,w,s,metodoP,parameter,...
                kmap,nflagno,benchmark,fonte,auxflag,gamma,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded)
global elem
nelem=size(elem,1);
J=sparse(nelem,nelem);
pj=p_old1;

for ielem=1:nelem
    % atualiza a i-ésima coluna com xi+h
    pj(ielem)=p_new(ielem);
    
    % Interpolação das pressões nas faces
    [auxpinterp]=pressureinterp(pj,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
    
    % Calculo da matriz global
    [auxM,auxRHS]=globalmatrix(pj,auxpinterp,gamma,nflagface,nflagno...
        ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
    
    % Calculo do residou
    Fkk= auxM*pj - auxRHS;
    
    % Montagem da matriz Jacobiano por método Diferencia Finita
    J(1:nelem,ielem)=(Fkk(:)-Fk(:))./(p_new(ielem)-p_old1(ielem));
    
    % Atualiza o vetor "pj"
    pj=p_old1;
    
end

end