% Escoamento oleo - agua (2-D) utilizando metodo IMPES
% Saturacao resolvida pelo metodo Upwind de primeira ordem
%% Este codigo somente roda MONOFASICO
clear all
clc
format short
global coord centelem elem esurn1 esurn2 nsurn1 nsurn2 bedge inedge ...
    normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens ...
    visc satlimit pormap bcflag courant totaltime filepath foldername;
%%========================================================================%

[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,normals,...
    esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,filepath,foldername,kmap,...
    wells] = preprocessor;
% bedge(217,4)=102;
% bedge(223,4)=102;
% bedge(229,4)=102;
% bedge(235,4)=102;
%% só em malha com furo
% x=bedge(217:240,1);
% y=bedge(217:240,2);
% bedge(217:240,1)=y;
% bedge(217:240,2)=x;
%% distorcao de malhas estruturadas
%[auxcoord]=distortedramd;
%% só em malha com furo
% x=bedge(129:144,1);
% y=bedge(129:144,2);
% bedge(129:144,1)=y;
% bedge(129:144,2)=x;
%% somente use nas malhas "Tipo1malha1", "Tipo1malha2", "Tipo1malha3" e "Tipo1malha4"

% Historial para malha "Tipo1malha0" não habilite nada, el ya viene
% ordenado en sentido antihorario tanto en elemento como no contorno
% para malha "Tipo1malha1" "Tipo1malha2 " "Tipo1malha3" e "Tipo1malha4"
% vamos habilitar na linha 803-806 do preprocessador, 
% e na linha, caso contrario vai dar erro no calculo do esurn1 esurn2 e no
% nsurn1 e no nsurn2

%   x=bedge(:,1);
%   y=bedge(:,2);
%   bedge(:,1)=y;
%   bedge(:,2)=x;
%   x1=elem(:,1);
%   x2=elem(:,3);
%   elem(:,1)=x2;
%   elem(:,3)=x1;   
%% use se deseja resolver o problema de Queiroz et al 2014
% unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
 % unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
%  x=bedge(71:78,1);
%  y=bedge(71:78,2);
%  bedge(71:78,1)=y;
%  bedge(71:78,2)=x;
%  bedge(71:78,4:5)=102; % benchmark23_3_GROSSO 18x18
%  bcflag(2,1)=102;
%  bcflag(2,2)=2;

%-----------------------------
%x=bedge(135:150,1);
%  y=bedge(135:150,2);
%  bedge(135:150,1)=y;
%  bedge(135:150,2)=x;
%  bedge(135:150,4:5)=102; % 36x36 cuidado
%  bcflag(2,1)=102;
%  bcflag(2,2)=2;

% unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
% x=bedge(289:320,1);
% y=bedge(289:320,2);
% bedge(289:320,1)=y;
% bedge(289:320,2)=x;
% bedge(289:320,4:5)=102; % benchmark23_3 72x72
% bcflag(2,1)=102;
% bcflag(2,2)=2;

% unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
%  x=bedge(577:640,1);
%  y=bedge(577:640,2);
%  bedge(577:640,1)=y;
%  bedge(577:640,2)=x;
%  bedge(577:640,4:5)=102; % benchmark23_3 144x144
%  bcflag(2,1)=102;
%  bcflag(2,2)=2;
%% Modificação Malha Kershaw
%bedge(:,4:5)=101;
%----------------------------
%% tratamento malha Hermeline
%bedge(:,4:5)=101;
% malha 16x16
 %x=bedge(16:24,1);
 %y=bedge(16:24,2);
 %bedge(16:24,1)=y;
 %bedge(16:24,2)=x;
% malha 32x32
% x=bedge(33:48,1);
% y=bedge(33:48,2);
% bedge(33:48,1)=y;
% bedge(33:48,2)=x;
% malha 64x64
% x=bedge(64:96,1);
% y=bedge(64:96,2);
% bedge(64:96,1)=y;
% bedge(64:96,2)=x;
% malha 128x128
%x=bedge(128:192,1);
%y=bedge(128:192,2);
%bedge(128:192,1)=y;
%bedge(128:192,2)=x;
%% calculo o flag do elemento que deseja 
%   a=6287;
%   b=445;
%   c=5740;
%   d=0;
%   [elemento]=searchelement(a,b,c,d)
%%
%tic
%% escolha o método de interpolação 
% 1-->LPEW1 Gao e Wu 2010
% 2-->LPEW2 Gao e Wu 2010
interptype=2;
% numérodo de Courant
CFL=courant;
%% exponente das permeabilidades relativas cuidado...! que varia com o benchmark
nw=2;
no=2;
%% Digite 
% monofasico ---> quando deseje rodar um problema de escoamento monofásico ou
% bifasico   ---> quando deseja rodar um problema de ecoamento bifásico ;
% quando rodar bifásico verifique linha 46-50 do Ks_Interp_LPEW2.m 
simu='monofasico';
%% escolha o tipo de erro discreto que deseja usar
% erromethod1 ---> erro utilizado por Gao e Wu 2010
% erromethod2 --->  ''     ''     por Lipnikov et al 2010
% erromethod3 --->  ''     ''     por Eigestad et al 2005
% erromethod4 --->  ''     ''     por Shen e Yuan 2015
erromethod='erromethod1';
%% Defina o tipo de solver de pressão
% tpfa      --> método Linear dos volumes finito TPFA
% mpfad     --> (MPFA-D)--> Alvaro
% lfvLPEW   --> método linear basedo no método não linear usando LPEW (MPFA-HD), ter cuidado linha 52 e 54 do preNLFV
% lfvHP     --> (MPFA-H)--> Allan 
% lfvEB     --> método completamente baseado na face (MPFA-BE).
% nlfvLPEW  --> (NLFV-PP)--> Abdiel
% nlfvDMPSY --> método não linear que preserva DMP baseado no artigo (Gao e Wu, 2013) e (Sheng e Yuan, 20...)

pmetodo='nlfvLPEW';
%% metodo de interação: iterpicard, iternewton, iterbroyden, itersecant,
% método de itereção proprio de métodos não lineares
%iterfreejacobian,iterdiscretnewton, JFNK
%iteration='iterdiscretnewton';
%iteration='iterbroyden';
%iteration='JFNK';
% iteration='iterpicard';
% iteration='fullpicard';
%iteration='RRE';
iteration='RRE';
%iteration='iterhybrid';
%% defina o ponto de interpolação
interpol='LPEW';
%interpol='LPEW';
%% Defina o tipo de solver de saturação
%1. escreve ['FOU'] se deseja executar FOU method
%2. escreve 'HOFV-E' se deseja executar método de alta por simple
%extrapolação
%3. escreve 'HOMFV' se deseja executar método de alta ordem com limitador
%de nó, na esqueça de definir upsilon (default é 0.2 como indica o artigo)
%2. escreve ['GOE'] se deseja executar o método GOE-free;
smetodo='FOU'; % cuidado no alterar !
order=1;       % cuidado no alterar !
%% ordem de Runge Kutta
tordem=1;      % cuidado no alterar !
% uso especial quando esta ativado o limitador de Woodfield, forneça um
% valor no intervalo [0,1]
upsilon=0.2;   % cuidado no alterar !
% Quando
% kapp=0--> método de Fromm
% kappa=1--> método de diferenças centradas com três pontos
% kappa=-1--> método de ponderação a montante de segunda ordem
% kappa=1/3--> método de ponderação a montante de terceira ordem
kappa=1/3;     % cuidado no alterar !
%% digite segundo o benchmark
% benchmark
% nikitin
% lamine
% durlofsky
% shuec
% buckley
benchmark='benchmar5_7';
%% nome do arquivo é unico para cada exemplo
% se utiliza quando se desea simular problemas bifásicos
namefile='Report_Production_Mesh_lamine_LFVHP.dat';
% escreve sobre o arquivo criado
fid3 = fopen(namefile,'w');
%% adequação dos flags
%nflag= calflag(pmetodo);
%% este flag só use quando o problema é Buckley-Leverett com fluxo imposto
% na face com flag 202, porém a saturação será imposto na face
auxflag=202;
%% adequação das permeabilidades, otros entes fisico-geometricos segundo o bechmark
[elem,kmap,normKmap,solanal,bedge,fonte,velanal]=adequapermeab(benchmark, kmap,elem,bedge);
% pre-aloca quando simulamos problemas monofásicos
 sat=ones(size(elem,1),1);
%postprocessor(solanal,sat,1);
% Condiçao inicial saturaçao
S_old = zeros(size(elem,1),1);

S_old(:)=satlimit(2);
% Condiçao de contorno da saturaçao

S_cont=1-satlimit(1);

% adequação dos poços
if max(max(wells))~=0
    for i=1:size(wells,1)
        if wells(i,2)==1
            S_old(wells(i,1))=S_cont;
        end
    end
else
    for ifacont=1:size(bedge,1)
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        rr=bcflag(r,2);
        
        if rr~=0 && bcflag(r,1)>200
            S_cont=1-satlimit(1);
        end
        
    end
end
%% pre-método-multidimensional
% não habilite, exclusivo para métodos multidimensionais
%[elemphant,elembedge,face_bedge,peso,nu]=premultidimensional(N);

%% pre-processador local
N=1;
[pointarmonic,parameter,gamma,p_old,tol,nit,er,nflagface,nflagno,...
    weightDMP,Hesq,Kde,Kn,Kt,Ded,auxface,calnormface]=preNLFV(kmap,N,pmetodo,benchmark,bedge);
nflag=nflagno;
% não habilite
%[aroundface]=aroundfacelement(F,pointarmonic);
%% IMPES
% inicializando as variáveis
cont = 1;
vpi_old = 0;
VPI(1) = 0;
cumulateoil(1) = 0;
oilrecovery(1) = 1;
watercut(1) = 0;
vpi = totaltime(2);
t_old = 0;
step=0;
time2=0;
%% calculo das mobilidades 
    [mobility] = mobilityface(S_old,nw,no,auxflag,S_cont,simu);
tic
%%  Calculo da Pressão Implicita pelos método lineares e não-lineares
    [p,errorelativo,flowrate,flowresult]=solverpressure(kmap,nflagface,nflagno,fonte, tol,...
        nit,p_old,mobility,gamma,wells,S_old,nw,no,parameter,...
        pmetodo,auxflag,interptype,Hesq, Kde, Kn, Kt, Ded,weightDMP,auxface,...
        benchmark,iteration,nflag,simu,calnormface);
toc

%step = 1000*time2;
%postprocessor(p,S_old,step)
%postprocessor(p,S_old,cont+1)
%% calculo do erro
% habilite se deseja calcular o erro entre a solução exata e a solução
% numérica
%[erropressure,errovelocity]=errorateconv(solanal, p, velanal,flowrate,erromethod)
%% calculo das pressões máximas e minimas
panalmax= max(p)
panalmin= min(p)


% if max(wells)==0
% plotBL(S_old,h)
% end
% S_old;

