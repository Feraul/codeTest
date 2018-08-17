clear all
clc

%   Sample input:
        F = @(x) [...
        % equation 1
        0.1309 * 1 / sqrt(3) * ...
        (1 / sqrt(1 / 3 + x(1) ^ 2 + x(1) ^ 2) + ...
         1 / sqrt(1 / 3 + x(2) ^ 2 + x(4) ^ 2) + ...
         1 / sqrt(1 / 3 + x(3) ^ 2 + x(3) ^ 2) + ...
         1 / sqrt(1 / 3 + x(4) ^ 2 + x(2) ^ 2)) - 0.4352;
        % equation 2
        0.1309 * ...
        (x(1) / sqrt(1 / 3 + x(1) ^ 2 + x(1) ^ 2) + ...
         x(2) / sqrt(1 / 3 + x(2) ^ 2 + x(4) ^ 2) + ...
         x(3) / sqrt(1 / 3 + x(3) ^ 2 + x(3) ^ 2) + ...
         x(4) / sqrt(1 / 3 + x(4) ^ 2 + x(2) ^ 2)) - 0.1751;
        % equation 3
        0.1309 * ...
        (x(1) / sqrt(1 / 3 + x(1) ^ 2 + x(1) ^ 2) + ...
         x(4) / sqrt(1 / 3 + x(2) ^ 2 + x(4) ^ 2) + ...
         x(3) / sqrt(1 / 3 + x(3) ^ 2 + x(3) ^ 2) + ...
         x(2) / sqrt(1 / 3 + x(4) ^ 2 + x(2) ^ 2)) - 0.1751;
        % equation 4
        0.1309 * 1 / sqrt(3) * ...
        (x(1) / (1 / 3 + x(1) ^ 2 + x(1) ^ 2) + ...
         x(2) / (1 / 3 + x(2) ^ 2 + x(4) ^ 2) + ...
         x(3) / (1 / 3 + x(3) ^ 2 + x(3) ^ 2) + ...
         x(4) / (1 / 3 + x(4) ^ 2 + x(2) ^ 2)) - 0.1395; ];
        
        %x0 = [0.4330, 0.4330, 0.1443, 0.1443];
         x0 = [0.3, 0.3, 0.1, 0.1];

        epsilon = 1e-8;

        max_iter = 10;

 % Output:
        %x = [0.4151, 0.3863, 0.0976, 0.0818]

        %R = 1.0e-6 * [-0.0367, -0.2264, -0.2264, -0.2882] 

%% calculo do residuo Inicial
p_old=x0;
R0=F(x0);
R=R0;
%% inicializando dados para iteração Picard
step=0;
er=1;
while (1e-8<er || 1e-5==er) && (step<10)
    
    %% atualiza iterações
    step=step+1
    %########################################################
    p_extra(:,step)=p_old';
    if 1e-3<er
       j_v_approx = @(v)JV_APPROX(v, F, p_old);
        
       v = gmres(j_v_approx, R); % solve for Krylov vector
    
       p_new = p_old - v'; % updated solution
    else
        %p_rrrrr = extrapMPE(p_extra);

        p_rrrrr = extrapRRE(p_extra);
        p_new=p_rrrrr;
    end
    R = F(p_new); % new residual
        
    if (R0 ~= 0.0)
        er = abs(norm(R)/norm(R0))
        %er = R/norm(RHS_new)
    else
        er = 0.0; %exact
    end
        
    %% atualizar
    p_old=p_new;  % inversão sem pivotamento
    
end
p_old