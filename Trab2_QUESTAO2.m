%% RE-Inicializa A MEMORIA
clear
clc
close all


%METODOS DE DIRECAO

%Univariante        - 1
%1.87840	 20.23754	 12.81299 	 -4.22e+00 	 -6.40e+00 	 3.9308e-04 Barreira (OK)
%1.87745	 20.25115	 12.80924 	 -6.34e-02 	 1.40e+02 	 2.0100e-06 Penalidade (OK)

%Powell             - 2 
%1.87848	 20.23889	 12.81384 	 -1.33e+01 	 -2.04e+01 	 1.2398e-03 Barreira
%1.87836	 20.23676	 12.81256 	  5.61e-01 	  2.40e-01 	 1.8648e-05 Penalidade

%Steepest Descent   - 3 
%1.93026	 19.45994	 13.01069 	 -2.70e+00 	 -8.14e+03 	 3.7050e-04 Barreira (OK)
%1.87836	 20.23676	 12.81256 	  5.61e-01 	  2.40e-01 	 1.8650e-05 Penalidade (ok)

%Fletcher-Reeves    - 4 
%1.93026	 19.45994	 13.01069 	 -2.70e+00 	 -8.14e+03 	 3.7050e-04 Barreira
%1.87836	 20.23676	 12.81256 	  5.61e-01 	  2.40e-01 	 1.8650e-05 Penalidade (ok)

%Newton-Raphson     - 5 
%1.87840	 20.23754	 12.81302 	 -4.42e+00 	 -7.00e+00 	 3.6909e-04 Barreira (OK)
%2.03164	 18.12059	 13.42179 	  1.55e+00 	 -2.47e+04 	 1.1986e-04 Penalidade

%BFGS               - 6
%1.93026	 19.45994	 13.01069 	 -2.70e+00 	 -8.14e+03 	 3.7050e-04 Barreira
%1.87836	 20.23676	 12.81256 	 5.61e-01 	 2.40e-01 	 1.8650e-05 Penalidade (ok)


%% SELECAO DE PARAMETROS E METODOS

%PARAMETROS

TOL = 5e-5;
TOL_BuscaLinear = 1e-3;
Max_iter = 10;

direction_meth = 1;
meth = 1;

%M�todos de Busca Linear
%Univariante        - 1
%Powell             - 2
%Steepest Descent   - 3
%Fletcher-Reeves    - 4
%Newton-Raphson     - 5
%BFGS               - 6
%____________________________
%METODOS OCR
%Barreira           - 1
%Penalidade         - 2 

%% 
switch direction_meth %APENAS PARA MUDAR TITULO DO GRAFICO
    case 1
        titulo = 'Univariante';
        alfa = 1e-4;
        if meth == 1
            TOL = 2e-4;
            alfa = 1e-6;
        end
    case 2
        titulo = 'Powell';
        alfa = 0.0001;
        if meth == 1
            TOL = 2e-3;
        end
    case 3
        titulo = 'Steepest Descent';
        alfa = 1e-7;
        TOL = 6e-5;
    case 4
        titulo = 'Fletcher-Reeves';
        alfa = 1e-7;
        TOL = 1e-4;
    case 5 
        titulo = 'Newton-Raphson';
        alfa = 1e-7;
        TOL = 4e-4;
    case 6
        titulo = 'BFGS';
        alfa =  1e-7;
        TOL = 1e-4;
end

%% DEFINICOES DA FUNCAO E PONTO INICIAL
fprintf('Iter \t xi_x  \t \t xi_y  \t\t xk_x \t\t xk_y \t\t F(xk)  \t c1(xk) \t c2(xk)  \t Crit_parada \t rp \t \n')
    
rho     = 	0.3;
B       =	30;
P       =   33e3;
t       =   0.1;
E       =   3e7;
sigma_y =	1e5;
constants = [rho, B , P, t, E, sigma_y];
alfa_inicial = alfa;

f = @(d,H) (2*rho*pi.*d.*t.*((H.^2 + B^2).^0.5)); %esta certo

c1 = @(d,H) (P*((H^2 + B^2)^0.5)/(pi*d*t*H) - sigma_y);
c2 = @(d,H) (P*((H.^2 + B^2)^0.5)/(pi*d*t*H) - (pi^2*E*(d^2+t^2))/(8*(H^2+B^2)));

switch meth        
    case 1 %Barreira
        
        x0 = [4;25];
        rp = 1e7;
        beta = 0.1;
        
        f_obj = @(d,H,rp)  ((2*rho*pi.*d.*t.*((H.^2+B.^2)^0.5))  - rp./(P*((H.^2+B.^2)^0.5)/(pi*d*t*H) - sigma_y)          - rp./((P*((H.^2 + B^2)^0.5)/(pi*d*t*H) - (pi^2*E*(d^2+t^2)/(8*(H.^2+B.^2))))));
        [x,y] = meshgrid(-1:0.5:10, -5:0.1:45);
        drestricao = linspace(1,4);
        
    case 2 %Penalidade
        x0 = [1;15];
        rp = 1e-7;
        beta = 10;      
                
        f_obj = @(d,H,rp_b) ((2*rho*pi.*d.*t.*((H.^2+B.^2)^0.5)) + 1/2*rp_b(1).*(P*((H.^2+B.^2)^0.5)/(pi*d*t*H) - sigma_y).^2 + 1/2*rp_b(2).*((P*((H.^2 + B^2)^0.5)/(pi*d*t*H) - (pi^2*E*(d^2+t^2)/(8*(H.^2+B.^2)))).^2));       
        [x,y] = meshgrid(-1:0.5:10, -5:0.1:45);
        drestricao = linspace(1,4);
end

z  = f(x,y);
figure(3)
contour(x,y,z,250)
line([0,5],[35,35],[f(1,35), f(5,35)],'Color','red','LineStyle','-') %Plot Restricao
line([0,5],[10,10],[f(0,10), f(5,10)],'Color','red','LineStyle','-') %Plot Restricao
line([1,1],[5,45],[f(1,5), f(1,45)],'Color','red','LineStyle','-') %Plot Restricao
line([3,3],[5,45],[f(3,5), f(3,45)],'Color','red','LineStyle','-') %Plot Restricao
ylim([5 40])
xlim([0 10])
ylabel('H')
xlabel('d')

iter = 1;
xk = x0;

%% TRABALHO 2
rp_b = zeros(1,2);
tic

while true    
    xbefore = xk; %para possibilitar plot
    if meth == 2 %Rp para o METODO Penalidade
        c = [c1(xk(1),xk(2)),c2(xk(1),xk(2))];
        for i = 1:length(c)
            if c(i) >= 0     %RESTRICAO ENTRA em PHI
                rp_b(i) = rp;
            else             %RESTRICAO NAO ENTRA em PHI
                rp_b(i) = 0;               
            end 
        end
       
        
    else %Rb para o METODO da Barreira
        rp_b = rp;
    end
    
    xk = trab1(f_obj, xk, rp_b, alfa, TOL_BuscaLinear, direction_meth, meth,constants);
    
    c = [c1(xk(1),xk(2)),c2(xk(1),xk(2))];
    H = xk(2);
    d = xk(1);
    
    if meth == 1   
        
        if c(1) <= 0 && c(2) <= 0  %VERIFICACAO DE RESTRICAO para METODO da barreira
            %restricao satisfeita
            crit_parada = (- rp./(P*((H.^2+B.^2)^0.5)/(pi*d*t*H) - sigma_y)          - rp./((P*((H.^2 + B^2)^0.5)/(pi*d*t*H) - (pi^2*E*(d^2+t^2)/(8*(H.^2+B.^2))))));
            func = f(xk(1),xk(2));
            
            figure(3)
            hold on
            plot3(xbefore(1),xbefore(2),f(xbefore(1),xbefore(2)),'ro')     % Plot de ponto antes da movimentacao
            plot3(xk(1),xk(2),f(xk(1),xk(2)),'rx')                         % Plot depois da movimentacao
            line([xbefore(1), xk(1)],[xbefore(2), xk(2)],[f(xbefore(1),xbefore(2)) f(xk(1),xk(2))],'Color','red','LineStyle','-.'); %Plot do Caminho
            title(titulo)
            xlabel('d')
            ylabel('H')
            
            hold off
               
           fprintf('%0.1f \t %0.5f \t %0.5f \t %0.5f\t %0.5f\t %0.5f \t %0.2d \t %0.2d \t %0.4d \t %0.2d \n',iter,xbefore(1),xbefore(2),xk(1),xk(2),func,c(1), c(2), crit_parada, rp)
    
           rp = rp*beta;                %Caso nao haja convergencia, atualizacao do rp
           alfa = alfa_inicial;

        else                                %Restricao Violada
            alfa = alfa/2;
            xk = xbefore;                   %Retorna o ponto anterior a busca linear
            iter = iter-1;                  %iteracao descontada
        end
        
    else %CASO 2 PENALIDADE
        
        figure(3)
        hold on
        plot3(xbefore(1),xbefore(2),f(xbefore(1),xbefore(2)),'ro')     % Plot de ponto antes da movimentacao
        plot3(xk(1),xk(2),f(xk(1),xk(2)),'rx')                         % Plot depois da movimentacao    
        line([xbefore(1), xk(1)],[xbefore(2), xk(2)],[f(xbefore(1),xbefore(2)) f(xk(1),xk(2))],'Color','red','LineStyle','-.'); %Plot do Caminho
        title(titulo)
        xlabel('d')
        ylabel('H')
        legend('Xmin', 'x0')
        hold off
        
        crit_parada = 1/2*rp_b(1)*c(1)^2 + 1/2*rp_b(2)*c(2)^2; 
        func = f(xk(1),xk(2));

        fprintf('%0.1f \t %0.5f \t %0.5f \t %0.5f\t %0.5f\t %0.5f \t %0.2d \t %0.2d \t %0.4d \t %0.2d \t %0.2d \n',iter,xbefore(1),xbefore(2),xk(1),xk(2),func,c(1), c(2), crit_parada, rp_b(1),rp_b(2))
        rp = rp*beta;                %Caso nao haja convergencia, atualizacao do rp        
      
    end
    
        if abs(crit_parada) < TOL || iter > Max_iter  %Verificacao de Convergencia
            break
        end
        
    iter = iter + 1;
    
end
toc

%% FUNCTIONS

function Xmin = trab1(f_obj,x0,rp,alfa,TOL,direct_meth,meth,constants)
    Xmin = x0; 
    if direct_meth == 2 %CASO ESPECIAL DE POWELL
        n_iter_powell = 3;
        Step = 1;
        X_inicio_ciclo = x0;
        e1 =  [1;0];
        e2 =  [0;1];
        dk = zeros(2,1); %Apenas pra reservar espaﾃｧo no vetor
        DirectionVector = [e1; e2 ; dk]';
    else
        n_iter_powell = 1;
    end
    
    n_iter = 50;
    for Ciclo = 1:n_iter 
        for k = 1:n_iter_powell %ESTE 'FOR' AUXILIA NO METODO DE POWELL
            switch direct_meth %DEFINICAO DA DIRECAO 
                case 1
                    % Univariante 
                     if mod(Ciclo,2) == 1
                          dk=[1;0];
                     else
                          dk=[0;1];
                     end
                case 2
                    % Powell
                    if Step == 1
                        if Ciclo == 1
                            dk = DirectionVector(1:2)';
                        else
                            dk = DirectionVector(3:4)';
                            X_inicio_ciclo = Xmin;
                        end
                        Step = Step + 1;
                    elseif Step == 2
                        Step = Step + 1;
                        if Ciclo == 1
                            dk = DirectionVector(3:4)';
                        else
                            dk = DirectionVector(5:6)';     
                        end 
                    else
                        Step = 1;            
                        dk = Xmin - X_inicio_ciclo;
                        if DirectionVector(5:6) == [0,0] %atualizaﾃδｧﾃδ｣o do vetor no primeiro ciclo ﾃδｩ diferente
                            DirectionVector(5:6) = dk;
                        else
                            DirectionVector(1:2) = DirectionVector(3:4);
                            DirectionVector(3:4) = DirectionVector(5:6);
                            DirectionVector(5:6) = dk;  
                        end    
                    end
                case 3
                    % Steepest Descent
                    dk = Gradiente(Xmin,rp,meth,constants,f_obj);
                case 4
                    % Fletcher-Reeves
                    if k ==1
                        dk = -Gradiente(Xmin,rp, meth,constants);
                        g0 = -dk;
                    else
                        g1 = Gradiente(Xmin,rp, meth,constants);
                        bk = g1.'*g1/(g0.'*g0);
                        dk = -g1 + bk*dk;
                        g0 = g1;
                    end
                case 5
                    % Newton-Raphson
                    dk = -(Hesiana(Xmin, rp, meth,constants,f_obj))^-1*Gradiente(Xmin,rp, meth,constants);
                case 6
                    % BFGS
                    if k == 1
                        S = 1;
                        dk = -(S*Gradiente(Xmin,rp,meth,constants));
                    else
                        SigmaX    =  xk - xbefore;
                        SigmaGrad =  Gradiente(Xmin,rp,meth,constants) - Gradiente(xbefore,rp,meth,constants);
                        A = ((SigmaX'*SigmaGrad)+(SigmaGrad'*S*SigmaGrad))*(SigmaX*SigmaX')/((SigmaX'*SigmaGrad))^2;
                        B = (S*(SigmaGrad*SigmaX')+(SigmaX*(S*SigmaGrad)'))/(SigmaX'*SigmaGrad);
                        S = S + A - B;
                        dk = -(S*Gradiente(Xmin,rp,meth,constants)); 
                    end               
            end 
            
            [xL, xD] = PassoCte(Xmin,dk,alfa,rp,f_obj);
            Xmin = Aurea(xL,xD,dk,TOL,rp,f_obj);
        end
        
        if norm(Gradiente(Xmin,rp,meth,constants,f_obj)) < TOL %criterio de parada
            break
        end
    end
end


function [xL, xD] = PassoCte(x0,dk,alfa,rp,f_obj)
    Deltax = dk*alfa*0.1;
    x0dir = x0 + Deltax;
    x0esq = x0 - Deltax; 
    if f_obj(x0dir(1),x0dir(2),rp) > f_obj(x0esq(1),x0esq(2),rp) %VERIFICACAO DE DIRECAO
        direction = -dk;
    else
        direction = dk;
    end
    xj = x0 + alfa*direction;
    while f_obj(x0(1),x0(2),rp) > f_obj(xj(1),xj(2),rp) %Busca pelo primeiro MINIMO     
        x0 = x0 + alfa*direction;
        xj = xj + alfa*direction;                
    end
    xL = x0; %Limite POSICAO a "esquerda"
    xD = xj; %Limite POSICAO a "direita" 
end



function [xmin] = Aurea(xL,xD,dk,TOL,rp,f_obj)
    RA = (sqrt(5)-1)/2;
    beta = pdist2(xL',xD');
    dk = dk/norm(dk);
    
    while beta > TOL
        
        alfaEsq = xL + dk*(1-RA)*beta;
        alfaDir = xL + dk*RA*beta;

        if f_obj(alfaEsq(1),alfaEsq(2),rp) > f_obj(alfaDir(1),alfaDir(2),rp)
            xL = alfaEsq;     
        else
            xD = alfaDir;
        end        
        beta = pdist2(xL',xD');
    end
    xmin = (xL+xD)/2;
end

function [D] = Gradiente(x,rp,meth,constants,f_obj)
    
    %USADO FUNCAO GRADIENT DO MATLAB COM VARIAVEL SYMS
    d = x(1);
    H = x(2);
    
    D    =   zeros(2,1);

    switch meth %CONTRIBUICAO
        case 1    % Barreira
%             syms d H rp
%             G = subs(gradient(f_obj(d,H,rp)),[d,H,rp],[x(1),x(2),rp_bp]); 
%             D(1,1) = G(2,1);
%             D(2,1) = G(1,1);
            D(1) = ((3*pi*(H^2 + 900)^(1/2))/50) + (((16*H*((1241883636537201*d^2)/4194304 + 1241883636537201/419430400))/(8*H^2 + 7200)^2 + 330000/(d*pi*(H^2 + 900)^(1/2)) - (330000*(H^2 + 900)^(1/2))/(H^2*d*pi))/(((1241883636537201*d^2)/4194304 + 1241883636537201/419430400)/(8*H^2 + 7200) - (330000*(H^2 + 900)^(1/2))/(H*d*pi))^2 + (rp*(330000/(d*pi*(H^2 + 900)^(1/2)) - (330000*(H^2 + 900)^(1/2))/(H^2*d*pi)))/((330000*(H^2 + 900)^(1/2))/(H*d*pi) - 100000)^2);
            D(2) = (3*pi*H*d)/(50*(H^2 + 900)^(1/2)) - ((1241883636537201*d)/(2097152*(8*H^2 + 7200)) + (330000*(H^2 + 900)^(1/2))/(H*d^2*pi))/(((1241883636537201*d^2)/4194304 + 1241883636537201/419430400)/(8*H^2 + 7200) - (330000*(H^2 + 900)^(1/2))/(H*d*pi))^2 - (330000*rp*(H^2 + 900)^(1/2))/(H*d^2*pi*((330000*(H^2 + 900)^(1/2))/(H*d*pi) - 100000)^2);
        case 2 %Penalidade
            D(1) =  rp(1)*((330000*(H^2 + 900)^(1/2))/(H*d*pi) - 100000)*(330000/(d*pi*(H^2 + 900)^(1/2)) - (330000*(H^2 + 900)^(1/2))/(H^2*d*pi)) - rp(2)*(((1241883636537201*d^2)/4194304 + 1241883636537201/419430400)/(8*H^2 + 7200) - (330000*(H^2 + 900)^(1/2))/(H*d*pi))*((16*H*((1241883636537201*d^2)/4194304 + 1241883636537201/419430400))/(8*H^2 + 7200)^2 + 330000/(d*pi*(H^2 + 900)^(1/2)) - (330000*(H^2 + 900)^(1/2))/(H^2*d*pi)) + (3*pi*H*d)/(50*(H^2 + 900)^(1/2));
            D(2) = (3*pi*(H^2 + 900)^(1/2))/50 + rp(2)*(((1241883636537201*d^2)/4194304 + 1241883636537201/419430400)/(8*H^2 + 7200) - (330000*(H^2 + 900)^(1/2))/(H*d*pi))*((1241883636537201*d)/(2097152*(8*H^2 + 7200)) + (330000*(H^2 + 900)^(1/2))/(H*d^2*pi)) - (330000*rp(1)*(H^2 + 900)^(1/2)*((330000*(H^2 + 900)^(1/2))/(H*d*pi) - 100000))/(H*d^2*pi);
    end
end

function [Hessian] = Hesiana(x,rp_bp,meth,constants,f_obj)
    Hessian = zeros(2,2);
    
    rho     = 	constants(1);
    B       =	constants(2);
    P       =   constants(3);
    t       =   constants(4);
    E       =   constants(5);
    sigma_y =	constants(6);
    
     switch meth %Contribui????o da Restri????o
        case 1    
            syms d H rp  
            G = subs(hessian(f_obj(d,H,rp)),[d,H,rp],[x(1),x(2),rp_bp]); 
            Hessian(1,1) = G(2,2);
            Hessian(2,1) = G(2,1);
            Hessian(1,2) = G(1,2);
            Hessian(2,2) = G(1,1);
            
         case 2 
            f_obj = @(d,H,rp_a,rp_b) ((2*rho*pi.*d.*t.*((H.^2+B.^2)^0.5)) + 1/2*rp_a.*(P*((H.^2+B.^2)^0.5)/(pi*d*t*H) - sigma_y).^2 + 1/2*rp_b.*((P*((H.^2 + B^2)^0.5)/(pi*d*t*H) - (pi^2*E*(d^2+t^2)/(8*(H.^2+B.^2)))).^2));
            
            syms d H rp_a rp_b
            G = subs(hessian(f_obj(d,H,rp_a,rp_b)),[d,H,rp_a,rp_b],[x(1),x(2),rp_bp(2),rp_bp(2)]); 
            Hessian(1,1) = G(1,1);
            Hessian(2,1) = G(2,1);
            Hessian(1,2) = G(1,2);
            Hessian(2,2) = G(2,2);
    end

end