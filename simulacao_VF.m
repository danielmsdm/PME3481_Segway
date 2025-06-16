% Código principal do projeto

% Andre Hatano Feris - 11382644 
% Daniel Marins Silva D'Oliveira Martins - 12554571
% Felipe Almeida Ribeiro - 12566564 
% Lucas Junji Koreeda - 11347262

close all;
clc;
clear;

% Definição das constantes
% Características do sistema
Mp = 100;      % [kg]: massa do conjunto (Segway + Pessoa)
Mrr = 7.5;      % [kg]: massa da roda
R = 0.2;      % [m]: raio da roda
D = 0.6;      % [m]: distância entre centros das rodas
L = 0.6;      % [m]: distância entre o eixo das rodas e o centro de massa
g = 9.81;     % [m/s²]: aceleração da gravidade

% Cálculos de propriedades
Jrr = Mrr*(R^2)/2;      % [kg·m²]: momento de inércia da roda
Jpt = Mp*(L^2);         % [kg·m²]: momento de inércia do eixo de arfagem
Jpy = Mp*((0.5*D)^2)/2; % [kg·m²]: momento de inércia do eixo de guinada

% Características do motor
N = 24;     % [-]: razão de redução do motor
ke = 0.05;  % [Vs/rad]: constante eletromotriz do motor
kt = ke;    % [Nm/A]: constante de torque do motor
Ra = 1;     % [OHm]: Resistência do enrolamento do motor

% Modelo matemático
k1 = (N * kt) / Ra;
k2 = ((N^2) * kt * ke) / (R * Ra);

m = 2 * (Mp*Jrr*(L^2) + Mp*Mrr*(R^2)*(L^2) + Jrr*Jpt + Jpt*Mrr*(R^2));
n = 2*Jpy*(R^2) - Jrr*(D^2) - Mrr*(R^2)*(D^2);

A22 = -2*k2 * (Mp*R*(L^2) + Mp*(R^2)*L + Jpt*R) / m;
A23 = -((Mp^2)*(R^2)*(L^2)*g) / m;
A42 = -2*k2/(Jpt + (Mp*L));
A43 = Mp*g*L / (Jpt + (Mp*L));
A66 = -k2*(D^2)*R / n;
B21 = (k1 * (Mp * R * (L^2) + Mp * (R^2) * L + Jpt * R)) / m;
B41 = k1 / (Jpt + (Mp*L));
B62 = -k1*D*R / n;

B22 = B21;
B42 = B41;
B61 = B62;
B62_n = -B62;

% Espaço de estados
A = [0 1   0    0    0   0;
     0 A22 A23  0    0   0;
     0 0   0    1    0   0;
     0 A42 A43  0    0   0;
     0 0   0    0    0   1;
     0 0   0    0    0  A66];

B = [0   0;
     B21 B22;
     0   0;
     B41 B42;
     0   0;
     B61  B62_n];

C = eye(6);
C(2,2) = 0;
C(4,4) = 0;
C(6,6) = 0;

D = zeros(6,2);

% Sistema em espaço de estados (malha aberta)
sys_MA = ss(A, B, C, D);

% Funções de transferência (malha aberta)
n = size(A,1);
[m, p] = size(D);
den = poly(A);

for j = 1:p % loop para obter e exibir as FTs
  [numj, ~] = ss2tf(A, B, C, D, j);
  fprintf('--- Funções de Transferência da entrada %d ---\n\n', j);
  for i = 1:m
    b = numj(i,:);
    if all(b==0)
      fprintf('G_{%d,%d}(s) = 0\n\n', i, j);
      continue
    end

    num_str = poly2str(b, 's');
    den_str = poly2str(den, 's');

    w = max(length(num_str), length(den_str));
    pad_num = floor((w-length(num_str))/2);
    pad_den = floor((w-length(den_str))/2);

    fprintf('G_{%d,%d}(s) =\n', i, j);
    fprintf('%s%s\n', repmat(' ',1,pad_num), num_str);
    fprintf('%s%s\n\n', repmat('-',1,w), repmat(' ',1,pad_den)),
    fprintf('%s%s\n\n', repmat(' ',1,pad_den), den_str);
  end
end

% Diagrama de polos e zeros (malha aberta)
disp('Polos malha aberta:');
disp(pole(sys_MA));
figure;
pzmap(sys_MA);
grid on;

% Controlabilidade e observabilidade
Co = ctrb(A, B);       
rank_Co = rank(Co);
disp(['Posto da Matriz de Controlabilidade: ', num2str(rank_Co)]);

Ob = obsv(A, C);           
rank_Ob = rank(Ob);        
disp(['Posto da Matriz de Observabilidade: ', num2str(rank_Ob)]);

% Controle por LQR
% Matrizes de Ponderação
Q_LQR = diag([500, 0.01, 10, 0.01, 10, 0.01]); 
R_LQR = diag([10, 10]);

[K_LQR, S, e] = lqr(A, B, Q_LQR, R_LQR);
disp('Matriz de ganhos (LQR): ');
disp(K_LQR);

% Malha fechada com LQR
A_MF_LQR = A - B*K_LQR;
sys_MF_LQR = ss(A_MF_LQR, zeros(size(B)),C ,D);
disp('Polos em malha fechada (LQR): ');
disp(eig(A - B*K_LQR));

% Comparação entre polos malha aberta x malha fechada (LQR)
figure;
pzmap(sys_MA, sys_MF_LQR);
grid on;
title('Comparação dos Polos: Malha Aberta x LQR');
legend('Malha Aberta', 'LQR');

% Simulação com condições iniciais
v_0 = 1;
thetap_0 = deg2rad(-1);
thetay_0 = deg2rad(45);

x0 = [0; v_0; thetap_0; 0; thetay_0; 0];
t = 0:0.01:10;
[y_LQR, t, x_lqr_ideal] = initial(sys_MF_LQR, x0, t);
u_lqr_ideal = -K_LQR * x_lqr_ideal';

% Malha fechada com alocação de polos
% Polos desejados
p1 = -0.8;
p2 = -1.0;
p3 = -2.5;
p4 = -2.7;
p5 = -3.0;
p6 = -3.5;
p_pp = [p1, p2, p3, p4, p5, p6];

K_pp = place(A, B, p_pp);
disp('Matriz de ganhos (alocação de polos): ');
disp(K_pp);

% Sistema em malha fechada (alocação de polos)
A_MF_pp = A - B*K_pp;
sys_MF_pp = ss(A_MF_pp, zeros(size(B)),C ,D);
disp('Polos em malha fechada (alocação de polos): ');
disp(eig(A - B*K_pp));

% Comparação entre polos LQR x alocação de polos
figure;
pzmap(sys_MF_pp, sys_MF_LQR);
grid on;
title('Comparação dos Polos: Alocação de Polos x LQR');
legend('Alocação de Polos', 'LQR');

% Simulação com a mesma condição inicial
[y_pp, ~, x_pp_ideal] = initial(sys_MF_pp, x0, t);
u_pp_ideal = -K_pp * x_pp_ideal';

% Observadores identidade
p_obs_LQR = [-0.8, -1.0, -1.2, -1.6, -1.2, -1.4];
L_lqr = place(A', C', p_obs_LQR)';
disp('Matriz de ganhos do observador (LQR): ');
disp(L_lqr);

p_obs_pp = [-0.8, -1.0, -1.2, -1.6, -1.2, -1.4];
L_pp = place(A', C', p_obs_pp)';
disp('Matriz de ganhos do observador (alocação de polos): ');
disp(L_pp);

% Simulação com observador de estados
dt = 0.01;
t_final = 10;
t_sim = 0:dt:t_final; 
n_pontos = length(t_sim);
n_states = size(A,1);

% Posição inicial do sistema real
x0_real = x0;
% O observador começa sabendo apenas o que podemos medir
x0_estimated = [x0(1); 0; x0(3); 0; x0(5); 0];

% Para o sistema com LQR + Observador
x_real_lqr_obs = zeros(n_states, n_pontos);
x_est_lqr_obs  = zeros(n_states, n_pontos);
u_lqr_obs      = zeros(size(B,2), n_pontos);
x_real_lqr_obs(:,1) = x0_real;
x_est_lqr_obs(:,1)  = x0_estimated;

% Para o sistema com Alocação de Polos + Observador
x_real_pp_obs = zeros(n_states, n_pontos);
x_est_pp_obs  = zeros(n_states, n_pontos);
u_pp_obs      = zeros(size(B,2), n_pontos);
x_real_pp_obs(:,1) = x0_real;
x_est_pp_obs(:,1)  = x0_estimated;

% Loop de simulação dos 2 métodos de controle
for k = 1:(n_pontos - 1)
    % Simulação para LQR + Observador
    y_medido_lqr = C * x_real_lqr_obs(:,k);
    u_lqr_obs(:, k) = -K_LQR * x_est_lqr_obs(:,k);
    dx_real_lqr = A * x_real_lqr_obs(:,k) + B * u_lqr_obs(:, k);
    x_real_lqr_obs(:, k+1) = x_real_lqr_obs(:,k) + dx_real_lqr * dt;
    erro_med_lqr = y_medido_lqr - C * x_est_lqr_obs(:,k);
    dx_hat_lqr = A * x_est_lqr_obs(:,k) + B * u_lqr_obs(:, k) + L_lqr * erro_med_lqr;
    x_est_lqr_obs(:, k+1) = x_est_lqr_obs(:,k) + dx_hat_lqr * dt;

    % Simulação para Alocação de Polos + Observador
    y_medido_pp = C * x_real_pp_obs(:,k);
    u_pp_obs(:, k) = -K_pp * x_est_pp_obs(:,k);
    dx_real_pp = A * x_real_pp_obs(:,k) + B * u_pp_obs(:, k);
    x_real_pp_obs(:, k+1) = x_real_pp_obs(:,k) + dx_real_pp * dt;
    erro_med_pp = y_medido_pp - C * x_est_pp_obs(:,k);
    dx_hat_pp = A * x_est_pp_obs(:,k) + B * u_pp_obs(:, k) + L_pp * erro_med_pp;
    x_est_pp_obs(:, k+1) = x_est_pp_obs(:,k) + dx_hat_pp * dt;
end

% Calcula o último ponto do sinal de controle
u_lqr_obs(:,n_pontos) = -K_LQR * x_est_lqr_obs(:,n_pontos); 
u_pp_obs(:,n_pontos) = -K_pp * x_est_pp_obs(:,n_pontos);

% Gráfico Comparativo dos estados
figure('Name', 'Comparação Final: Ideal vs. Observador');
sgtitle('Comparação dos Estados: Ideal vs. Observador', 'FontSize', 14, 'FontWeight', 'bold');
titles = {'x_1: posição [m]', 'x_2: velocidade [m/s]', 'x_3: ângulo de arfagem [°]', ...
          'x_4: vel. angular de arfagem [°/s]', 'x_5: ângulo de guinada [°]', 'x_6: vel. angular de guinada [°/s]'};

for i = 1:6
    subplot(3,2,i);
    
    if ismember(i, [3 4 5 6]) % converter ângulos para graus
        data_lqr_ideal = rad2deg(x_lqr_ideal(:,i));
        data_pp_ideal  = rad2deg(x_pp_ideal(:,i));
        data_lqr_obs   = rad2deg(x_real_lqr_obs(i,:));
        data_pp_obs    = rad2deg(x_real_pp_obs(i,:));
    else
        data_lqr_ideal = x_lqr_ideal(:,i);
        data_pp_ideal  = x_pp_ideal(:,i);
        data_lqr_obs   = x_real_lqr_obs(i,:);
        data_pp_obs    = x_real_pp_obs(i,:);
    end
    
    % Plotando os gráficos
    plot(t, data_lqr_ideal, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, data_pp_ideal, 'r-', 'LineWidth', 1.5); hold on;
    plot(t_sim, data_lqr_obs, 'c--', 'LineWidth', 1.5); hold on;
    plot(t_sim, data_pp_obs, 'm--', 'LineWidth', 1.5);
    
    title(titles{i});
    xlabel('Tempo [s]');
    grid on;
    if i == 1
        legend('LQR (Ideal)', 'Alocação Polos (Ideal)', 'LQR (Observador)', 'Alocação Polos (Observador)', 'Location', 'southeast');
    end
end

figure('Name', 'Comparação dos Sinais de Controle');
sgtitle('Comparação dos Sinais de Controle: Ideal vs. Observador');

% Controle u1
subplot(2,1,1);
plot(t, u_lqr_ideal(1,:), 'b-', t, u_pp_ideal(1,:), 'r-'); hold on;
plot(t_sim, u_lqr_obs(1,:), 'c--', t_sim, u_pp_obs(1,:), 'm--');
ylabel('u_1 [Torque]');
legend('LQR (Ideal)', 'Alocação (Ideal)', 'LQR (Obs)', 'Alocação (Obs)');
grid on;

% Controle u2
subplot(2,1,2);
plot(t, u_lqr_ideal(2,:), 'b-', t, u_pp_ideal(2,:), 'r-'); hold on;
plot(t_sim, u_pp_obs(2,:), 'm--', t_sim, u_lqr_obs(2,:), 'c--');
xlabel('Tempo [s]'); ylabel('u_2 [Torque]');
legend('LQR (Ideal)', 'Alocação (Ideal)', 'LQR (Obs)', 'Alocação (Obs)');
grid on;

% Observador de ordem reduzida
% Matriz de medição C: são medidos x1 (posição), x3 (arfagem) e x5 (guinada)
C_medicao = [1 0 0 0 0 0;  
             0 0 1 0 0 0;  
             0 0 0 0 1 0];

% Matriz V: estados não medidos são: x2, x4, x6.
V = [0 1 0 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 0 1];

% Matriz de transformação P
P = [V; C_medicao];
P_inv = inv(P);

% Matrizes do sistema transformado
A_til = P * A * P_inv;
B_til = P * B;

% Matrizes A_til e B_til
n_est = 3; % estados estimados: x2, x4, x6
A11 = A_til(1:n_est, 1:n_est);
A12 = A_til(1:n_est, n_est+1:end);
A21 = A_til(n_est+1:end, 1:n_est);
A22 = A_til(n_est+1:end, n_est+1:end);
B1 = B_til(1:n_est, :);
B2 = B_til(n_est+1:end, :);

% Escolha dos polos do observador
p_obs_red = [-0.9, -1, -1.1]; % mais rápidos que os polos do controlador
L = place(A11', A21', p_obs_red)';
disp('Ganho do observador de ordem reduzida: ');
disp(L)

% Matriz M usada para no loop de simulação
M = (A11*L - L*A21*L + A12 - L*A22);

% Comparação entre polos malha fechada (alocação de polos) e polos do
% observador de ordem reduzida
A_sys_obs = diag(p_obs_red);                % matriz diagonal com os polos
B_sys_obs = zeros(length(p_obs_red), 1);    % matriz vetor coluna de zeros
C_sys_obs = eye(length(p_obs_red));         % matriz identidade
D_sys_obs = zeros(length(p_obs_red), 1);    % matriz nula
sys_obs_p = ss(A_sys_obs, B_sys_obs, C_sys_obs, D_sys_obs);

figure; 
pzmap(sys_MF_pp, sys_obs_p);
grid on;
title('Comparação dos Polos: Controlador (Alocação de Polos) vs. Observador');
legend('Controlador (Alocação de Polos)', 'Observador (Ordem Reduzida)');

% Gráfico para convergência do observador
figure('Name', 'Convergência do Observador de Ordem Reduzida');
sgtitle('Convergência dos Estados Não Medidos: Real vs. Estimado');

% Plotando x2 (velocidade longitudinal)
subplot(3,1,1);
plot(t_sim, x_real_pp_obs(2,:), 'b-', 'LineWidth', 1.5); hold on;
plot(t_sim, x_est_pp_obs(2,:), 'r--', 'LineWidth', 1.5);
title('x_2: velocidade [m/s]');
xlabel('Tempo [s]');
ylabel('Velocidade [m/s]');
legend('Real', 'Estimado');
grid on;

% Plotando x4 (velocidade angular de arfagem)
subplot(3,1,2);
plot(t_sim, rad2deg(x_real_pp_obs(4,:)), 'b-', 'LineWidth', 1.5); hold on;
plot(t_sim, rad2deg(x_est_pp_obs(4,:)), 'r--', 'LineWidth', 1.5);
title('x_4: vel. angular de arfagem [°/s]');
xlabel('Tempo [s]');
ylabel('Velocidade Angular [°/s]');
legend('Real', 'Estimado');
grid on;

% Plotando x6 (velocidade angular de guinada)
subplot(3,1,3);
plot(t_sim, rad2deg(x_real_pp_obs(6,:)), 'b-', 'LineWidth', 1.5); hold on;
plot(t_sim, rad2deg(x_est_pp_obs(6,:)), 'r--', 'LineWidth', 1.5);
title('x_6: vel. angular de guinada [°/s]');
xlabel('Tempo [s]');
ylabel('Velocidade Angular [°/s]');
legend('Real', 'Estimado');
grid on;

% Seguidor de referência e rejeição de perturbação
% Cálculo da pré-alimentação
C_track = [1 0 0 0 0 0;   % referências: x1 (posição) e x5 (guinada)
           0 0 0 0 1 0];
ganho_dc = C_track * inv(-(A - B*K_pp)) * B;
Nr = pinv(ganho_dc);

% Simulação dos diferentes casos
for caso = 1:2 
    dt = 0.01;
    t_final = 10;
    t_sim = 0:dt:t_final;
    n_pontos = length(t_sim);
    x_real = zeros(6,1);
    
    % Inicialização do observador
    y_medido_inicial = C_medicao * x_real;
    w_est = zeros(3,1); % estimativa inicial dos estados não medidos é zero
    v_est = w_est - L*y_medido_inicial;
    
    x_til_est = [w_est; y_medido_inicial];
    x_estimado = P_inv * x_til_est;
    
    % Históricos
    x_hist = zeros(6, n_pontos);
    x_hist(:,1) = x_real;

    if caso == 1 % acompanhamento de referência
        r = [2; deg2rad(30)]; % ir para x=2m e virar 30 graus
    else 
        % rejeição de perturbações
        r = [0; 0]; % ficar parado
    end
    
    % Loop de simulação principal
    for k = 1:(n_pontos - 1)
        if caso == 2 && abs(t_sim(k) - 2) < dt/2
            x_real(2) = x_real(2) + 1; % perturbação na velocidade de 1m/s
        end

        u = -K_pp * x_estimado + Nr * r; % lei de controle
        
        % Simula a planta
        dx_real = A*x_real + B*u;
        x_real = x_real + dx_real*dt;
        
        % Simula o observador de ordem reduzida
        y_medido = C_medicao * x_hist(:,k);
        dv_est = (A11 - L*A21)*v_est + M*y_medido + (B1 - L*B2)*u;
        v_est = v_est + dv_est*dt;
        
        % Reconstitui a estimativa do estado completo
        w_est = v_est + L*y_medido;
        x_til_est = [w_est; y_medido];
        x_estimado = P_inv * x_til_est;
        
        % Salva o estado real no histórico
        x_hist(:,k+1) = x_real;
    end

    % Plot dos resultados
    figure;
    if caso == 1
        sgtitle('Caso 1: Acompanhamento de Referência');
    else
        sgtitle('Caso 2: Rejeição de Perturbações');
    end
    
    titles = {'x_1: posição', 'x_2: velocidade', 'x_3: ângulo arfagem', 'x_4: vel. ang. arfagem', 'x_5: ângulo guinada', 'x_6: vel. ang. guinada'};
    if caso == 1; refs_plot = {r(1), [], [], [], r(2), []}; else; refs_plot = {0,0,0,0,0,0}; end

    for i = 1:6
        subplot(3,2,i);
        data_to_plot = x_hist(i,:);
        ref_val = refs_plot{i};
        
        if ismember(i, [3, 4, 5, 6]); data_to_plot = rad2deg(data_to_plot); if ~isempty(ref_val); ref_val = rad2deg(ref_val); end; end

        plot(t_sim, data_to_plot, 'b', 'LineWidth', 1.5); hold on;
        
        if caso == 1 && ~isempty(ref_val)
            plot(t_sim, ref_val*ones(size(t_sim)), 'r--', 'LineWidth', 1.2);
            legend('Resposta', 'Referência', 'Location', 'best');
        elseif caso == 2
             plot(t_sim, zeros(size(t_sim)), 'r--', 'LineWidth', 1.2);
             xline(2, 'm--', 'Perturbação');
        end

        title(titles{i}); xlabel('Tempo [s]'); grid on;
    end
end

