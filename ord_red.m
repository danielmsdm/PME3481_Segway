% Código para comparar os observadores

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

% Controle por LQR
% Matrizes de Ponderação
Q_LQR = diag([500, 0.01, 10, 0.01, 10, 0.01]); 
R_LQR = diag([10, 10]);
[K_LQR, S, e] = lqr(A, B, Q_LQR, R_LQR);

% Malha fechada com LQR
A_MF_LQR = A - B*K_LQR;
sys_MF_LQR = ss(A_MF_LQR, zeros(size(B)),C ,D);

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

% Sistema em malha fechada (alocação de polos)
A_MF_pp = A - B*K_pp;
sys_MF_pp = ss(A_MF_pp, zeros(size(B)),C ,D);

% Simulação com a mesma condição inicial
[y_pp, ~, x_pp_ideal] = initial(sys_MF_pp, x0, t);
u_pp_ideal = -K_pp * x_pp_ideal';

% Observador de ordem reduzida
% Vamos definir V para que P seja inversível
V = [0 1 0 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 0 1];

% A matriz de medições C deve ser simplificada:
C = [1 0 0 0 0 0;  
     0 0 1 0 0 0;  
     0 0 0 0 1 0]; 

% Definindo P
P = [V; C];
P_inv = inv(P);

% Aplicando as transformações no sistema
A_til = P * A * P_inv;
B_til = P * B;

% Separando as matrizes 6x6 em blocos 3x3
n_states = 6;   % todos os estados
n_est = 3;      % estados estimados
n_med = 3;      % estados medidos

A11 = A_til(1:n_est, 1:n_est);
A12 = A_til(1:n_est, n_est+1:end);
A21 = A_til(n_est+1:end, 1:n_est);
A22 = A_til(n_est+1:end, n_est+1:end);

B1 = B_til(1:n_est, :);
B2 = B_til(n_est+1:end, :);

% Polos do observador
p_obs_red = [-0.9, -1, -1.1];

% Ganho do observador
L = place(A11', A21', p_obs_red)';
disp('Ganho L do observador de ordem reduzida:');
disp(L)

% Observador identidade (para comparação)
observer_poles_id = [-0.8, -1.0, -1.2, -1.6, -1.2, -1.4];
% O cálculo do ganho L_id usa as matrizes originais A e C (3x6)
L_id = place(A', C', observer_poles_id)';

dt = 0.01;
t_sim = 0:dt:max(t);
n_pontos = length(t_sim);
n_states = 6;
n_est = 3;

% Vetores para o caso com Observador de Ordem Reduzida
x_real_ro = zeros(n_states, n_pontos);
x_est_ro  = zeros(n_states, n_pontos);
u_ro      = zeros(size(B,2), n_pontos);
x_est_v       = zeros(n_est, 1);

% Vetores para o caso Ideal (agora simulado com o mesmo método)
x_ideal_loop = zeros(n_states, n_pontos); 
u_ideal_loop = zeros(size(B,2), n_pontos);
x_real_id = zeros(n_states, n_pontos);
x_est_id = zeros(n_states, n_pontos);
u_id      = zeros(size(B,2), n_pontos);

% Condições iniciais
x_real_ro(:,1) = x0;
x_ideal_loop(:,1)  = x0; % O caso ideal também parte de x0

x0_estimated = [x0(1); 0; x0(3); 0; x0(5); 0];

x_real_id(:,1) = x0;
x_est_id(:,1)  = x0_estimated;

x_est_v(:,1)       = zeros(n_est, 1);
y0 = C * x0;
w_est_0 = x_est_v(:,1) + L*y0;
x_til_est_0 = [w_est_0; y0];
x_est_ro(:,1) = P_inv * x_til_est_0;

% Loop de Simulação Unificado
for k = 1:(n_pontos - 1)
    
    % Opção 1: Simulação do Sistema com Observador
    u_ro(:, k) = -K_LQR * x_est_ro(:, k);
    dx_real = A * x_real_ro(:,k) + B * u_ro(:, k);
    x_real_ro(:, k+1) = x_real_ro(:,k) + dx_real * dt;
    
    y_k = C * x_real_ro(:,k);
    M = (A11*L - L*A21*L + A12 - L*A22);
    dv_hat = (A11 - L*A21)*x_est_v(:,k) + M*y_k + (B1 - L*B2)*u_ro(:,k);
    x_est_v(:,k+1) = x_est_v(:,k) + dv_hat * dt;
    
    y_k_plus_1 = C * x_real_ro(:,k+1);
    w_hat_k_plus_1 = x_est_v(:,k+1) + L*y_k_plus_1;
    x_tilde_hat_k_plus_1 = [w_hat_k_plus_1; y_k_plus_1];
    x_est_ro(:, k+1) = P_inv * x_tilde_hat_k_plus_1;
    
    % Opção 2: Simulação do Sistema Ideal (usando o mesmo método de Euler)
    u_ideal_loop(:, k) = -K_LQR * x_ideal_loop(:,k); 
    dx_ideal = A * x_ideal_loop(:,k) + B * u_ideal_loop(:,k);
    % Forma mais simples: dx_ideal = (A - B*K_LQR) * x_ideal_loop(:,k);
    x_ideal_loop(:, k+1) = x_ideal_loop(:,k) + dx_ideal * dt;

    % Opção 3: Simulação do Sistema com Observador de Ordem Completa
    u_id(:, k) = -K_LQR * x_est_id(:, k);
    dx_real_id = A * x_real_id(:, k) + B * u_id(:, k);
    x_real_id(:, k+1) = x_real_id(:, k) + dx_real_id * dt;
    
    % Atualização do observador de ordem completa
    y_id_k = C * x_real_id(:, k);
    dx_est_id = A * x_est_id(:, k) + B * u_id(:, k) + L_id * (y_id_k - C * x_est_id(:, k));
    x_est_id(:, k+1) = x_est_id(:, k) + dx_est_id * dt;


end
% Calcula último ponto dos controles
u_ro(:,n_pontos) = -K_LQR * x_est_ro(:, n_pontos);
u_ideal_loop(:,n_pontos) = -K_LQR * x_ideal_loop(:, n_pontos);
u_id(:,n_pontos) = -K_LQR * x_est_id(:, n_pontos);

% Comparação de Desempenho (Ideal vs. Observadores)
figure('Name', 'Comparação Final de Desempenho');
sgtitle('Desempenho: Ideal vs. Obs. Reduzida vs. Obs. Completa', 'FontSize', 14);
titles = {'x_1: posição [m]', 'x_2: velocidade [m/s]', 'x_3: ângulo de arfagem [°]','x_4: vel. angular de arfagem [°/s]', 'x_5: ângulo de guinada [°]', 'x_6: vel. angular de guinada [°/s]'};

for i = 1:6
    subplot(3,2,i);
    
    % Plotando as 3 curvas
    if ismember(i, [3, 4, 5, 6])
        plot(t_sim, rad2deg(x_ideal_loop(i,:)), 'r-', 'LineWidth', 2); hold on;
        plot(t_sim, rad2deg(x_real_ro(i,:)), 'm--', 'LineWidth', 1.5);
        plot(t_sim, rad2deg(x_real_id(i,:)), 'g:', 'LineWidth', 2.5); 
    else
        plot(t_sim, x_ideal_loop(i,:), 'r-', 'LineWidth', 2); hold on; 
        plot(t_sim, x_real_ro(i,:), 'm--', 'LineWidth', 1.5);
        plot(t_sim, x_real_id(i,:), 'g:', 'LineWidth', 2.5); 
    end
    
    title(titles{i}); xlabel('Tempo [s]'); grid on;
    
    % Legenda atualizada com as 3 entradas
    if i==1
        legend('Ideal', 'Obs. Ordem Reduzida', 'Obs. Ordem Completa', 'Location', 'best');
    end
end

figure('Name', 'Comparação dos Sinais de Controle');
sgtitle('Sinais de Controle: Ideal vs. Observadores', 'FontSize', 14);

%  u1: torque do primeiro motor
subplot(2,1,1);
plot(t_sim, u_ideal_loop(1,:), 'r-', 'LineWidth', 2); 
hold on;
plot(t_sim, u_ro(1,:), 'm--', 'LineWidth', 1.5);
plot(t_sim, u_id(1,:), 'g:', 'LineWidth', 2.5);
hold off;
ylabel('u_1 [N.m]');
grid on;
legend('Ideal', 'Obs. Ordem Reduzida', 'Obs. Ordem Completa', 'Location', 'northeast');
title('Sinal de Controle u_1');

%  u2: torque do segundo motor
subplot(2,1,2);
plot(t_sim, u_ideal_loop(2,:), 'r-', 'LineWidth', 2); 
hold on;
plot(t_sim, u_ro(2,:), 'm--', 'LineWidth', 1.5);
plot(t_sim, u_id(2,:), 'g:', 'LineWidth', 2.5);
hold off;
ylabel('u_2 [N.m]');
xlabel('Tempo [s]');
grid on;
legend('Ideal', 'Obs. Ordem Reduzida', 'Obs. Ordem Completa', 'Location', 'northeast');
title('Sinal de Controle u_2');