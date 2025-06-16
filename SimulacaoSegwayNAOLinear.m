% Código para simulação do modelo NÃO linearizado do Segway

% Andre Hatano Feris - 11382644 
% Daniel Marins Silva D'Oliveira Martins - 12554571
% Felipe Almeida Ribeiro - 12566564 
% Lucas Junji Koreeda - 11347262

clear; clc; close all;

% Parâmetros
params.Mp  = 100;      % Massa do conjunto (chassi + pessoa) [kg]
params.Mr  = 7.5;       % Massa de cada roda [kg]
params.R   = 0.2;    % Raio da roda [m]
params.D   = 0.60;    % Distância entre rodas [m]
params.L   = 0.60;    % Distância vertical do CM ao eixo [m]
params.g   = 9.81;    % Aceleração da gravidade [m/s^2]
params.Jr  = 0.15;    % Momento de inércia da roda [kg*m^2]
params.Jpt = 36;   % Momento de inércia do corpo (pitch) [kg*m^2]
params.Jpy = 4.5;    % Momento de inércia do corpo (yaw) [kg*m^2]
params.kt  = 0.05;    % Constante de torque do motor [Nm/A]
params.ke  = 0.05;    % Constante eletromotriz do motor [V*s/rad]
params.Ra  = 1;     % Resistência do enrolamento do motor [Ohm]
params.N   = 24;      % Razão de redução do motor

params.k1 = (params.N * params.kt) / params.Ra;
params.k2 = (params.N^2 * params.kt * params.ke) / (params.R * params.Ra);

% Configurando a Simulação
% Tempo simulado
t_span = [0 25];
% Condições Iniciais [x_RM, x_dot, theta_p, theta_p_dot, theta_y, theta_y_dot]'
x0 = [0; 0; 0.3; 0; 0; 0];

% Simulando
[t, x] = ode45(@(t, x_vec) ode_segway_nao_linear(t, x_vec, params), t_span, x0); % simula


% Gerando os gráficos
% Extraindo cada estado da matriz de resultados 'x'
x_pos       = x(:,1);
x_vel       = x(:,2);
theta_p     = x(:,3);
theta_p_vel = x(:,4);
theta_y     = x(:,5);
theta_y_vel = x(:,6);

% Criando as figura para os gráficos
figure;
sgtitle('Simulação do Segway Não Linear');

% Gráfico do Ângulo de Pitch
subplot(2, 3, 1);
plot(t, rad2deg(theta_p), 'b-', 'LineWidth', 1.2);
title('Ângulo de Pitch');
%title('($\theta_P$)', 'Interpreter','latex');
xlabel('Tempo (s)');
ylabel('Ângulo (graus)');
grid on;
legend('$\theta_P$', 'Interpreter','latex');

% Gráfico da Posição Longitudinal
subplot(2, 3, 2);
plot(t, x_pos, 'r-', 'LineWidth', 1.2);
title('Posição Longitudinal');
xlabel('Tempo (s)');
ylabel('Posição (m)');
grid on;
legend('$x_{RM}$', 'Interpreter','latex');

% Gráfico da Velocidade Angular de Pitch
subplot(2, 3, 4);
plot(t, rad2deg(theta_p_vel), 'b-', 'LineWidth', 1.2);
title('Velocidade Angular de Pitch');
xlabel('Tempo (s)');
ylabel('Velocidade (graus/s)');
grid on;
legend('$\dot{\theta}_P$', 'Interpreter','latex');

% Gráfico da Velocidade Longitudinal
subplot(2, 3, 5);
plot(t, x_vel, 'r-', 'LineWidth', 1.2);
title('Velocidade Longitudinal');
xlabel('Tempo (s)');
ylabel('Velocidade (m/s)');
grid on;
legend('$\dot{x}_{RM}$', 'Interpreter','latex');

% Gráfico do Ângulo de Yaw
subplot(2, 3, 3);
plot(t, rad2deg(theta_y), 'g-', 'LineWidth', 1.2);
title('Ângulo de Yaw');
xlabel('Tempo (s)');
ylabel('Ângulo (graus)');
grid on;
legend('$\theta_y$', 'Interpreter','latex');

% Gráfico da Velocidade Angular de Yaw
subplot(2, 3, 6);
plot(t, rad2deg(theta_y_vel), 'g-', 'LineWidth', 1.2);
title('Velocidade Angular de Yaw');
xlabel('Tempo (s)');
ylabel('Velocidade (graus/s)');
grid on;
legend('$\dot{\theta}_y$', 'Interpreter','latex');


% Definindo as EDOs do sistema não linear
function dxdt = ode_segway_nao_linear(t, x, params) % <-- Adicionado 'params' aqui
    
    % Os parâmetros agora são passados explicitamente para esta função local.
    
    % Desempacotar o vetor de estados
    Mp  = params.Mp;
    Mr  = params.Mr;
    R   = params.R;
    D   = params.D;
    L   = params.L;
    g   = params.g;
    Jr  = params.Jr;
    Jpt = params.Jpt;
    Jpy = params.Jpy;
    k1  = params.k1; % Constante derivada
    k2  = params.k2; % Constante derivada

    x_rm        = x(1);
    x_dot       = x(2);
    theta_p     = x(3);
    theta_p_dot = x(4);
    theta_y     = x(5);
    theta_y_dot = x(6);
    
    % Cálculo das acelerações
    C_p = cos(theta_p);
    S_p = sin(theta_p);
    
    % Matriz de massa M(q)
    M11 = 2*(Mr + Jr/R^2) + Mp; 
    M12 = Mp * L * C_p;         
    M21 = Mp * L * C_p;
    M22 = Jpt + Mp * L^2;       
    M_matrix = [M11, M12; M21, M22];
    
    % Entradas de controle (malha aberta)
    uL = 0;
    uR = 0;
    
    % Vetor de forças F(q, q_dot)
    f1 = Mp * L * S_p * theta_p_dot^2 + (k1/R)*(uL + uR) - (2*k2/R^2)*x_dot;
    f2 = Mp * g * L * S_p + (k1)*(uL + uR) - (2*k2/R)*x_dot;
    F_vector = [f1; f2];
    
    % Resolvendo para as acelerações
    accel_pitch_trans = M_matrix \ F_vector;
    x_ddot       = accel_pitch_trans(1);
    theta_p_ddot = accel_pitch_trans(2);
    
    % Dinâmica de Yaw
    theta_y_ddot = (D/(2*R*(Jpy*R^2 + 2*Jr + 2*Mr*R^2))) * (k1*(uR - uL) - (k2*D/R)*theta_y_dot);
    
    % Montar o vetor derivada do espaço de estados
    dxdt = zeros(6,1);
    dxdt(1) = x_dot;
    dxdt(2) = x_ddot;
    dxdt(3) = theta_p_dot;
    dxdt(4) = theta_p_ddot;
    dxdt(5) = theta_y_dot;
    dxdt(6) = theta_y_ddot;
end

% Preparando a animação
h_fig_anim = figure;
set(h_fig_anim, 'Color', 'w'); % Fundo branco
set(h_fig_anim, 'DoubleBuffer', 'on'); % Melhora a suavidade da animação
axis equal; % Garante que as proporções dos eixos sejam iguais (círculos parecem círculos)
grid on;
hold on; % Permite que múltiplas coisas sejam plotadas no mesmo eixo sem apagar o anterior

% Definindo os limites dos eixos para a animação
x_min_anim = min(x_pos) - params.L - params.R; % Margem para o corpo inclinado + roda
x_max_anim = max(x_pos) + params.L + params.R; % Margem para o corpo inclinado + roda

y_max_anim = params.R + params.L + 0.1; % Altura máxima do CM + margem
y_mix_anim = params.R - params.L - 0.1; % Altura máxima do CM + margem

xlim([x_min_anim, x_max_anim]);
ylim([y_mix_anim, y_max_anim]);
xlabel('Posição Longitudinal (m)');
ylabel('Altura (m)');
title('Simulação não linear do modelo do Segway');

% Definindo os objetos gráficos
% Roda (círculo vazado)
theta_circulo = linspace(0, 2*pi, 100); % 100 pontos para um círculo suave

x_roda_circulo_initial = 0 + params.R * cos(theta_circulo); % A roda começa em x=0
y_roda_circulo_initial = params.R + params.R * sin(theta_circulo); % A roda fica a R de altura

% Plotando a roda
h_roda = plot(x_roda_circulo_initial, y_roda_circulo_initial, 'k-', 'LineWidth', 1.5);

% Ponto de conexão do pêndulo (centro da roda)
h_ponto_conexao = plot(0, params.R, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

% Corpo do Segway (linha do centro da roda ao CM)
% Inicialmente vertical, partindo do centro da roda
h_corpo = plot([0 0], [params.R params.R+params.L], 'b-', 'LineWidth', 3);

% Ponto do Centro de Massa (CM)
h_cm = plot(0, params.R+params.L, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

% Linha do chão (para referência)
plot([x_min_anim, x_max_anim], [0 0], 'k--', 'LineWidth', 0.5);

% Loop de animação -> calcula quadro por quadro
animation_speed_factor = 2; % Atualiza a animação a cada x iterações em t

for k = 1:animation_speed_factor:length(t)
    % Parâmetros do Segway no tempo atual
    current_x_pos = x_pos(k); % Posição longitudinal do centro da roda
    current_theta_p = theta_p(k); % Ângulo de Pitch do corpo

    % Calculo dos movimentos
    % Centro da Roda: A roda se move horizontalmente com x_pos
    roda_x_centro = current_x_pos;
    roda_y_centro = params.R; % O centro da roda está sempre a R acima do chão

    % Coordenadas do CM (Centro de Massa) do corpo
    % coordenadas da roda + angulo que L faz no momento
    cm_x = roda_x_centro + params.L * sin(current_theta_p); 
    cm_y = roda_y_centro + params.L * cos(current_theta_p);

    % Atualiza a posição da Roda (círculo)
    % Re-calcula os pontos do círculo transladados para a nova posição do centro da roda
    x_circulo_atualizado = roda_x_centro + params.R * cos(theta_circulo);
    y_circulo_atualizado = roda_y_centro + params.R * sin(theta_circulo);
    set(h_roda, 'XData', x_circulo_atualizado, 'YData', y_circulo_atualizado);

    % Atualiza o Ponto de Conexão (centro da roda)
    set(h_ponto_conexao, 'XData', roda_x_centro, 'YData', roda_y_centro);

    % Atualiza o corpo (linha do centro da roda ao CM)
    set(h_corpo, 'XData', [roda_x_centro, cm_x], 'YData', [roda_y_centro, cm_y]);

    % Atualiza a posição do CM
    set(h_cm, 'XData', cm_x, 'YData', cm_y);

    % Desenha o frame
    drawnow; % Força o MATLAB a atualizar o gráfico
end