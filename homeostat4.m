clearvars();

% パラメータの設定
h = 1;
j = 0.5;
a11 = -1; a12 = 0.5; a13 = 0.2; a14 = 0;
a21 = 0.3; a22 = -2; a23 = 0; a24 = 0.1;
a31 = 0; a32 = 0.1; a33 = -1.5; a34 = 0.2;
a41 = 0.2; a42 = 0; a43 = 0.3; a44 = -1;

% 係数行列 A の定義
A = [a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44];

% 収束する例
A1 = [-1 -0.2  0.1  0.3;
       0.2 -2 -0.1  0.2;
      -0.1  0.1 -1.5 -0.3;
       0.3  0  0.2 -1];

A2 = [-0.5  0.1 -0.2  0.1;
         0.2 -1 -0.1  0.2;
        -0.1  0.2 -0.8 -0.2;
         0.1 -0.1  0.1 -1.2];

A3 = [-1.2  0.1 -0.1  0.2;
          0.1 -1.1  0.2 -0.1;
          0.2 -0.2 -0.9  0.1;
         -0.1  0.1  0.1 -1.3];

A4 = [-0.8  0.2 -0.1  0.1;
          0.1 -1.5  0.1 -0.2;
         -0.1  0.1 -0.7  0.2;
          0.2 -0.1 -0.2 -1.1];

A = A2;

% 時間範囲の設定
time_span = [0 20];

% 初期値の設定 (収束する例)
x0 = [1; 0; 0; 0];
x0_dot = [0; 0; 0; 0];
initial_conditions = [x0; x0_dot];

% 微分方程式を定義する関数
function dx_dt = coupled_system(t, x, A, j)
    x_1 = x(1); x_2 = x(2); x_3 = x(3); x_4 = x(4);
    x_1_dot = x(5); x_2_dot = x(6); x_3_dot = x(7); x_4_dot = x(8);
    dx_dt = zeros(8, 1);
    dx_dt(1) = x_1_dot;
    dx_dt(2) = x_2_dot;
    dx_dt(3) = x_3_dot;
    dx_dt(4) = x_4_dot;
    dx_dt(5) = A(1,1)*x_1 + A(1,2)*x_2 + A(1,3)*x_3 + A(1,4)*x_4 - j*x_1_dot;
    dx_dt(6) = A(2,1)*x_1 + A(2,2)*x_2 + A(2,3)*x_3 + A(2,4)*x_4 - j*x_2_dot;
    dx_dt(7) = A(3,1)*x_1 + A(3,2)*x_2 + A(3,3)*x_3 + A(3,4)*x_4 - j*x_3_dot;
    dx_dt(8) = A(4,1)*x_1 + A(4,2)*x_2 + A(4,3)*x_3 + A(4,4)*x_4 - j*x_4_dot;
end

% 時間発展をode45で計算
[t, x_solution] = ode45(@(t, x) coupled_system(t, x, A, j), time_span, initial_conditions);

% flow fieldの作成 (x1-x2平面への射影)
x1_range = linspace(-1.2, 1.2, 20);
x2_range = linspace(-1.2, 1.2, 20);
dx1_dt = zeros(length(x1_range), length(x2_range));
dx2_dt = zeros(length(x1_range), length(x2_range));

for i = 1:length(x1_range)
    for j = 1:length(x2_range)
         x = [x1_range(i); x2_range(j); 0; 0; 0; 0; 0; 0];
         dx = coupled_system(0, x, A, j);
         dx1_dt(i,j) = dx(5);
         dx2_dt(i,j) = dx(6);
    end
end

% flow fieldの作成 (x3-x4平面への射影)
x3_range = linspace(-2, 2, 20);
x4_range = linspace(-2, 2, 20);
dx3_dt = zeros(length(x3_range), length(x4_range));
dx4_dt = zeros(length(x3_range), length(x4_range));
for i = 1:length(x3_range)
    for j = 1:length(x4_range)
        x = [0; 0; x3_range(i); x4_range(j); 0; 0; 0; 0];
        dx = coupled_system(0, x, A, j);
        dx3_dt(i,j) = dx(7);
        dx4_dt(i,j) = dx(8);
    end
end

% グラフ描画
fig = figure();
fig.PaperOrientation = 'landscape';
fig.PaperPosition = [0.6345 0.6345 28.4310 19.7310];

% x1-x2 平面への射影
subplot(2, 2, 1);
quiver(x1_range, x2_range, dx1_dt', dx2_dt', 'AutoScaleFactor', 1);
hold on;
plot(x_solution(:, 1), x_solution(:, 2), 'b-', 'LineWidth', 1, 'DisplayName', 'Trajectory');
plot(x0(1), x0(2), 'bo', 'MarkerFaceColor', 'b')
xlabel('x1');
ylabel('x2');
title('Flow Field (x1-x2)');
grid on;
axis square
axis([-1.2 1.2 -1.2 1.2]);
hold off;

% x3-x4 平面への射影
subplot(2, 2, 3);
quiver(x3_range, x4_range, dx3_dt', dx4_dt', 'AutoScaleFactor', 1);
hold on;
plot(x_solution(:, 3), x_solution(:, 4), 'b-', 'LineWidth', 1, 'DisplayName', 'Trajectory');
plot(x0(3), x0(4), 'bo', 'MarkerFaceColor', 'b')
xlabel('x3');
ylabel('x4');
title('Flow Field (x3-x4)');
grid on;
axis square
axis([-1.2 1.2 -1.2 1.2]);
hold off;

%% 時間発展のグラフ
subplot(4, 2, 2);
hold on
plot([-5 0 0], [0 0 x0(1)], 'k', 'LineWidth', 1)
plot(t, x_solution(:, 1), 'k-', 'LineWidth', 1, 'DisplayName', 'x1');
axis([-Inf Inf -1.2 1.2])
subplot(4, 2, 4);
hold on
plot([-5 0 0], [0 0 x0(2)], 'k', 'LineWidth', 1)
plot(t, x_solution(:, 2), 'k-', 'LineWidth', 1, 'DisplayName', 'x2');
axis([-Inf Inf -1.2 1.2])
subplot(4, 2, 6);
hold on
plot([-5 0 0], [0 0 x0(3)], 'k', 'LineWidth', 1)
plot(t, x_solution(:, 3), 'k-', 'LineWidth', 1, 'DisplayName', 'x3');
axis([-Inf Inf -1.2 1.2])
subplot(4, 2, 8);
hold on
plot([-5 0 0], [0 0 x0(4)], 'k', 'LineWidth', 1)
plot(t, x_solution(:, 4), 'k-', 'LineWidth', 1, 'DisplayName', 'x4');
axis([-Inf Inf -1.2 1.2])
