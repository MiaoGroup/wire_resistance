close all;
% resistivity 
% Cu
rho = 1.7e-8;
% wire parameters
% area
S = 5e-6*200e-9;
% length
L = 5e-6;

% 线路（相邻两个格点）电阻
% wire resistance
R0 = rho*L/S;
% test, 10 has huge effect
R0 =  10;

% 阵列大小，会被测试用例覆盖
% array N X N size
N = 32;

% 随机生成
% random memoristor resistance and input vector
R = round((rand(N, N))*36000) + 4000
V_input = rand(1, N)*20;

% 此处是测试用
% test example
% R = [7000 3000 3000 3000; 4000 6000 15000 15000; 9000 6000 6000 6000; 9000 6000 6000 6000];
% V_input = [3 6 7 7];
% [N, ~] = size(R);


[I_actual, conductance] = calculate(R, R0, N, V_input);

% 未调节时的电流输出
% I when no tune
I_no_tune = I_actual;
% 
R_equivalent_current = 1./conductance*1000;
R_equivalent_no_tune = R_equivalent_current

% 训练次数
% train times
epoch = 80;

% 学习速率
% learning rate
eta = 0.4;

% ideal output
I_ideal = V_input*(1./R)*1000;   

% 记录每次循环的电流偏移情况
% record  deviation in each circle
I_deviation = zeros(epoch, N);
I_deviation(1, :) = deviation(I_actual, I_ideal);

R_resistor_tune = R;
for i=2:epoch
    delta_R =  R_equivalent_current - R;
    R_resistor_tune = R_resistor_tune - eta*delta_R;
    % 随机输入
    % randome input
    V_input = rand(1, N)*20;
    I_ideal = V_input*(1./R)*1000;
    [I_actual, conductance] = calculate(R_resistor_tune, R0, N, V_input);
    R_equivalent_current = 1./conductance*1000;
    % 随机输入情况下的偏移
    % deviation for random input
    I_deviation(i, :) = deviation(I_actual, I_ideal);
end
plot(1:epoch, I_deviation);
xlim([1 epoch]);
title('mean deviation by tune conductance');
xlabel('times');
ylabel('deviation');
R_equivalent_current

function [A] = coefficient_matrix(R, R0, N)
% 格点总数
% all memoristor
nodenum = N*N;
% 系数矩阵
% coefficient matrix for A*I = V
A = sparse(nodenum, nodenum);
temp = zeros(N, nodenum);

% 系数矩阵第k行
% the kth row of the coefficient matrix
k = 0;
m = 0;
for i =1:N
    for j =1:N
        k = k + 1;
        m = m + 1;
        temp(m,:) = find_coefficient(R, R0, i, j, N);
        if m == N
            A(k-N+1:k, :) = sparse(temp);
        m = 0;
        end
%         A(k,:) = find_coefficient(R, R0, i, j, N);
    end
end
end

function [a_k] = find_coefficient(R, R0, i0, j0, N)
a_ij = zeros(N, N);
% 计算阵列中这一行上的电流Iij对应的系数
% calculate coefficient in selected row
for j=1:N
    if j < j0
        a_ij(i0, j) = R0*j;
    end    
end
a_ij(i0, j0+1:end) = R0*j0;
a_ij(i0, j0) = (N + 1 - i0 + j0)*R0 + R(i0, j0);

% 计算选中的列上的
% calculate column
for i=(i0+1):N
    a_ij(i, j0) = (N + 1 - i)*R0;
end
a_ij(1:i0-1, j0) = (N + 1 - i0)*R0;

% 降维为一列，按行展开，所以要转置一下
% need transpose, because of expanding by row
a_k = reshape(a_ij', [1, N*N]);
end

function [matrix] = find_dot_matrix(A, N)
matrix = zeros(N, N);
k = 0;
for i=1:N:N*N
    V = zeros(N*N, 1);
    V(i:i+N - 1) = 1;
    I = (A \ V)*1000;
    % reshape I matrix for analysising
    I = reshape(I, [N, N])';
    % output
    k = k + 1;
    matrix(k, :) = sum(I, 1);
end
end

function [I_actual, matrix] = calculate(R, R0, N, V_input)
    % generate V verctor from orginal V_INPUT
    V = zeros(N*N, 1);
    for i=1:N:N*N
        V(i:i+N - 1) = V_input(floor(i / N + 1));
    end
    % calculate A*I = V
    A = coefficient_matrix(R, R0, N);
    I = (A \ V)*1000;
    % reshape I matrix for analysising
    I = reshape(I, [N, N])';
    % output
    I_actual = sum(I);
    matrix = find_dot_matrix(A, N);
end

function [value] = deviation(I_acutal, I_ideal)
I_ideal_mean = mean(I_ideal);
value = abs((I_acutal - I_ideal)./I_ideal_mean);
end
