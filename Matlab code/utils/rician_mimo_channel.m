function H = rician_mimo_channel(R, T, K)
    % R: 接收天线数
    % T: 发射天线数
    % K: 莱斯K因子（视距分量功率与散射分量功率的比）

    % 1. 生成视距（LoS）分量（假设为全1或定向排列）
    H_los = ones(R, T);  % 可替换为 exp(1j*theta) 等定向模型

    % 2. 生成瑞利散射分量（复高斯分布）
    H_nlos = (randn(R, T) + 1j * randn(R, T)) / sqrt(2);

    % 3. 计算莱斯通道矩阵
    H = sqrt(K / (K + 1)) * H_los + sqrt(1 / (K + 1)) * H_nlos;
end
