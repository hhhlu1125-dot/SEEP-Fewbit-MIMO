function [y, Delta, eta_b] = few_bit_quantizer(x, b, input_power)

    % Optimal Delta table
    Delta_table = [sqrt(8.0/pi), 0.9957, 0.5860, 0.3352, 0.1890, 0.1050];
    eta_b_table = [1.0 - 2.0 / pi, 0.1188, 0.0374, 0.0115];
    Delta = Delta_table(b) * sqrt(input_power);
    eta_b = eta_b_table(b);

    % 阈值
    tau = ((-2^(b-1) + 1):(2^(b-1) - 1)) * Delta;   % 1 x (2^b-1)

    % 展平
    x_vec = x(:);

    % 等价于 np.digitize(x, tau)，idx ∈ {0,1,...,2^b-1}
    idx = sum(x_vec > tau, 2);

    % 初始化
    y_vec = zeros(size(x_vec));

    % 左饱和
    y_vec(idx == 0) = -(2^b - 1) * Delta / 2;

    % 右饱和
    y_vec(idx == 2^b - 1) = (2^b - 1) * Delta / 2;

    % 中间区间
    mask = (idx > 0) & (idx < 2^b - 1);
    y_vec(mask) = tau(idx(mask)) + Delta / 2;

    % 恢复原尺寸
    y = reshape(y_vec, size(x));

end