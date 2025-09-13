% --- 1. 参数设置 ---

% 定义横坐标范围
x_min = 0;
x_max = 3600;

% 定义纵坐标波动范围
y_min = -0.009;
y_max = 0.009;

% 定义数据点的数量 (从0到3600，共3601个点)
num_points = x_max - x_min + 1;


% --- 2. 生成模拟数据 ---

% 创建横坐标向量
x = linspace(x_min, x_max, num_points);

% 生成一组基础随机数据，使其看起来不规则但又连续
% 使用randn生成正态分布的随机数，然后通过移动平均使其平滑，模拟更真实的波动
smoothing_factor = 15; % 可以调整这个值来改变曲线的“剧烈”程度，值越小越剧烈
raw_y = randn(1, num_points);
y_processed = movmean(raw_y, smoothing_factor);

% 将处理后的数据通过线性缩放，精确地映射到您指定的[y_min, y_max]区间
y_scaled = y_min + (y_processed - min(y_processed)) * (y_max - y_min) / (max(y_processed) - min(y_processed));


% --- 3. 绘图 ---

% 创建一个新的图形窗口
figure;

% --- 颜色选择 ---
% 请从下面三个选项中选择一个您想要的颜色，取消对应行的注释即可。
% 其他两行请确保是注释状态 (行首有 '%' 符号)。

% 选项 1: 蓝色
plot(x, y_scaled, 'b', 'LineWidth', 1.2, 'DisplayName', 'xchange');

% 选项 2: 红色
% plot(x, y_scaled, 'r', 'LineWidth', 1.2, 'DisplayName', 'xchange');

% 选项 3: 黄色
% plot(x, y_scaled, 'y', 'LineWidth', 1.2, 'DisplayName', 'xchange');


% --- 4. 美化图形 ---

% 设置坐标轴范围 (在y轴上下留出一点空白，让图形更美观)
xlim([x_min, x_max]);
ylim([y_min - 0.001, y_max + 0.001]);

% 添加标题和坐标轴标签
title('X方向定位结果图', 'FontSize', 14);
xlabel('历元', 'FontSize', 12);
ylabel('change（m)', 'FontSize', 12);

% 显示图例
legend('show');

% 添加网格线
grid on;

% 设置图形背景为白色
set(gcf, 'color', 'w');