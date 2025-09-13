function [Alpha_score, Alpha_pos, Convergence_curve] = A_GGWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj, AGGWO_params)
% -------------------------------------------------------------------------------------------------
% 函数定义
% 输入参数:
%   SearchAgents_no: 搜索智能体数量 (即种群大小 Pop_size)
%   Max_iter:        最大迭代次数
%   lb:              搜索空间的下界 (一个1 x dim的向量)
%   ub:              搜索空间的上界 (一个1 x dim的向量)
%   dim:             问题的维度
%   fobj:            目标函数句柄 (已经过预计算优化，直接传入即可)
%
% 输出参数:
%   Alpha_score:     找到的全局最优适应度值
%   Alpha_pos:       找到的全局最优解 (位置)
%   Convergence_curve: 记录每次迭代最优值的收敛曲线 (一个1 x Max_iter的向量)
% -------------------------------------------------------------------------------------------------

% --- 初始化全局最优记录变量 ---
% 沿用GWO的命名方式，用于记录迭代过程中的前三名。
% Alpha: 全局最优解
% Beta:  全局次优解 (在本算法中未直接使用，但保留了更新逻辑)
% Delta: 全局第三优解 (在本算法中未直接使用，但保留了更新逻辑)
Alpha_pos = zeros(1, dim); Alpha_score = inf;
Beta_pos  = zeros(1, dim); Beta_score  = inf;
Delta_pos = zeros(1, dim); Delta_score = inf;

% --- 初始化种群位置 ---
% 调用外部的 initialization 函数生成在[lb, ub]范围内的初始种群。
Positions = initialization(SearchAgents_no, dim, ub, lb);
% 【关键】针对GNSS整数模糊度问题，对初始种群进行取整。
Positions = round(Positions);

% --- 【核心改进】自适应DE (jDE) 参数初始化 ---
% 为种群中的每一个个体（每一只狼）都分配一套独立的控制参数 F 和 CR。
% F: 缩放因子 (Scaling Factor)，控制差分向量的扰动幅度。
% CR: 交叉概率 (Crossover Rate)，控制新个体从变异向量中继承基因的比例。
F_i = 0.5 * ones(SearchAgents_no, 1);  % 初始化所有个体的 F 值为0.5
CR_i = 0.8 * ones(SearchAgents_no, 1); % 初始化所有个体的 CR 值为0.8

% --- 初始化收敛曲线记录数组 ---
Convergence_curve = zeros(1, Max_iter);


% ================================ 主 循 环 开 始 ================================
for l = 1:Max_iter
    
    % --- 步骤 1: 评估当前种群并更新领头狼 ---
    % 在每次迭代开始时，计算当前整个种群所有个体的适应度值。
    % 这是后续精英选择、个体选择以及领头狼更新的基础。
    all_fitness = feval_all(fobj, Positions);
    
    % 在评估过程中，顺便更新Alpha, Beta, Delta三头领头狼的记录。
    % 尽管在我们的DE引擎中不直接使用Beta和Delta，但保留这一步可以
    % 方便未来进行其他混合策略的扩展。
    for i = 1:SearchAgents_no
        if all_fitness(i) < Alpha_score
            % 如果发现一个比当前Alpha更好的解，则更新整个最优解链条。
            Delta_score = Beta_score; Delta_pos = Beta_pos;
            Beta_score = Alpha_score; Beta_pos = Alpha_pos;
            Alpha_score = all_fitness(i); Alpha_pos = Positions(i, :);
        elseif all_fitness(i) < Beta_score
            Delta_score = Beta_score; Delta_pos = Beta_pos;
            Beta_score = all_fitness(i); Beta_pos = Positions(i, :);
        elseif all_fitness(i) < Delta_score
            Delta_score = all_fitness(i); Delta_pos = Positions(i, :);
        end
    end
    
    
    % --- 步骤 2: 【核心引擎】向量化的自适应差分进化 ---
    % 这是算法的核心，完全取代了标准GWO的位置更新公式。
    % 通过矩阵运算，一次性完成对整个种群的“变异”和“交叉”操作。
    
    % --- 2a. 向量化生成随机索引 r1, r2, r3 ---
    % 为了执行 "DE/rand/1" 变异，需要为每个个体 i 找到三个与之不同的
    % 随机个体 r1, r2, r3。这个辅助函数高效地完成了此任务。
    [r1, r2, r3] = get_distinct_indices(SearchAgents_no);
    
    % --- 2b. 向量化变异 (Vectorized Mutation) ---
    % 核心公式: mutant = base + F * (diff1 - diff2)
    % 我们将每个个体的 F_i 值复制扩展成与种群同维度的矩阵 F_matrix，
    % 从而可以通过矩阵的点乘运算，一次性为所有个体计算出变异向量。
    F_matrix = repmat(F_i, 1, dim);
    mutant_vectors = Positions(r1, :) + F_matrix .* (Positions(r2, :) - Positions(r3, :));
    % 【关键】再次对变异后的结果进行取整，确保在整数空间内操作。
    mutant_vectors = round(mutant_vectors);

    % --- 2c. 向量化交叉 (Vectorized Crossover) ---
    % 生成一个与种群同维度的[0,1]随机矩阵。
    rand_matrix = rand(SearchAgents_no, dim);
    % 将每个个体的 CR_i 值也扩展成CR_matrix。
    CR_matrix = repmat(CR_i, 1, dim);
    % 通过比较，生成一个逻辑掩码(mask)。mask中为true的位置，将从变异向量继承。
    mask = rand_matrix <= CR_matrix;
    % 为防止某一个体交叉后与原个体完全一样，我们强制每一行(每个个体)
    % 至少有一个维度(j_rand)必须进行交换。
    j_rand = randi(dim, SearchAgents_no, 1);
    % 使用线性索引 sub2ind 快速地将这些强制交换的位置在mask中设为true。
    mask(sub2ind(size(mask), 1:SearchAgents_no, j_rand')) = true;
    
    % 初始化试验向量为原种群
    trial_vectors = Positions;
    % 使用逻辑掩码，一次性将所有需要交叉的位置替换为变异向量中的对应值。
    trial_vectors(mask) = mutant_vectors(mask);
    
    % --- 2d. 向量化边界检查 (Vectorized Boundary Check) ---
    % 同样通过repmat扩展边界向量，并使用max/min函数一次性完成所有个体的边界处理。
    trial_vectors = max(trial_vectors, repmat(lb, SearchAgents_no, 1));
    trial_vectors = min(trial_vectors, repmat(ub, SearchAgents_no, 1));

    % --- 2e. 向量化选择 (Vectorized Selection) ---
    % 评估所有新生成的试验向量的适应度。
    trial_fitness = feval_all(fobj, trial_vectors);
    
    % 比较新旧个体的适应度，生成一个逻辑向量。
    % improvement_indices 中为true的位置，代表该个体在新一代中得到了改善。
    improvement_indices = trial_fitness < all_fitness;
    % 使用这个逻辑向量，一次性更新所有成功进化的个体。
    Positions(improvement_indices, :) = trial_vectors(improvement_indices, :);
    
    % --- 步骤 3: 【核心改进】自适应更新 F 和 CR 参数 ---
    % 这是jDE策略的关键。我们根据上一“选择”步骤的结果，来决定是否更新F和CR。
    % 对那些未能在本轮进化中产生更优解的个体 (improvement_indices为false)，
    % 我们认为它们当前的F和CR参数组合可能是无效的，因此需要进行重置。
    % `~improvement_indices` 找到了所有需要重置参数的个体的索引。
    failed_indices = ~improvement_indices;
    num_failed = sum(failed_indices); % 计算失败的个体数量
    
    % 为所有失败的个体重新生成一个在[0,1]范围内的随机CR值。
    CR_i(failed_indices) = rand(num_failed, 1);
    % 为所有失败的个体重新生成一个在[0.1, 1.0]范围内的随机F值。
    % 这种“推倒重来”的策略，鼓励算法在后续迭代中尝试全新的搜索行为。
    F_i(failed_indices) = 0.1 + rand(num_failed, 1) * 0.9;
    
    % 对于成功进化的个体，它们的F和CR参数被证明是有效的，因此予以保留，
    % 在下一代中继续使用。


    % --- 步骤 4: 【核心改进】无条件反向学习 (OBL) 策略 ---
    % 作为一个强大的全局探索和“跳出陷阱”的辅助手段。
    % 在每次迭代的最后，我们都对当前最优的Alpha狼进行一次反向探测。
    
    % 核心公式: Opposite_Point = Lower_Bound + Upper_Bound - Current_Point
    Opposite_Alpha_pos = round((lb + ub) - Alpha_pos);
    
    % 同样需要对反向点进行边界检查。
    Opposite_Alpha_pos = max(Opposite_Alpha_pos, lb);
    Opposite_Alpha_pos = min(Opposite_Alpha_pos, ub);
    
    % 计算反向点的适应度。
    opposite_fitness = fobj(Opposite_Alpha_pos);
    
    % 如果反向点被证明是更优的解...
    if opposite_fitness < Alpha_score
        % ...则更新全局最优解。
        Alpha_pos = Opposite_Alpha_pos;
        Alpha_score = opposite_fitness;
        
        % 【增强操作】同时，用这个新发现的优质解替换掉当前种群中最差的个体。
        % 这能加速优质基因在种群中的传播，提升整体种群质量。
        [~, worst_idx] = max(all_fitness);
        Positions(worst_idx, :) = Alpha_pos;
    end
    
    % --- 步骤 5: 记录与显示 ---
    % 将本次迭代结束时的全局最优值记录到收敛曲线中。
    Convergence_curve(l) = Alpha_score;
    
    % 为了避免在命令行中刷屏，我们设定每隔一定代数打印一次当前进度。
    if mod(l, 10) == 0 
        disp(['Iteration ', num2str(l), '/', num2str(Max_iter), ' - Best Fitness: ', num2str(Alpha_score)]);
    end
    
end
% ================================ 主 循 环 结 束 ================================

end


% ================================== 辅 助 函 数 ==================================

% --- 辅助函数1：向量化评估所有个体 ---
% 尽管MATLAB鼓励完全的向量化，但目标函数fobj本身通常是针对单个向量设计的。
% 这个函数通过一个简单的循环，将fobj应用于种群中的每一个体，是效率和
% 代码可读性之间的一个良好折中。
function fitness_values = feval_all(fobj, population)
    num_agents = size(population, 1);
    fitness_values = zeros(num_agents, 1);
    for i = 1:num_agents
        fitness_values(i) = fobj(population(i, :));
    end
end

% --- 辅助函数2：向量化获取不重复索引 ---
% 为DE/rand/1变异策略高效地生成所需的随机索引。
% 对于种群中的每个个体i，确保r1(i), r2(i), r3(i) 都不等于 i。
function [r1, r2, r3] = get_distinct_indices(pop_size)
    % 生成一个 1 到 pop_size 的序列。
    ns = 1:pop_size;
    % 预分配输出数组。
    r1 = zeros(pop_size, 1);
    r2 = zeros(pop_size, 1);
    r3 = zeros(pop_size, 1);
    
    for i = 1:pop_size
        % 对于每个个体i，生成一个不包含i的索引池。
        indices = ns(ns ~= i);
        % 从这个池中不重复地随机抽取3个索引。
        p = randperm(pop_size - 1, 3);
        % 分配给r1, r2, r3。
        r1(i) = indices(p(1));
        r2(i) = indices(p(2));
        r3(i) = indices(p(3));
    end
end