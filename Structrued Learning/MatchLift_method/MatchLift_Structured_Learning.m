
% GMMCP的结构化学习框架
clear; close all

diary('learning_log.txt')
diary on 

%% 1、特征提取
dataaddr = 'E:\人脸识别数据集\ChokePoint\cropped_face\P1E_cropped_face\P1E\';
if exist('feature.mat', 'file')
    load('feature.mat');
else
    feature = get_feature_wrapper(dataaddr);
    save('feature.mat', 'feature');
end

%% 2、特征合成

% 将一个人的视频序列合成一个特征
i_cam = 1;
[ feature_in_1_cam_comb ] = feature_merge( feature, i_cam );
% 滤掉其中的nan，并获得人脸标签
[ feature_in_1_cam_comb, labels ] = filter_nan_and_get_label( feature_in_1_cam_comb );
n_p_all = cellfun(@(x) size(x,1), feature_in_1_cam_comb,'un',1)';

% 设定要使用的seq数目
n_seq = 3;
n_person_in_seq = n_p_all(1:n_seq); % 每个seq中人数
n_person_cum = [0, cumsum(n_person_in_seq)]; % 人数的累加（用与矩阵分块）
labels = labels(1:n_seq);
% 获取gt的连接方式
[ gt_connect, truth_num ] = get_gt_match_result( labels );

% 对选定的seq进行pca降维
im_mat = cell2mat(feature_in_1_cam_comb(1:n_seq));
[X, mean_face, COEFF] = pca_eig( im_mat, 0.975 ); % imshow(reshape(mean_face,96,96))
feature_in_1_cam_comb_pca = mat2cell(X', n_person_in_seq, size(X',2));

fe_to_use = feature_in_1_cam_comb_pca(1:n_seq);

%% 使用朴素匹配的结果来估算m
% 计算矩阵度量矩阵W （越大月相似）
W  = cal_W( fe_to_use, n_person_in_seq );
Kval = mat2bolckmat(W, n_person_cum);
pre_connect = naive_match( Kval );
% 利用朴素匹配的结果来估算m
[ m_hat, X_in_mat ] = estimate_m_by_inputX( pre_connect, n_person_in_seq );

%% 分配变量、计算约束、phi
N = 1; % 暂时只有一个样本（4个场景）

vars_ADN_F = cal_MatchLift_constraints( X_in_mat, m_hat, n_person_in_seq );
[ phi_x_z, phi_x_zstar, loss ] = MatchLift_cal_phi_and_loss( vars_ADN_F, fe_to_use, gt_connect, truth_num, n_person_cum );

%% 分配变量空间，定义参数
rng(0)
w = (rand(size(phi_x_z))-0.5)*10;
w = zeros(size(phi_x_z));

iter = 500;
W = cell(iter,1); % W 存放综合权值w
Wavg = cell(iter,1);
Wi = cell(iter,N); % Wi存放样本权值w
Wavg{1} = w;
W{1} = w; % 全体样本的W需要设定初值w
for i=1:N
    Wi{1,i} = w; % Wi 存放每次循环中特定样本更新后的Wi
end
L = zeros(iter,1); % L 存放综合L
Li = zeros(iter,N);
for i=1:N
    Li(1,i) = 0;
end

lambda = 1e0;
pre_connect = cell(iter,1); % 每轮预测结果
gamma = zeros(iter,1); % 步长gamma
time_consume = zeros(iter,1); % 记录每次循环所用的时间
aver_loss = zeros(iter,1); % 记录每一轮中样本损失函数均值
linesearch = 0; % 是否使用线搜索策略
t = 0;


%% 开始循环
while t<iter
    t = t + 1;
    tic;
    disp('  ==========================');
    disp(['  开始第 ', num2str(t), ' 轮循环...']);
    % 选中第ind个样本作为训练样本
    ind = randi(N); 

    % 怀疑速度慢是由dot(w,phi)引起？？？因此不使用这种形式，而临时计算<w,feature>，再组合为目标函数
    [ pre_connect{t}, phi_x_zhat, delta_zstar_zhat ] = CXSL_cal_MatchLift_object( Wavg{t}, vars_ADN_F, phi_x_z, loss, n_person_cum );
    disp('      更新权向量 w...');
    
    U_x_zstar_zhat = phi_x_zstar - phi_x_zhat;
    Ws = U_x_zstar_zhat/(lambda*N);
    ls = delta_zstar_zhat/N;
    aver_loss(t) = delta_zstar_zhat;
    fprintf('      当前样本损失函数△(z*,z^):\t%f\n', delta_zstar_zhat);
    
    % 计算gap，gap的值随样本、lambda都会变化，无法确定下来，因此还是用loss做gap比较合适
    if linesearch
        tmp = lambda*(Wi{t,ind}- Ws)'*W{t}- Li(t,ind)+ ls;
        gamma(t) = tmp/(lambda*norm(Wi{t,ind}- Ws)^2);
        gamma(t) = max([0, min([gamma(t),1])]);
    else
        gamma(t) = 2*N/(2*N + t-1);
    end

    % 更新 wi和Li，将更新后的w保存在 W{t+1,ind}中
    % 对于Wi来说，由于上一次不一定选得是i，因此Wi{t,ind}可能为空的，会导致更新的时候似乎有问题？
    % 2015.10.23 发现其实没有被选中的样本是直接将Wi带入了下轮！
    % 直接将所有样本的Wi带入下一轮
    for ii=1:N
        Wi{t+1,ii} = Wi{t,ii}; 
        Li(t+1,ii) = Li(t,ii);
    end
    % 再更新单个样本的Wi
    Wi{t+1,ind} = (1- gamma(t))*Wi{t,ind}+ gamma(t)*Ws; 
    Li(t+1,ind) = (1- gamma(t))*Li(t,ind)+ gamma(t)*ls;

    % 更新 w和L，将更新后的w保存在 W{t+1,N+1}中
    W{t+1} = W{t} + Wi{t+1,ind} - Wi{t,ind};
    L(t+1) = L(t) + Li(t+1,ind) - Li(t,ind); % 计算L似乎没什么用？
    
    % 更新Wavg
    Wavg{t+1} = (t-1)/(t+1)*Wavg{t} + 2/(t+1)*W{t+1};
    
    % 记录时间
    time_consume(t) = toc;
    fprintf('      时间花费:\t%1.2f s\n', time_consume(t)); 
end

%% 后续处理
% 绘制损失曲线
plot(aver_loss, '-*');
t_best = find(aver_loss==min(aver_loss));
[ ACC, TP_cell ] = cal_precision( pre_connect{t_best}, gt_connect, truth_num );


save('aver_loss.mat', 'aver_loss')
diary off;
system('shutdown /s');


