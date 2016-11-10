
% GMMCP的结构化学习框架
clear; close all

%% 1、特征提取
baseaddr = 'E:\人脸识别数据集\ChokePoint\cropped_face\';
savename = 'feature.mat';
if exist(savename, 'file')
    load(savename);
else
    feature = struct;
    is_front = struct;
    for ds={'P1E','P1L','P2E','P2L'}
        disp(['process daatset ', ds{1}]);
        dataaddr = [baseaddr, ds{1}, '_cropped_face\', ds{1}, '\'];
        feature = setfield( feature, ds{1}, get_feature_wrapper(dataaddr) );
        is_front = setfield( is_front, ds{1}, check_if_front_face(dataaddr) );
    end
    save(savename, 'feature', 'is_front');
end

% 去掉非正面的脸
feature = remove_non_front_face(feature, is_front);

%% 2、特征合成

% 合并好几个场景下的特征
feature_train = {feature.P1E, feature.P1L};
f_train_final = make_f_train_final(feature_train);

savename = 'cluster_centers-1000.mat';
i_cam = 1;
numClusters = 1000;
if exist(savename, 'file')
    load(savename);
else
    [ ~, cluster_centers ] = feature_merge( f_train_final, i_cam, numClusters );
    save(savename, 'cluster_centers');
end

[ feature_in_1_cam_comb, ~ ] = feature_merge( f_train_final, i_cam, numClusters, cluster_centers );
% 滤掉其中的nan，并获得人脸标签
[ feature_in_1_cam_comb, labels ] = filter_nan_and_get_label( feature_in_1_cam_comb );
% 每个seq中的人数
n_p_all = cellfun(@(x) size(x,1), feature_in_1_cam_comb,'un',1)';
% 设定要使用的seq数目
n_seq = 3;
n_person_in_seq = n_p_all(1:n_seq); % 每个seq中人数
n_person_cum = [0, cumsum(n_person_in_seq)]; % 人数的累加（用与矩阵分块）
labels = labels(1:n_seq);

% ----------------------------------------------------------------------- %
% 对选定的seq进行pca降维
if 0
    im_mat = cell2mat(feature_in_1_cam_comb(1:n_seq));
    [X, mean_face, COEFF] = pca_eig( im_mat, 0.975 ); % imshow(reshape(mean_face,96,96))
    feature_in_1_cam_comb_pca = mat2cell(X', n_person_in_seq, size(X',2));
    fe_to_use = feature_in_1_cam_comb_pca(1:n_seq);
else
    fe_to_use = feature_in_1_cam_comb(1:n_seq);
end
    
%% 获取gt的连接方式
% 设置样本数目，规定样本切分方式
N = 5;
nn = cell(1,n_seq);
nn_cum = cell(1,n_seq);
for ii=1:n_seq
    % 将每个n_person_in_seq分为N段
    nn{ii} = [];
    cut = n_person_in_seq(ii);
    for i_N=1:N-1
        nn{ii} = [nn{ii}, floor(cut/N)];
    end
    nn{ii} = [nn{ii}, cut-sum(nn{ii})];
    nn_cum{ii} = [0, cumsum(nn{ii})];
end

% 分配各个样本对应的变量
N_p_in_seq = cell(N,1);
N_fe_to_use = cell(N,1);
N_labels = cell(N,1);
N_gt_connect = cell(N,1);
N_truth_num = cell(N,1);
for ii=1:N
    for jj=1:n_seq
        N_p_in_seq{ii}(jj) = nn{jj}(ii);
        fe_s = nn_cum{jj}(ii)+1;
        fe_e = nn_cum{jj}(ii+1);
        N_fe_to_use{ii}{jj,1} = fe_to_use{jj}(fe_s:fe_e,:);
        N_labels{ii}{jj,1} = labels{jj}(fe_s:fe_e,:);
    end
    [ N_gt_connect{ii}, N_truth_num{ii} ] = get_gt_match_result( N_labels{ii} );
end

%% 分配变量、计算约束、phi

N_vars_ADN_F = cell(N,1);
N_phi_x_z = cell(N,1);
N_phi_x_zstar = cell(N,1);
N_loss = cell(N,1);

for ii=1:N
    fprintf('计算样本 %d 的约束...\n', ii);
    N_vars_ADN_F{ii} = cal_GMMCP_constraints( N_p_in_seq{ii} );  
    [ N_phi_x_z{ii}, N_phi_x_zstar{ii}, N_loss{ii} ] = ...
    CXSL_cal_phi_and_loss( N_vars_ADN_F{ii}, N_fe_to_use{ii}, N_gt_connect{ii}, N_truth_num{ii} );
end

%% 分配变量空间，定义参数
rng(0)
iter = 200;
w = zeros(size(N_phi_x_z{1}));
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

lambda = 1e-3;
gamma = zeros(iter,1); % 步长gamma
time_consume = zeros(iter,1); % 记录每次循环所用的时间
sample_loss = zeros(iter,N); % 记录每一轮中样本损失函数均值
Cd = 0;
linesearch = 0; % 是否使用线搜索策略
t = 0;
random = 0;
ind = 0;

%% 开始循环
while t<iter
    
    %% 求解目标函数
    t = t + 1;
    tic;
    disp('  ==========================');
    disp(['  开始第 ', num2str(t), ' 轮循环...']);
    % 选中第ind个样本作为训练样本
    if random
        ind = randi(N);
    else
        ind = ind + 1;
        ind(ind==(N+1)) = 1;
    end
        
    disp(['  选择样本 ', num2str(ind)]);
    % 怀疑速度慢是由dot(w,phi)引起？？？因此不使用这种形式，而临时计算<w,feature>，再组合为目标函数
    [ pre_connect, ADN_value, phi_x_zhat, delta_zstar_zhat ] = ...
        CXSL_cal_GMMCP_object( Wavg{t}, N_vars_ADN_F{ind}, Cd, N_phi_x_z{ind}, N_loss{ind} );
    
    %% 更新w过程
    disp('      更新权向量 w...'); 
    U_x_zstar_zhat = N_phi_x_zstar{ind} - phi_x_zhat;
    Ws = U_x_zstar_zhat/(lambda*N);
    ls = delta_zstar_zhat/N;
    sample_loss(t, ind) = delta_zstar_zhat;
    fprintf('      当前样本损失函数△(z*,z^):\t%f\n', delta_zstar_zhat);
    
    % 计算gap，gap的值随样本、lambda都会变化，无法确定下来，因此还是用loss做gap比较合适
    if linesearch
        tmp = lambda*(Wi{t,ind}- Ws)'*W{t}- Li(t,ind)+ ls;
        gamma(t) = tmp/(lambda*norm(Wi{t,ind}- Ws)^2);
        gamma(t) = max([0, min([gamma(t),1])]);
    else
        gamma(t) = 2*N/(2*N + t-1);
    end

    %% 更新 wi和Li，将更新后的w保存在 W{t+1,ind}中
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
aver_loss = sum(sample_loss,2);
plot(aver_loss, '-*');
fprintf('min loss = %f\n', min(aver_loss));
t_best = find(aver_loss==min(aver_loss));
t_best
w_best = Wavg{t_best};






