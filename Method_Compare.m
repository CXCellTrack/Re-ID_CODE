clear, close all


% 比较在不学习/pspsp半监督学习下，GMMCP/MATCHLIST/NAIVE方法的性能
% ------------------------- Main file ------------------------- % 

%% 1、特征提取
% 特征在 GMMCP_Structured_Learning.m 中计算
load('feature.mat');

%% 2、特征合成

i_cam = 1;
numClusters = 1000;
load('cluster_centers-1000.mat');

% 选择测试集特征
feature_test = feature.P1L;
[ feature_in_1_cam_comb, ~ ] = feature_merge( feature_test, i_cam, numClusters, cluster_centers );
% 滤掉其中的nan，并获得人脸标签
[ feature_in_1_cam_comb, labels ] = filter_nan_and_get_label( feature_in_1_cam_comb );
n_p_all = cellfun(@(x) size(x,1), feature_in_1_cam_comb,'un',1)';

% 设定要使用的seq数目
n_seq = 3;
fe_to_use = feature_in_1_cam_comb(1:n_seq);
n_person_in_seq = n_p_all(1:n_seq); % 每个seq中人数
n_person_cum = [0, cumsum(n_person_in_seq)]; % 人数的累加（用与矩阵分块）
labels = labels(1:n_seq);
[ gt_connect, truth_num ] = get_gt_match_result( labels );
% 计算矩阵度量矩阵W （越大月相似）
W  = cal_W( fe_to_use, n_person_in_seq );

% 是否在最后的GMMCP方法中设置每个clique为一个人（若设1会对PCPSP的约束产生影响）
global each_person_a_clique;
each_person_a_clique = 0; % 从效果上看GMMCP中每个clique只有一个人是不好的！

%% 3、使用PCPSP方法求解核内积Kval（作为距离度量）
if 1
    % 直接使用欧氏距离作为Kval
    Kval = mat2bolckmat(W, n_person_cum);
else
    max_must_link = 0;
    max_cannot_link = 0; % 必须至少有一个mustlink和cannotlink，否则kval全为1
    [Kval, seq_mustlink, seq_cannotlink] = PCPSP( W, n_person_in_seq, labels, max_must_link, max_cannot_link );
end

if each_person_a_clique
    W = log(W./(1-W)); % 将W变换为有正有负，否则会变为全链接
%     W(W==-Inf) = 0; % W对角线设为0
    Kval = num2cell(W);
end

%% 4、使用GMMCP方法进行轨迹连接

% 为了便于运行不同数目link下的结果，将计算约束与最终求解分开
% [vars, ADN] = GMMCP(Kval, n_person_in_seq);
vars_ADN_F = cal_GMMCP_constraints( n_person_in_seq );  

%% 5、求解GMMCP并进行精度计算
Cd = 0.00; % Cd为对ADN的赋值，越大说明越鼓励ADN的出现（each_person_a_clique=1时应将Cd增大）
[ pre_connect, ADN_value ] = cal_GMMCP_object( Kval, vars_ADN_F, Cd );
% 计算每2个clique匹配的精度
[ ACC, TP_cell ] = cal_precision( pre_connect, gt_connect, truth_num );

%% 6、朴素求解法：利用权值22进行匹配（匈牙利算法等）
pre_connect = naive_match( Kval );
% 简单22匹配的效果比全局的要差一些
[ ACC, TP_cell ] = cal_precision( pre_connect, gt_connect, truth_num );

%% 7、使用matchlift（5和7都是全局匹配的方法，应该择1使用）
% 假设 vars_value 为输入矩阵Xin，每个场景作为一个对象，2个场景的匹配关系作为一个map
% 构成一个分块矩阵
% 以距离度量矩阵W作为输入
[ X_in_mat, X_out_val_int ] = MatchLift( Kval, n_person_in_seq );
X_out = mat2bolckmat( X_out_val_int, n_person_cum );
% 有效果！
[ ACC, TP_cell ] = cal_precision( X_out, gt_connect, truth_num );
















