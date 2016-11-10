function [ feature_in_1_cam_comb, cluster_centers ] = feature_merge( feature, i_cam, numClusters, cluster_centers)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% 先使用简单的相加合成方法
% i_cam = 1; % 只使用这一个视角下的相机
feature_in_1_cam = feature(:,i_cam);
feature_in_1_cam_comb = cell(size(feature_in_1_cam));

%% BOW词袋模型

if ~exist('cluster_centers', 'var')
    tmp_feature_in_1_cam = feature_in_1_cam;
    for seq=1:numel(tmp_feature_in_1_cam)
        this_fe = reshape(tmp_feature_in_1_cam{seq}, 1, []);
        this_fe = this_fe(~isemptycell(this_fe));
        this_fe = cell2mat(this_fe);
        tmp_feature_in_1_cam{seq} = this_fe;
    end
    % 合成一个128*N维向量，N为所有描述子的个数（共有14w特征）
    all_feature = cell2mat( reshape(tmp_feature_in_1_cam, 1, []) );

    % numClusters = 1000;% 词袋个数，不能太大，否则聚类会outofmemory
    % idx = kmeans(im2double(all_feature‘), N_bow); % 点数目太多，Kmean很慢（不用这个）
    % vl_feat中的kmeans很好用
    tic;
    fprintf('doing K-means with %d clusters...\n', numClusters);
    [cluster_centers, ~] = vl_kmeans(im2single(all_feature), numClusters, 'Initialization', 'plusplus');
    toc;
end

%% 特征融合方法
for seq=1:size(feature_in_1_cam,1)
    for i_p=1:size(feature_in_1_cam{seq},1)
        % 某个人物p的图像序列
        p_frames = feature_in_1_cam{seq}(i_p,:);
        % 后面有一些空矩阵，需要去除
        p_frames = p_frames(~isemptycell(p_frames));
        if isempty(p_frames)
            continue
        end
        
        tmp = im2single(cell2mat(p_frames));
        [~, assign] = min(vl_alldist(tmp, cluster_centers), [], 2);
        combine_fe = zeros(numClusters, 1);
        for x=assign' % 按个数增加相应位置
            combine_fe(x) = combine_fe(x) + 1;
        end
        % 必须要进行归一化
        combine_fe = combine_fe/sum(combine_fe);

        % 对这个人的序列进行特征融合
        % ----------------------------------------- %
        if 0
            % 简单相加、取均值合成combine_fe
            combine_fe = zeros(size(p_frames(1)));
            for i_f=1:p_frames_num
                combine_fe = combine_fe + p_frames{i_f};
            end
            % 可以显示一下合成后的图片
            combine_fe = combine_fe/p_frames_num; % imshow(combine_fe)
        end
        
        % ----------------------------------------- %
        % 保存在feature_in_1_cam_comb中
        feature_in_1_cam_comb{seq}{i_p,1} = combine_fe;
    end
end


end

