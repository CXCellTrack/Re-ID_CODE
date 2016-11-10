function [ feature_in_1_cam_comb, cluster_centers ] = feature_merge( feature, i_cam, numClusters, cluster_centers)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% ��ʹ�ü򵥵���Ӻϳɷ���
% i_cam = 1; % ֻʹ����һ���ӽ��µ����
feature_in_1_cam = feature(:,i_cam);
feature_in_1_cam_comb = cell(size(feature_in_1_cam));

%% BOW�ʴ�ģ��

if ~exist('cluster_centers', 'var')
    tmp_feature_in_1_cam = feature_in_1_cam;
    for seq=1:numel(tmp_feature_in_1_cam)
        this_fe = reshape(tmp_feature_in_1_cam{seq}, 1, []);
        this_fe = this_fe(~isemptycell(this_fe));
        this_fe = cell2mat(this_fe);
        tmp_feature_in_1_cam{seq} = this_fe;
    end
    % �ϳ�һ��128*Nά������NΪ���������ӵĸ���������14w������
    all_feature = cell2mat( reshape(tmp_feature_in_1_cam, 1, []) );

    % numClusters = 1000;% �ʴ�����������̫�󣬷�������outofmemory
    % idx = kmeans(im2double(all_feature��), N_bow); % ����Ŀ̫�࣬Kmean���������������
    % vl_feat�е�kmeans�ܺ���
    tic;
    fprintf('doing K-means with %d clusters...\n', numClusters);
    [cluster_centers, ~] = vl_kmeans(im2single(all_feature), numClusters, 'Initialization', 'plusplus');
    toc;
end

%% �����ںϷ���
for seq=1:size(feature_in_1_cam,1)
    for i_p=1:size(feature_in_1_cam{seq},1)
        % ĳ������p��ͼ������
        p_frames = feature_in_1_cam{seq}(i_p,:);
        % ������һЩ�վ�����Ҫȥ��
        p_frames = p_frames(~isemptycell(p_frames));
        if isempty(p_frames)
            continue
        end
        
        tmp = im2single(cell2mat(p_frames));
        [~, assign] = min(vl_alldist(tmp, cluster_centers), [], 2);
        combine_fe = zeros(numClusters, 1);
        for x=assign' % ������������Ӧλ��
            combine_fe(x) = combine_fe(x) + 1;
        end
        % ����Ҫ���й�һ��
        combine_fe = combine_fe/sum(combine_fe);

        % ������˵����н��������ں�
        % ----------------------------------------- %
        if 0
            % ����ӡ�ȡ��ֵ�ϳ�combine_fe
            combine_fe = zeros(size(p_frames(1)));
            for i_f=1:p_frames_num
                combine_fe = combine_fe + p_frames{i_f};
            end
            % ������ʾһ�ºϳɺ��ͼƬ
            combine_fe = combine_fe/p_frames_num; % imshow(combine_fe)
        end
        
        % ----------------------------------------- %
        % ������feature_in_1_cam_comb��
        feature_in_1_cam_comb{seq}{i_p,1} = combine_fe;
    end
end


end

