function [ feature_in_1_cam_comb, labels ] = filter_nan_and_get_label( feature_in_1_cam_comb )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% 从 feature_in_1_cam_comb 中去掉nan，并抽取出label
n_seq = numel(feature_in_1_cam_comb); % 序列的个数
labels = cell(n_seq,1);
for seq=1:n_seq
    n_person = numel(feature_in_1_cam_comb{seq,1});
    notnan_flag = [];
    seq_label = (1:n_person)';
    for p=1:n_person
        % 判断是否为Nan
        if any([0 1]==numel(feature_in_1_cam_comb{seq,1}{p,1}))
            notnan_flag = [notnan_flag, 0];
        else
            notnan_flag = [notnan_flag, 1];
        end
    end
    % 去除nan，label是person真实编号
    feature_in_1_cam_comb{seq,1} = feature_in_1_cam_comb{seq,1}(logical(notnan_flag));
    labels{seq,1} = seq_label(logical(notnan_flag));
    
    % 将feature_comb{seq,1}变为矩阵
    feature_in_1_cam_comb{seq,1} = cell2mat( cellfun(@(x) reshape(x,1,[]),feature_in_1_cam_comb{seq,1},'un',0) );
end


end

