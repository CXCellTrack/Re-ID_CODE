function [ feature_in_1_cam_comb, labels ] = filter_nan_and_get_label( feature_in_1_cam_comb )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% �� feature_in_1_cam_comb ��ȥ��nan������ȡ��label
n_seq = numel(feature_in_1_cam_comb); % ���еĸ���
labels = cell(n_seq,1);
for seq=1:n_seq
    n_person = numel(feature_in_1_cam_comb{seq,1});
    notnan_flag = [];
    seq_label = (1:n_person)';
    for p=1:n_person
        % �ж��Ƿ�ΪNan
        if any([0 1]==numel(feature_in_1_cam_comb{seq,1}{p,1}))
            notnan_flag = [notnan_flag, 0];
        else
            notnan_flag = [notnan_flag, 1];
        end
    end
    % ȥ��nan��label��person��ʵ���
    feature_in_1_cam_comb{seq,1} = feature_in_1_cam_comb{seq,1}(logical(notnan_flag));
    labels{seq,1} = seq_label(logical(notnan_flag));
    
    % ��feature_comb{seq,1}��Ϊ����
    feature_in_1_cam_comb{seq,1} = cell2mat( cellfun(@(x) reshape(x,1,[]),feature_in_1_cam_comb{seq,1},'un',0) );
end


end

