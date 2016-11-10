function [phi_x_z, phi_x_zstar, loss] = CXSL_cal_phi_and_loss( vars_ADN_F, fe_to_use, gt_connect, truth_num )

% 解除vars
% vars = vars_ADN_F.vars; % 这样应该是无效的??
n_seq = numel(fe_to_use);

% 1、计算损失（先用0-1损失，只算recall）
% 使用trace计算相同点数目（同为1）
% TP_cell = cell(n_seq);
% TP_num = 0;
% for i=1:n_seq-1
%     for j=i+1:n_seq
%         TP_cell{i,j} = trace(gt_connect{i,j}'*vars_ADN_F.vars{i,j});
%         TP_num = TP_num + TP_cell{i,j};
%     end
% end
% loss = 1 - TP_num/truth_num;

% 换一种loss计算方式（2种无区别）
FN_num = 0;
for i=1:n_seq-1
    for j=i+1:n_seq
        FN_num = FN_num + sum(sum(gt_connect{i,j}.*(1-vars_ADN_F.vars{i,j})));
    end
end
loss = FN_num/truth_num;


% 2、计算phi_x_z = Σf*z，f为2张图的匹配特征
feature_diff = cell(n_seq);
phi_x_z = 0;
phi_x_zstar = 0; % 标准答案
for i=1:n_seq-1
    for j=i+1:n_seq
        n_p1 = size(fe_to_use{i},1);
        n_p2 = size(fe_to_use{j},1);
        feature_diff{i,j} = cell(n_p1, n_p2);
        for i_p1=1:n_p1
            for i_p2=1:n_p2
                % 对应位置2个特征之差
                feature_diff{i,j}{i_p1,i_p2} = abs(fe_to_use{i}(i_p1,:) - fe_to_use{j}(i_p2,:));
                phi_x_z = phi_x_z + vars_ADN_F.vars{i,j}(i_p1,i_p2) * feature_diff{i,j}{i_p1,i_p2};
                phi_x_zstar = phi_x_zstar + gt_connect{i,j}(i_p1,i_p2) * feature_diff{i,j}{i_p1,i_p2};
            end
        end
    end
end

phi_x_z = phi_x_z';
phi_x_zstar = phi_x_zstar';




















