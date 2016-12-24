function [ ACC, TP_cell ] = cal_precision( pre_connect, gt_connect, truth_num )

% 去除vars的右下角（改为在生成时去除）
% for seq2 = 1:n_seq-1
%     for seq1=seq2+1:n_seq
%         vars_value{seq1,seq2} = [];
%     end
% end

% 使用trace计算相同点数目（同为1）
TP_cell = cellfun(@(x,y) trace(x'*y), gt_connect, pre_connect, 'un',0);
pre_cell = cellfun(@(x) sum(sum(x)), pre_connect, 'un',0);

% 计算recall和precision
pre_num = sum(sum((cell2mat(pre_cell))));
TP_num = sum(sum((cell2mat(TP_cell))));
recall = TP_num/truth_num;
precision = TP_num/pre_num;
f_m = 2/(1/precision+1/recall);

% 合成最后结果
ACC.precision = precision;
ACC.recall = recall;
ACC.f_m = f_m;

disp(['truth num： ',num2str(truth_num), '    ', 'pre num： ',num2str(pre_num)])
disp(['匹配precision为： ', num2str(precision)]);
disp(['匹配recall为：    ', num2str(recall)]);
disp(['匹配f_m为：       ', num2str(f_m)]);

%% 计算一致性约束违背情况
if_debug = false;
n_seq = size(pre_connect, 1);
all_zuhe = nchoosek(1:n_seq, 3);
n_violations = 0;
for col=1:size(all_zuhe,1)
    i = all_zuhe(col, 1); 
    j = all_zuhe(col, 2); 
    k = all_zuhe(col, 3);
    [x_ij, y_ij] = find(pre_connect{i,j});
    [x_jk, y_jk] = find(pre_connect{j,k});
    [x_ik, y_ik] = find(pre_connect{i,k});
    for ii=1:numel(x_ij)
        tmp_i = x_ij(ii);
        tmp_j = y_ij(ii);
        tmp_k = y_jk(x_jk==tmp_j);
        if isempty(tmp_k)
            continue
        end
        tmp_i_from_k = x_ik(y_ik==tmp_k);
        % 查找i->j,j->k,k->i是否一致
        if tmp_i~=tmp_i_from_k
            n_violations = n_violations + 1;
            if if_debug
                fprintf(['find a violation %d:\n{%d, %d}=%d,%d\t', ...
                        '{%d, %d}=%d,%d\t', ...
                        '{%d, %d}=%d,%d\n'], ...
                        n_violations, i,j,tmp_i,tmp_j, j,k,tmp_j,tmp_k, k,i,tmp_k,tmp_i_from_k);
            end
        end
    end
end
fprintf('find %d violations\n', n_violations);


end

