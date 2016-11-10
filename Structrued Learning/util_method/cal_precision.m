function [ ACC, TP_cell ] = cal_precision( pre_connect, gt_connect, truth_num )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



% ȥ��vars�����½ǣ���Ϊ������ʱȥ����
% for seq2 = 1:n_seq-1
%     for seq1=seq2+1:n_seq
%         vars_value{seq1,seq2} = [];
%     end
% end

% ʹ��trace������ͬ����Ŀ��ͬΪ1��
TP_cell = cellfun(@(x,y) trace(x'*y), gt_connect, pre_connect, 'un',0);
pre_cell = cellfun(@(x) sum(sum(x)), pre_connect, 'un',0);

% ����recall��precision
pre_num = sum(sum((cell2mat(pre_cell))));
TP_num = sum(sum((cell2mat(TP_cell))));
recall = TP_num/truth_num;
precision = TP_num/pre_num;
f_m = 2/(1/precision+1/recall);

% �ϳ������
ACC.precision = precision;
ACC.recall = recall;
ACC.f_m = f_m;

disp(['truth num�� ',num2str(truth_num), '    ', 'pre num�� ',num2str(pre_num)])
disp(['ƥ��precisionΪ�� ', num2str(precision)]);
disp(['ƥ��recallΪ��    ', num2str(recall)]);
disp(['ƥ��f_mΪ��       ', num2str(f_m)]);




end

