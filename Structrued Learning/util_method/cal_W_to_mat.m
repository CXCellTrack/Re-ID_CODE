function [ W ] = cal_W_to_mat( fe_to_use, n_person_in_seq )


%��W�Ǵ���2����֮������ƶȣ���ô�϶�ֵԽ�����ƶ�Խ�ߣ���WΪ����ķ��ȣ���ʽ24

n_all = sum(n_person_in_seq);
% �ϳ����� feature_all ���󣬲��������
feature_all = cell2mat(fe_to_use);
distance = squareform( pdist(feature_all, 'euclidean') ); % ��||xi-xj||

% ��r��r�������нڵ㵽������ĵ�10���ڵ�֮���ƽ������
n_tenth = 10; % ������Ϊ��Ƚ��٣����ȡ4
tenth = zeros(n_all,1);
for i=1:n_all
    [Y,~]= sort(distance(i,:));
    tenth(i) = Y(n_tenth+1);
end
r = mean(tenth);

% S(0.1r, r, 5) �� S(r, 10r, 5)
S_r5 = linspace(0.1*r,r,5);
S_10r5 = linspace(r,10*r,5);
sigma_range = [S_r5(1:4), S_10r5]; % ���sigma�ķ�Χ

sigma = sigma_range(5); % ��ȡһ����Ϊsigma��ֵ
W = zeros(n_all);
for i=1:n_all-1
    for j=i+1:n_all
        W(i,j) = exp(-(distance(i,j)^2/(2*sigma^2))); % W�Ĺ�ʽ�����Լ���ƣ�����wԽ��2��Խ���Ƽ���
        W(j,i) = W(i,j); % ���öԳ�����W
    end
end

end

