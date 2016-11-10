clear
close all

%% ��������
n1 = 5;
n2 = 10;
n3 = 40;
r1 = 10;
r2 = 30;
[data_n1, data_n2, data_n3] = make_data(n1, n2, n3, r1, r2, 1);

%% 1����W����W�Ǵ���2����֮������ƶȣ���ô�϶�ֵԽ�����ƶ�Խ�ߣ���WΪ����ķ��ȣ���ʽ24
n_all = n1 + n2 + n3;
data_all = [data_n1; data_n2; data_n3];
distance = squareform( pdist(data_all, 'euclidean') ); % ��||xi-xj||

% ��r��r�������нڵ㵽������ĵ�10���ڵ�֮���ƽ������
n_tenth = 4; % ������Ϊ��Ƚ��٣����ȡ4
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
        W(i,j) = exp(-(distance(i,j)^2/(2*sigma^2)));
        W(j,i) = W(i,j); % ���öԳ�����W
    end
end

%% 2������L����normalized graph Laplacian��
D = zeros(n_all);
for i=1:n_all
    D(i,i) = sum(W(i,:));
end
L = D - W;
L_ = D^(-0.5)*L*D^(-0.5); % L_ ��ʽ��8�� 
% L_�ǶԳư������ģ�����ֵ��0-2֮�䣨������eig(L_)���һ�£�

%% ���������滮
label = [ones(n1,1); ones(n2,1)*2; ones(n3,1)*3]; % û���õ�label
% ��3�ֽ���ȡ��
n_sample = 5;
sample1 = randi(n1, n_sample,1);
sample2 = randi(n2, n_sample,1)+n1;
sample3 = randi(n3, n_sample,1)+n2+n1;
% ����must-link��cannot-link����
MustLink = [sample1, sample3]; % MustLink = []; % û��mustlink��cannotlinkʱ�����K���ж�Ϊ1
CannotLink = [sample1, sample2]; % CannotLink = [];
%
% ��������
for i=1:size(MustLink,1)
    p1 = data_all(MustLink(i,1),:);
    p2 = data_all(MustLink(i,2),:);
    line([p1(1),p2(1)], [p1(2),p2(2)], 'linestyle', '-','color','r');
end
for i=1:size(CannotLink,1)
    p1 = data_all(CannotLink(i,1),:);
    p2 = data_all(CannotLink(i,2),:);
    line([p1(1),p2(1)], [p1(2),p2(2)], 'linestyle', '--', 'color','k');
end
hold off

% ----------------------------------- %
% ���б������䡢���
K = sdpvar(n_all);
% kij =<phi(xi),phi(xj)>F. Then the matrix K =[kij ]n��n
%  is symmetric and positive semidefinite, denoted
% by K>=0��KΪ�Գư���������˿���Ϊһ���˺���
OBJ = sum(sum(L_.*K));

% ����3��Լ��
F1 = [];
F2 = [];
F3 = [];

disp('����Լ������...');
tic;
for i=1:n_all
%     E = zeros(n_all);
%     E(i,i) = 1;
%     F1 = [ F1, sum(sum(E.*K))==1 ]; % ����Ϊ1
    F1 = [ F1, abs( K(i,i)-1 )<=1e-6 ]; % ����Ϊ1
end
for ii=1:size(MustLink,1) % MustLinkΪ1
%     E = zeros(n_all);
%     E(MustLink(ii,1), MustLink(ii,2)) = 1;
%     F2 = [ F2, sum(sum(E.*K))==1 ];
    F2 = [ F2, abs( K(MustLink(ii,1), MustLink(ii,2))-1 )<=1e-6 ];
end
for ii=1:size(CannotLink,1) % CannotLinkΪ0
%     E = zeros(n_all);
%     E(CannotLink(ii,1), CannotLink(ii,2)) = 1;
%     F3 = [ F3, sum(sum(E.*K))==0 ];
    F3 = [ F3, abs( K(CannotLink(ii,1), CannotLink(ii,2)) )<=1e-6 ];
end

F = [ K>=0, F1, F2, F3 ]; 
toc;
% ��ʼ���
disp('��ʼ���...');
options = sdpsettings('verbose', 0, 'debug',1, 'solver','SEDUMI');
sol = solvesdp( F, OBJ, options ) % checkset(F)
Kval = value(K);













    
    
    
    