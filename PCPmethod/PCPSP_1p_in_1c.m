function [weight, seq_mustlink, seq_cannotlink] = PCPSP_1p_in_1c( W, n_person_in_seq, labels, max_must_link, max_cannot_link )

% �������� Pairwise Constraint Propagation by Semidefinite Programming
% feature Ϊ n*1 ����ʽ
% nΪ�����ĸ���
% feature{i,1}�д�ŵ��ǳ���i�������˸��Ե��������ںϽ����
% feature{i,1}{p}��������p����������Ϊ�գ�˵�������ڴ˴�������
% ��pҲ����������ı�ǩ��Ϣ

%% 1����W����W�Ǵ���2����֮������ƶȣ���ô�϶�ֵԽ�����ƶ�Խ�ߣ���WΪ����ķ��ȣ���ʽ24
n_seq = numel(n_person_in_seq); % ���еĸ���

% ÿ�������е�����
n_all = sum(n_person_in_seq);

%% 2������L����normalized graph Laplacian��
D = zeros(n_all);
for i=1:n_all
    D(i,i) = sum(W(i,:));
end
L = D - W;
L_ = D^(-0.5)*L*D^(-0.5); % L_ ��ʽ��8�� 
% L_�ǶԳư������ģ�����ֵ��0-2֮�䣨������eig(L_)���һ�£�

%% �趨must-link��cannot-linkԼ��
% �����ۼ�
n_person_cum = cumsum(n_person_in_seq);
n_person_cum = [0, n_person_cum]; % ��0��������һ��

% �涨ÿ2����must_link�ĸ�������������Ӿ�ȷ��
% max_must_link = 4;
% max_cannot_link = 4;
rng(0)

MustLink = [];
CannotLink = [];

for ii=1:n_seq-1
    for jj=ii+1:n_seq
        n_must_link = 0;
        n_cannot_link = 0;
        % -------- �������һЩindex��MustLink����֤label��ͬ��
        index = randsample(n_person_in_seq(ii), max_must_link)';
        for ii_ind=index 
            final_i = n_person_cum(ii) + ii_ind;
            % ��ii���ҵ���Щindex��Ӧ�ı�ǩ������jj���Ƿ�����ͬ��ǩ�ģ���������Ϊmustlink
            label_ii = labels{ii}(ii_ind);
            if any(labels{jj}==label_ii) && n_must_link<max_must_link
                jj_ind = find(labels{jj}==label_ii);
                % ����MustLink�е���������
                final_j = n_person_cum(jj) + jj_ind;
                MustLink = [MustLink; [final_i, final_j]];
                n_must_link = n_must_link + 1;
            end
        end
        
        % -------- �������һЩindex��CannotLink����֤label��ͬ��
        index = randsample(n_person_in_seq(ii), max_cannot_link)';
        for ii_ind=index 
            label_ii = labels{ii}(ii_ind);
            final_i = n_person_cum(ii) + ii_ind;
            not_equal_ind = find(labels{jj}~=label_ii);
            jj_ind = not_equal_ind( randsample(numel(not_equal_ind),1) );
            % ����MustLink�е���������
            final_j = n_person_cum(jj) + jj_ind;
            if n_cannot_link<max_cannot_link
                CannotLink = [CannotLink; [final_i, final_j]];
                n_cannot_link = n_cannot_link + 1;
            end
        end
    end
end

% ���mustlink��cannotlink�Ƿ�ì��
% һ������ì�ܣ��������Ϊ���ɽ�����
seq_mustlink = cell(n_seq);
for ii=1:size(MustLink,1)
    p1 = MustLink(ii,1);
    p2 = MustLink(ii,2);
    seq1 = find(p1-n_person_cum>0, 1, 'last');
    seq2 = find(p2-n_person_cum>0, 1, 'last');
    sp1 = p1-n_person_cum(seq1);
    sp2 = p2-n_person_cum(seq2);
    % tmp_mustlinkΪ����seq��mustlink���������ڼ��
    seq_mustlink{seq1,seq2} = [seq_mustlink{seq1,seq2}; sp1,sp2];
    if labels{seq1}(sp1)~=labels{seq2}(sp2)
        error('MustLink�ĵ�%d�д���%d��%d����mustlink��ϵ��', ii,p1,p2);
    end
end
seq_cannotlink = cell(n_seq);
for ii=1:size(CannotLink,1)
    p1 = CannotLink(ii,1);
    p2 = CannotLink(ii,2);
    seq1 = find(p1-n_person_cum>0, 1, 'last');
    seq2 = find(p2-n_person_cum>0, 1, 'last');
    sp1 = p1-n_person_cum(seq1);
    sp2 = p2-n_person_cum(seq2);
    % tmp_mustlinkΪ����seq��mustlink���������ڼ��
    seq_cannotlink{seq1,seq2} = [seq_cannotlink{seq1,seq2}; sp1,sp2];
    if labels{seq1}(sp1)==labels{seq2}(sp2)
        error('CannotLink�ĵ�%d�д���%d��%d����cannotlink��ϵ��', ii,p1,p2);
    end
end

%% ����������滮��Լ��
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
gap = 1e-6; % ʹ���ɳڵ�=1��=0����Ȼ����� numerical problem
exact_euqal = 1;
for i=1:n_all
    if exact_euqal
        E = zeros(n_all);
        E(i,i) = 1;
        F1 = [ F1, sum(sum(E.*K))==1 ]; % ����Ϊ1
    else
        F1 = [ F1, abs( K(i,i) -1 )<=gap]; % ����Ϊ1
    end
end
for ii=1:size(MustLink,1) % MustLinkΪ1
    if exact_euqal
        E = zeros(n_all);
        E(MustLink(ii,1), MustLink(ii,2)) = 1;
        F2 = [ F2, sum(sum(E.*K))==1 ];
    else
        F2 = [ F2, abs( K(MustLink(ii,1), MustLink(ii,2)) -1 )<=gap ];
    end
end
for ii=1:size(CannotLink,1) % CannotLinkΪ0
    if exact_euqal
        E = zeros(n_all);
        E(CannotLink(ii,1), CannotLink(ii,2)) = 1;
        F3 = [ F3, sum(sum(E.*K))==0 ];
    else
        F3 = [ F3, abs( K(CannotLink(ii,1), CannotLink(ii,2)) )<=gap ];
    end
end
% �������ͬһseq�и������cannotlink
for ii=1:n_seq
    for p1 = n_person_cum(ii)+1:n_person_cum(ii+1)-1
        for p2 = p1+1:n_person_cum(ii+1)
            F3 = [F3, K(p1,p2)==0];
        end
    end
end
         
toc;
F = [ K>=0, F1, F2, F3 ]; 

%% ���������滮

% ��ʼ���
disp('��ʼ���...');
options = sdpsettings('verbose', 0, 'debug',1, 'solver','mosek'); % ��ҵ����� mosek �� sedumi ��ܶ�
sol = solvesdp( F, OBJ, options ) % checkset(F)
Kval = value(K);


% ��kvalת�����������seq�ֿ������ڹ۲죩
weight = cell(n_seq);
for i1=1:n_seq-1
    for i2=i1+1:n_seq
        i1_s = n_person_cum(i1)+1;
        i1_e = n_person_cum(i1+1);
        i2_s = n_person_cum(i2)+1;
        i2_e = n_person_cum(i2+1);
        % ��Kval��ֵ��ΪȨֵ��������
        weight{i1,i2} = Kval(i1_s:i1_e, i2_s:i2_e);
    end
end










    
    
    
    