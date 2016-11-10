function vars_ADN_F = cal_MatchLift_constraints( X_in_mat, m_hat, n_person_in_seq )

disp('����Լ������...');
vars = sdpvar(size(X_in_mat,1));

%% ����Լ������
n_seq = numel(n_person_in_seq);
n_person_cum = [0, cumsum(n_person_in_seq)];
% ��ѡԼ����X��ÿ��ÿ�� 0<=��<=1
% ע�����Լ��������Ϊ0-1����ʱ����ʡȥ���Լӿ�����ٶ�
% ��������ΪȨֵʱ����Ҫ���ϣ����������Լ�����ģ���
F0 = [];
for ii=1:n_seq-1
    for jj=ii+1:n_seq
        ii_s = n_person_cum(ii)+1;
        ii_e = n_person_cum(ii+1);
        jj_s = n_person_cum(jj)+1;
        jj_e = n_person_cum(jj+1);
        for ih=ii_s:ii_e % �ֿ��к�Լ��
            F0 = [ F0, sum(vars(ih, jj_s:jj_e))<=1, sum(vars(ih, jj_s:jj_e))>=0 ];
        end
        for iw=jj_s:jj_e % �ֿ��к�Լ��
            F0 = [ F0, sum(vars(ii_s:ii_e, iw))<=1, sum(vars(ii_s:ii_e, iw))>=0 ];
        end
    end
end

% Լ��1��X_out�ĶԽǷֿ�Ϊ��λ��
F1 = [];
for ii=1:size(vars,1)
    F1 = [F1, vars(ii,ii)==1];
end

% Լ��2��X>=0������X�Գƣ�ֻ��Ū1�뼴��
F2 = [];
for i=1:size(vars,1)
    for j=i:size(vars,1)
        F2 = [F2, vars(i,j)>=0];
    end
end

% Լ��3��������
X_out_ex = [m_hat, ones(1,size(vars,1));
            ones(size(vars,1),1), vars];
F3 = [ X_out_ex>=0 ];

%%
F = [ F0, F1, F2, F3 ];
vars_ADN_F.vars = vars;
vars_ADN_F.F = F;



