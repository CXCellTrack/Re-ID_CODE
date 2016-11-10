function vars_ADN_F = cal_MatchLift_constraints( X_in_mat, m_hat, n_person_in_seq )

disp('构造约束条件...');
vars = sdpvar(size(X_in_mat,1));

%% 计算约束条件
n_seq = numel(n_person_in_seq);
n_person_cum = [0, cumsum(n_person_in_seq)];
% 可选约束：X的每行每列 0<=和<=1
% 注意这个约束在输入为0-1矩阵时可以省去，以加快求解速度
% 但在输入为权值时，需要加上（这种是我自己试验的！）
F0 = [];
for ii=1:n_seq-1
    for jj=ii+1:n_seq
        ii_s = n_person_cum(ii)+1;
        ii_e = n_person_cum(ii+1);
        jj_s = n_person_cum(jj)+1;
        jj_e = n_person_cum(jj+1);
        for ih=ii_s:ii_e % 分块行和约束
            F0 = [ F0, sum(vars(ih, jj_s:jj_e))<=1, sum(vars(ih, jj_s:jj_e))>=0 ];
        end
        for iw=jj_s:jj_e % 分块列和约束
            F0 = [ F0, sum(vars(ii_s:ii_e, iw))<=1, sum(vars(ii_s:ii_e, iw))>=0 ];
        end
    end
end

% 约束1：X_out的对角分块为单位阵
F1 = [];
for ii=1:size(vars,1)
    F1 = [F1, vars(ii,ii)==1];
end

% 约束2：X>=0，由于X对称，只需弄1半即可
F2 = [];
for i=1:size(vars,1)
    for j=i:size(vars,1)
        F2 = [F2, vars(i,j)>=0];
    end
end

% 约束3：半正定
X_out_ex = [m_hat, ones(1,size(vars,1));
            ones(size(vars,1),1), vars];
F3 = [ X_out_ex>=0 ];

%%
F = [ F0, F1, F2, F3 ];
vars_ADN_F.vars = vars;
vars_ADN_F.F = F;



