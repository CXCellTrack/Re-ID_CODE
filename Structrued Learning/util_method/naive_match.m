function pre_connect = naive_match( Kval )


n_seq = size(Kval,1);
pre_connect = cell(n_seq);

for ii=1:n_seq-1
    for jj=ii+1:n_seq
        % ���ö���ͼƥ���㷨
        fprintf('Bipartite graph matching between seq-%d and seq-%d...\n', ii,jj); 
        pre_connect{ii,jj} = Bipartite_graph_matching( Kval{ii,jj} );
    end
end



function matching_res = Bipartite_graph_matching( graph_mat )

[h,w] = size(graph_mat);
res_var = binvar(h, w, 'full');
% ָ��Լ�������к�Ϊ1��
F_b = [];
for i=1:h
    F_b = [ F_b, sum(res_var(i,:))<=1 ];
end
for j=1:w
    F_b = [ F_b, sum(res_var(:,j))<=1 ];
end
% Ŀ�꺯����ƥ�����
OBJ = sum(sum(res_var.*graph_mat));

% ��� Bipartite graph matching preoblem
options = sdpsettings('verbose',0,'solver','cplex');
sol = solvesdp( F_b, -OBJ, options );

matching_res = value(res_var);









