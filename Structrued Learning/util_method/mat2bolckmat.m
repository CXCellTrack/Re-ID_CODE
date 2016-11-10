function X_blockmat = mat2bolckmat( X_mat, n_person_cum )

% ������ת��Ϊcell�����ֿ����
n_seq = numel(n_person_cum) - 1;
X_blockmat = cell(n_seq);
for i1=1:n_seq-1
    for i2=i1+1:n_seq
        i1_s = n_person_cum(i1)+1;
        i1_e = n_person_cum(i1+1);
        i2_s = n_person_cum(i2)+1;
        i2_e = n_person_cum(i2+1);
        % ��Kval��ֵ��ΪȨֵ��������
        X_blockmat{i1,i2} = X_mat(i1_s:i1_e, i2_s:i2_e);
    end
end

end

