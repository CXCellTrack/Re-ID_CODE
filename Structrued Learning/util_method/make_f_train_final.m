function f_train_final = make_f_train_final( feature_train )

f_train_final = feature_train{1};
for i_f=2:numel(feature_train)
    for ii=1:numel(f_train_final)
        % ��2��cell������ͬ��Ⱥ���ܽ��д�ֱ�ϲ�
        ncol_1 = size(f_train_final{ii},2);
        % ����һ���߶Ȳ���50����ֹlabel����
        f_train_final{ii} = [f_train_final{ii}; cell(50*(i_f-1)-size(f_train_final{ii},1), ncol_1)];
        ncol_2 = size(feature_train{i_f}{ii},2);
        if ncol_2-ncol_1>=0
            f_train_final{ii} = [f_train_final{ii}, cell(size(f_train_final{ii},1), ncol_2-ncol_1)];
        else
            feature_train{i_f}{ii} = [feature_train{i_f}{ii}, cell(size(feature_train{i_f}{ii},1), ncol_1-ncol_2)];
        end
        % ���д�ֱ����
        f_train_final{ii} = [f_train_final{ii}; feature_train{i_f}{ii}];
    end
end