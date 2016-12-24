function feature_new = convert_feature(feature)
% ��featureת��һ�½ṹ
% item��P1E�ȱ�Ϊseq1-4

fields = fieldnames(feature);
feature_new = struct;
for i=1:4
    new_name = sprintf('seq%d', i);
    new_data = [];
    for j=1:numel(fields)
        data = getfield(feature, fields{j});
        new_data = [new_data; data(i,:)];
    end
    feature_new = setfield(feature_new, new_name, new_data);
end
        


