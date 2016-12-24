function feature_new = convert_feature(feature)
% 将feature转换一下结构
% item从P1E等变为seq1-4

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
        


