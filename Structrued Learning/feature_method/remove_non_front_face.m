function feature = remove_non_front_face(feature, is_front)

fnames = fieldnames(feature)';
for fn=fnames
    fe = getfield(feature, fn{1});
    is_f = getfield(is_front, fn{1});
    for i=1:numel(fe)
        for j=1:numel(fe{i})
            if is_f{i}{j}~=1
                fe{i}{j} = [];
            end
        end
    end
    feature = setfield(feature, fn{1}, fe);
end