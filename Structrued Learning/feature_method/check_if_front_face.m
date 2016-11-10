function is_front = check_if_front_face( dataaddr )

% dataaddr = 'E:\人脸识别数据集\ChokePoint\cropped_face\P1E_cropped_face\P1E\';
sub_dir = dir(dataaddr);

% is_front 存储是否为正面
is_front = cell(4,3);
for seq=1:4 % 序列编号
    for cam=1:3 % 相机编号
        % 文件的dir会有.和.. 所以要去掉2个
        scpath = [dataaddr, sub_dir(seq*cam+2).name, '\'];
        person_dir = dir(scpath); 
        is_front{seq,cam} = cell(numel(person_dir), 1);
        
        % 这个的dir有.和..和 .directory（可能有）
        for p=3:numel(person_dir)
            if strcmp(person_dir(p).name, '.directory')
                continue
            end
            fprintf('sequence %d, cam %d, person %s...\n', seq, cam, person_dir(p).name);
            scppath = [scpath, person_dir(p).name];
            % 对其中的每个人的图片
            person_id = str2double(person_dir(p).name);
            % --------------------------------- % 
            % 调用python opencv判断该人脸是不是正面
            cmd = ['python "D:\Eclipse Mars\pycharm workspace\FaceDeteion\check_if_front_face.py" ', scppath];
            [status, result] = system(cmd);
            if status~=0
                error('调用python opencv出错！');
            end
            % --------------------------------- %
            for i_p=1:2:numel(result)
                is_front{seq,cam}{person_id, (i_p+1)/2} = str2num(result(i_p));
            end
        end
    end
end