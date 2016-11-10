function is_front = check_if_front_face( dataaddr )

% dataaddr = 'E:\����ʶ�����ݼ�\ChokePoint\cropped_face\P1E_cropped_face\P1E\';
sub_dir = dir(dataaddr);

% is_front �洢�Ƿ�Ϊ����
is_front = cell(4,3);
for seq=1:4 % ���б��
    for cam=1:3 % ������
        % �ļ���dir����.��.. ����Ҫȥ��2��
        scpath = [dataaddr, sub_dir(seq*cam+2).name, '\'];
        person_dir = dir(scpath); 
        is_front{seq,cam} = cell(numel(person_dir), 1);
        
        % �����dir��.��..�� .directory�������У�
        for p=3:numel(person_dir)
            if strcmp(person_dir(p).name, '.directory')
                continue
            end
            fprintf('sequence %d, cam %d, person %s...\n', seq, cam, person_dir(p).name);
            scppath = [scpath, person_dir(p).name];
            % �����е�ÿ���˵�ͼƬ
            person_id = str2double(person_dir(p).name);
            % --------------------------------- % 
            % ����python opencv�жϸ������ǲ�������
            cmd = ['python "D:\Eclipse Mars\pycharm workspace\FaceDeteion\check_if_front_face.py" ', scppath];
            [status, result] = system(cmd);
            if status~=0
                error('����python opencv����');
            end
            % --------------------------------- %
            for i_p=1:2:numel(result)
                is_front{seq,cam}{person_id, (i_p+1)/2} = str2num(result(i_p));
            end
        end
    end
end