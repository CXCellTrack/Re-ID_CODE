function feature = get_feature_wrapper(dataaddr)

% dataaddr = 'E:\����ʶ�����ݼ�\ChokePoint\cropped_face\P1E_cropped_face\P1E\';
sub_dir = dir(dataaddr);

% feature �洢����
feature = cell(4,3);
for seq=1:4 % ���б��
    for cam=1:3 % ������
        % �ļ���dir����.��.. ����Ҫȥ��2��
        scpath = [dataaddr, sub_dir(seq*cam+2).name, '\'];
        person_dir = dir(scpath); 
        feature{seq,cam} = cell(numel(person_dir), 1);
        
        % �����dir��.��..�� .directory�������У�
        for p=3:numel(person_dir)
            if strcmp(person_dir(p).name, '.directory')
                continue
            end
            fprintf('sequence %d, cam %d, person %s...\n', seq, cam, person_dir(p).name);
            scppath = [scpath, person_dir(p).name];
            pgm_dir = dir([scppath, '\*.pgm']); 
            % �����е�ÿһ��ͼƬ
            for i_p=1:numel(pgm_dir)
                im = imread([scppath, '\', pgm_dir(i_p).name]);
                % --------------------------------- % 
                % ��im����������ȡ�Ĳ���
%                 im_feature = struct;
%                 % 1����ȡȫͼ���ص�����
%                 im_feature.raw = get_feature_raw(im2double(im));
                % 2��vl_sift������ 
                [~, im_feature] = vl_sift(im2single(im));
                
                % --------------------------------- %
                % ��im_feature����cell��
                person_id = str2double(person_dir(p).name);
                feature{seq,cam}{person_id, i_p} = im_feature;
                
            end
        end
    end
end
         


                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
        
    