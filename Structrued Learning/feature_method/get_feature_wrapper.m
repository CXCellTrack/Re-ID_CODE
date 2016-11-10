function feature = get_feature_wrapper(dataaddr)

% dataaddr = 'E:\人脸识别数据集\ChokePoint\cropped_face\P1E_cropped_face\P1E\';
sub_dir = dir(dataaddr);

% feature 存储特征
feature = cell(4,3);
for seq=1:4 % 序列编号
    for cam=1:3 % 相机编号
        % 文件的dir会有.和.. 所以要去掉2个
        scpath = [dataaddr, sub_dir(seq*cam+2).name, '\'];
        person_dir = dir(scpath); 
        feature{seq,cam} = cell(numel(person_dir), 1);
        
        % 这个的dir有.和..和 .directory（可能有）
        for p=3:numel(person_dir)
            if strcmp(person_dir(p).name, '.directory')
                continue
            end
            fprintf('sequence %d, cam %d, person %s...\n', seq, cam, person_dir(p).name);
            scppath = [scpath, person_dir(p).name];
            pgm_dir = dir([scppath, '\*.pgm']); 
            % 对其中的每一张图片
            for i_p=1:numel(pgm_dir)
                im = imread([scppath, '\', pgm_dir(i_p).name]);
                % --------------------------------- % 
                % 对im进行特征抽取的操作
%                 im_feature = struct;
%                 % 1、提取全图像素点特征
%                 im_feature.raw = get_feature_raw(im2double(im));
                % 2、vl_sift描述子 
                [~, im_feature] = vl_sift(im2single(im));
                
                % --------------------------------- %
                % 把im_feature填充进cell中
                person_id = str2double(person_dir(p).name);
                feature{seq,cam}{person_id, i_p} = im_feature;
                
            end
        end
    end
end
         


                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
        
    