clear 
% ʹ����� Closing all open AVI files
clear mex

% ��ͼƬ������avi
avi_name = 'P1E_S1_C3';
aviobj = avifile(['C:\Users\Administrator\Desktop\',avi_name, '.avi'],...
    'FPS',30, 'COMPRESSION','Cinepak');
dataaddr = ['E:\����ʶ�����ݼ�\ChokePoint\', avi_name(1:end-3), '\', avi_name];
frames = dir([dataaddr, '\*.jpg']);

tic
for i = 1:numel(frames)
    if mod(i,100)==0
        disp(['����ͼ ', num2str(i), '...'])
    end
    pic = [dataaddr, frames(i).name];
    picdata = imread(pic);
    aviobj = addframe(aviobj, picdata);
end
aviobj = close(aviobj);
toc
