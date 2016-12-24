
pwd
home = pwd;
subpath = genpath(home);
subpath = strsplit(';', subpath(1:end-1)); % ȥ��ĩβ�Ŀո�

newpath = [];
for i=1:numel(subpath)
    if ~isempty( strfind(subpath{i}, '.git') ) % ��Ҫ����git�ļ�
        subpath{i} = '';
        newpath = [ newpath, subpath{i} ];
    else
        newpath = [ newpath, ';', subpath{i} ];
    end
end

if newpath(1) == ';' % ȥ����ͷ���ܴ��ڵĿո�
    newpath(1) = '';
end

addpath(newpath);
savepath;
disp([home, ' ·����ӳɹ�']);
clear
