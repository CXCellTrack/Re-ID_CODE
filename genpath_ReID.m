
pwd
home = pwd;
subpath = genpath(home);
subpath = strsplit(';', subpath(1:end-1)); % 去掉末尾的空格

newpath = [];
for i=1:numel(subpath)
    if ~isempty( strfind(subpath{i}, '.git') ) % 不要包含git文件
        subpath{i} = '';
        newpath = [ newpath, subpath{i} ];
    else
        newpath = [ newpath, ';', subpath{i} ];
    end
end

if newpath(1) == ';' % 去掉开头可能存在的空格
    newpath(1) = '';
end

addpath(newpath);
savepath;
disp([home, ' 路径添加成功']);
clear
