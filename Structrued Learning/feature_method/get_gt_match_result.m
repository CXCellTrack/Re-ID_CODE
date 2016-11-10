function [ gt_connect, truth_num ] = get_gt_match_result( labels )

% 从labels中得到标准答案的连接方式
n_seq = numel(labels);
gt_connect = cell(n_seq);


% 如果每个观测当作一个clique，则需要对label进行转换
global each_person_a_clique;
if each_person_a_clique
    labels = num2cell(cell2mat(labels));
end

for ii = 1:n_seq-1
    for jj = ii+1:n_seq
        n_p1 = numel(labels{ii});
        n_p2 = numel(labels{jj});
        gt_connect{ii,jj} = zeros(n_p1, n_p2);
        for p1=1:n_p1
            p2 = find(labels{jj}==labels{ii}(p1));
            if ~isempty(p2)
                gt_connect{ii,jj}(p1,p2) = 1;
            end
        end        
    end
end

gt_cell = cellfun(@(x) sum(sum(x)), gt_connect, 'un',0);
truth_num = sum(sum((cell2mat(gt_cell))));

end

