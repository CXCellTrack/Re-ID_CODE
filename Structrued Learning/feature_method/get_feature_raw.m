function im_feature = get_feature_raw(im)

% ֱ��ͼ���⻯��ȥ�����ȵ�Ӱ��
im = histeq(im);
im_feature = im(:);
