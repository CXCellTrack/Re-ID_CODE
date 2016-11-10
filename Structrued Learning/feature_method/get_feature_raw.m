function im_feature = get_feature_raw(im)

% 直方图均衡化，去除亮度的影响
im = histeq(im);
im_feature = im(:);
