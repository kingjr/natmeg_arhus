function img = rotate_color(img, angle)
 img=rgb2hsv(img);
 img(:,:,1) = mod(img(:,:,1)+angle,1);
img=hsv2rgb(img);