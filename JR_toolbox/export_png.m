function img=export_png(name,xy)

drawnow;pause(.1);drawnow;
robo = java.awt.Robot;
t = java.awt.Toolkit.getDefaultToolkit();
rectangle = java.awt.Rectangle(t.getScreenSize());
image = robo.createScreenCapture(rectangle);
filehandle = java.io.File('tmp.png');
javax.imageio.ImageIO.write(image,'png',filehandle);
img=imread('tmp.png');delete('tmp.png');
if nargin>1
    img = img(max(1,xy(1)):min(end,xy(2)),max(1,xy(3)):min(end,xy(4)),:);
end
imwrite(img,name, 'png');
if nargout == 1
    img = imread(name);
end