image1 = im2double(imread("data/Rainier1.png"));
image2 = im2double(imread("data/Rainier2.png"));
hom = H;
homInv = inv(H);
imgResult = Stitch(image1,image2,hom,homInv);
function [stitchedImage] = Stitch(image1,image2,hom,homInv)
width2 = size(image2,2);
height2 = size(image2,1);

tl = [1,1];
tr = [1,width2];
bl = [height2,1];
br = [height2,width2];

tlp = zeros(2,1);
trp = zeros(2,1);
blp = zeros(2,1);
brp = zeros(2,1);

[tlp(1), tlp(2)] = Project(tl(1),tl(2), homInv);
[trp(1), trp(2)] = Project(tr(1),tr(2), homInv);
[blp(1), blp(2)] = Project(bl(1),bl(2), homInv);
[brp(1), brp(2)] = Project(br(1),br(2), homInv);

xs = [tlp(2), trp(2), blp(2), brp(2)];
ys = [tlp(1), trp(1), blp(1), brp(1)];
x_max = max(xs);
x_min = min(xs);
y_max = max(ys);
y_min = min(ys);

stitchedImage = zeros(uint32(abs(x_max) + abs(x_min)) ,uint32(abs(y_max) + abs(y_min)),  3);
minX = min([size(stitchedImage,2) , size(image1,2)]);
minY = min([size(stitchedImage,1) , size(image1,1)]);
for i = 1: size(image1,2)
    for j = 1 : size(image1,1)
        for d = 1:3
        stitchedImage(j,i,d) = image1(j,i,d);
        end
    end
end
for i = 1: size(stitchedImage,2)
    for j = 1 : size(stitchedImage,1)
        for d = 1:3
            [xCheck ,yCheck] = Project(i,j,hom);
            if((xCheck <= size(image2,2) && xCheck > 0) && (yCheck <= size(image2,1) && yCheck>0))
               stitchedImage(j,i,d) = image2(ceil(yCheck),ceil(xCheck),d);
            end
        end
    end
end
imshow(stitchedImage);
end

function [x2,y2] = Project(x1,y1,h)
d = [x1 ; y1 ; 1];
result = h*d;
result = result / result(3);

x2 = result(1);
y2 = result(2);
end