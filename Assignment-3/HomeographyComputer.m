
threshold = 0.0002;

img1 = imread("Rainier1.png");
img2 = imread("Rainier2.png");

width1 = size(img1,1);
height1 = size(img1,2);
width2 = size(img2,1);
height2 = size(img2,2);

tl1 = [1,1];
tr1 = [1,width1];
br1 = [height1,width1];
bl1 = [height1,1];

tl2 = [1,1];
tr2 = [1,width2];
br2 = [height2,width2];
bl2 = [height2,1];
a = zeros(2,1);
b = zeros(2,1);
c = zeros(2,1);
d = zeros(2,1);
[a(1),a(2)] = Project(tl2(1),tl2(2),inv(H));
[b(1),b(2)] = Project(tr2(1),tr2(2),inv(H));
[c(1),c(2)] = Project(br2(1),br2(2),inv(H));
[d(1),d(2)]= Project(bl2(1),bl2(2),inv(H));

stitchedImage = zeros(457,737);


function [x2,y2] = Project(x1,y1,h)
x2 = (h(1,1)*x1 + h(1,2)*y1 + h(1,3)) / (h(3,1)*x1 + h(3,2)*y1 + h(3,3));
y2 = (h(2,1)*x1 + h(2,2)*y1 + h(2,3)) / (h(3,1)*x1 + h(3,2)*y1 + h(3,3));
end

function count = ComputeInlierCount(h,Matches, numMatches, inlierThreshold)
count = 0;
for i = 1: numMatches
    
     TP = Matches(i,:);
     x1 = TP(1);
     y1 = TP(2);
     x2 = TP(3);
     y2 = TP(4);
     
    [x1p,y1p] = Project(x1,y1,h);
    %[x2p,y2p] = Project(x2,y2,h);
    d = pdist2([x1p,y1p],[x2,y2], 'euclidean');
   % d = sqrt(square(x2p-x1p) + square(y2p-y1p));
    
    if(d<inlierThreshold)
        count = count +1;
    end 
end
end

%inlierTreshold is distance
function [hom, homInv] = RANSAC(matches, numMatches, Iterations, inlierTreshold)
hom = [0,0,0 ; 1,0,1; 1,1,1];
firstMatches = matches(:, [1,2]);
secondMatches = matches(:, [3,4]);

for i = 1: Iterations
    k = randperm(numMatches);
    pts1 = firstMatches(k(1:4), :);
    pts2 = secondMatches(k(1:4), :);
    correspondingYValues = 4;
    hom = homography(pts1(1,1),pts1(1,2),pts1(2,1),pts1(2,2),pts1(3,1),pts1(3,1),pts1(4,1),pts1(4,1), pts2(1,1),pts2(1,2),pts2(2,1),pts2(2,2),pts2(3,1),pts2(3,2),pts2(4,1),pts2(4,2));
    c = ComputeInlierCount(hom, matches, numMatches, inlierTreshold);
end

%finds the best homography from inlier counts
for i = 1: Iterations
     FinalC = ComputeInlierCount(hom, matches, numMatches, inlierTreshold);
    if(c ~= FinalC)
        k = randperm(numMatches);
    pts1 = firstMatches(k(1:4), :);
    pts2 = secondMatches(k(1:4), :);
    correspondingYValues = 4;
    hom = homography(pts1(1,1),pts1(1,2),pts1(2,1),pts1(2,2),pts1(3,1),pts1(3,1),pts1(4,1),pts1(4,1), pts2(1,1),pts2(1,2),pts2(2,1),pts2(2,2),pts2(3,1),pts2(3,2),pts2(4,1),pts2(4,2));
    end
end

homInv = inv(hom);
end

function H=homography(x1,y1,x2,y2,x3,y3,x4,y4 , xp1,yp1,xp2,yp2,xp3,yp3,xp4,yp4)
%This function will find the homography betweeb 4 points using svd 
A=[
-x1  -y1  -1   0    0    0   x1*xp1   y1*xp1   xp1;
 0    0    0 -x1   -y1  -1   x1*yp1   y1*yp1   yp1;
-x2  -y2  -1   0    0    0   x2*xp2   y2*xp2   xp2;
 0    0    0 -x2   -y2  -1   x2*yp2   y2*yp2   yp2;
-x3  -y3  -1   0    0    0   x3*xp3   y3*xp3   xp3;
 0    0    0 -x3   -y3  -1   x3*yp3   y3*yp3   yp3;
-x4  -y4   -1  0    0    0   x4*xp4   y4*xp4   xp4;
 0    0    0  -x4  -y4  -1   x4*yp4   y4*yp4   yp4];
 
[U,S,V] = svd(A);
 
 
H=V(:,end);
H=reshape(H,3,3);
end

function [a, cnt] = HarrisCornerDetector(img, sigma, thres)
%%applying sobel edge detector in the horizontal direction
fx = [-1 0 1;-1 0 1;-1 0 1];
Ix = filter2(fx,img);
% applying sobel edge detector in the vertical direction
fy = [1 1 1;0 0 0;-1 -1 -1];
Iy = filter2(fy,img);
Ix2 = Ix.^2;
Iy2 = Iy.^2;
Ixy = Ix.*Iy;
clear Ix;
clear Iy;
%applying gaussian filter on the computed value
h = fspecial('gaussian',[5,5],sigma);
Ix2 = filter2(h,Ix2);
Iy2 = filter2(h,Iy2);
Ixy = filter2(h,Ixy);

height = size(img,1);
width = size(img,2);

result = zeros(height,width);
R = zeros(height,width);
Rmax = 0;

% Computing harris matrix
for i = 1:height
    for j = 1:width
        M = [Ix2(i,j) Ixy(i,j);Ixy(i,j) Iy2(i,j)];
        R(i,j) = det(M)-0.01*(trace(M))^2; %Corner response
        if R(i,j)> Rmax
        Rmax = R(i,j);
        end;
    end;
end;

%Thersholding
cnt = 0;
for i = 2:height-1
    for j = 2:width-1
        if R(i,j) > thres*100000 && R(i,j) > R(i-1,j-1) && R(i,j) > R(i-1,j) && R(i,j) > R(i-1,j+1) && R(i,j) > R(i,j-1) && R(i,j) > R(i,j+1) && R(i,j) > R(i+1,j-1) && R(i,j) > R(i+1,j) && R(i,j) > R(i+1,j+1)
        result(i,j) = 1;
        cnt = cnt+1;
        end;
    end;
end;
[posc, posr] = find(result == 1);
a = [posc, posr];
cnt ;

end

function [MNew] = newMatches(h,Matches, numMatches, inlierThreshold)
count = 0;
    for i = 1: size(Matches,1)
        if(i< size(Matches,1))
         TP = Matches(i,:);
         x1 = TP(1);
         y1 = TP(2);
         x2 = TP(3);
         y2 = TP(4);

        [x1p,y1p] = Project(x1,y1,h);
        [x2p,y2p] = Project(x2,y2,h);
        d = pdist2([x1p,y1p],[x2p,y2p], 'euclidean');
       % d = sqrt(square(x2p-x1p) + square(y2p-y1p));

        if(d<inlierThreshold)
            Matches(i, :) = [];
            i = i-1;
           count = count +1;
        end
        end
    end
    MNew = Matches;
end