%% Edge shape functions of rectangular elements
% https://www.sciencedirect.com/topics/physics-and-astronomy/shape-functions
%% Test Rectangle 1
w = 1; % width (x)
h = 1; % height (x)
xc = 0.5;
yc = 0.5;

% values on the sides
a = 1;
b = 2;
c = 3;
d = -4;

%%
Ax = @(x,y) (1/h)*(yc + h/2 -y)*a + (1/h)*(y - yc + h/2)*b;
Ay = @(x,y) (1/w)*(xc + w/2 -x)*c + (1/w)*(x - xc + w/2)*d;

%%
px = [];
py = [];
vx = [];
vy = [];
for ix = 0:0.05:1
    for iy = 0:0.05:1
        px = [px; ix];
        py = [py; iy];
        vx = [vx; Ax(ix, iy)];
        vy = [vy; Ay(ix, iy)];
    end
end
%%
clf
hold on
for ii = 1:length(px)
    px1 = px(ii)+vx(ii)*0.01;
    py1 = py(ii)+vy(ii)*0.01;
    plot([px(ii) px1], [py(ii) py1],'r')
end
%% pick a triangle
tria = ginput(3);
x1 = tria(1,1); y1 = tria(1,2);
x2 = tria(2,1); y2 = tria(2,2);
x3 = tria(3,1); y3 = tria(3,2);
%%
plot(tria([1 2 3 1],1), tria([1 2 3 1],2))
%% Node shape functions of triangle
% http://www.uniroma2.it/didattica/TCM/deposito/2-D_elements.pdf#page=4
A = 0.5*det([1 x1 y1; 1 x2 y2; 1 x3 y3]);
N1 = @(x,y) (1/(2*A))*((x2*y3 - x3*y2) + (y2 - y3)*x + (x3 - x2)*y);
N2 = @(x,y) (1/(2*A))*((x3*y1 - x1*y3) + (y3 - y1)*x + (x1 - x3)*y);
N3 = @(x,y) (1/(2*A))*((x1*y2 - x2*y1) + (y1 - y2)*x + (x2 - x1)*y);
%% plot Node shape functions
a = 1; b = 2; c = 3;
xyz = [];
for ix = linspace(min(tria(:,1)),max(tria(:,1)),20)
    for iy = linspace(min(tria(:,2)),max(tria(:,2)),20)
        xyz = [xyz; ix iy N1(ix,iy)*a + N2(ix,iy)*b + N3(ix,iy)*c];
    end
end
%% 
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.') 
hold on
plot3(tria(:,1), tria(:,2), [a;b;c]);
%% Using Sym
syms x1 x2 x3 y1 y2 y3
A = 0.5*det([1 x1 y1; 1 x2 y2; 1 x3 y3]);
syms N1(x,y) N2(x,y) N3(x,y)
N1(x,y) = (1/(2*A))*((x2*y3 - x3*y2) + (y2 - y3)*x + (x3 - x2)*y);
N2(x,y) = (1/(2*A))*((x3*y1 - x1*y3) + (y3 - y1)*x + (x1 - x3)*y);
N3(x,y) = (1/(2*A))*((x1*y2 - x2*y1) + (y1 - y2)*x + (x2 - x1)*y);
N1x = diff(N1,x); N1y = diff(N1,y);
N2x = diff(N1,x); N2y = diff(N2,y);
N3x = diff(N1,x); N3y = diff(N3,y);
l12 = ((x2-x1)^2 + (y2-y1)^2)^0.5;
l23 = ((x3-x2)^2 + (y3-y2)^2)^0.5;
l31 = ((x1-x3)^2 + (y1-y3)^2)^0.5;
KS1 = (N1*(N2x + N2y) + N2*(N1x + N1y))/l12;
KS2 = (N2*(N3x + N3y) + N3*(N2x + N2y))/l23;
KS3 = (N3*(N1x + N1y) + N1*(N3x + N3y))/l31;
%%
tria = ginput(3);
x1 = tria(1,1); y1 = tria(1,2);
x2 = tria(2,1); y2 = tria(2,2);
x3 = tria(3,1); y3 = tria(3,2);
A = 0.5*det([1 x1 y1; 1 x2 y2; 1 x3 y3]);
syms N1(x,y) N2(x,y) N3(x,y)
N1(x,y) = (1/(2*A))*((x2*y3 - x3*y2) + (y2 - y3)*x + (x3 - x2)*y);
N2(x,y) = (1/(2*A))*((x3*y1 - x1*y3) + (y3 - y1)*x + (x1 - x3)*y);
N3(x,y) = (1/(2*A))*((x1*y2 - x2*y1) + (y1 - y2)*x + (x2 - x1)*y);
N1x = diff(N1,x); N1y = diff(N1,y);
N2x = diff(N1,x); N2y = diff(N2,y);
N3x = diff(N1,x); N3y = diff(N3,y);
l12 = ((x2-x1)^2 + (y2-y1)^2)^0.5;
l23 = ((x3-x2)^2 + (y3-y2)^2)^0.5;
l31 = ((x1-x3)^2 + (y1-y3)^2)^0.5;
KS1 = (N1*(N2x + N2y) + N2*(N1x + N1y))/l12;
KS2 = (N2*(N3x + N3y) + N3*(N2x + N2y))/l23;
KS3 = (N3*(N1x + N1y) + N1*(N3x + N3y))/l31;
%%
n12 = (tria(2,:) - tria(1,:))/norm(tria(2,:) - tria(1,:));
n23 = (tria(3,:) - tria(2,:))/norm(tria(3,:) - tria(2,:));
n31 = (tria(1,:) - tria(3,:))/norm(tria(1,:) - tria(3,:));
%% 
a = 1; b = 1; c = 1;
xyz = [];
for ix = linspace(min(tria(:,1)),max(tria(:,1)),20)
    for iy = linspace(min(tria(:,2)),max(tria(:,2)),20)
        xyz = [xyz; ix iy ...
                n12*eval(subs(subs(KS1,x,ix),y,iy))*a +  ...
                n23*eval(subs(subs(KS2,x,ix),y,iy))*b + ...
                n31*eval(subs(subs(KS3,x,ix),y,iy))*c];
    end
end
%%
clf
hold on
px1 = xyz(:,1) + xyz(:,3)*0.01;
py1 = xyz(:,2) + xyz(:,4)*0.01;
for ii = 1:length(px1)
    if inpolygon(xyz(ii,1), xyz(ii,2),tria(:,1), tria(:,2))
        plot([xyz(ii,1) px1(ii)], [xyz(ii,2) py1(ii)],'r')
    end
end


%% Edge shape functions of triangular elements
% http://www.iue.tuwien.ac.at/phd/nentchev/node43.html

%% RT shape functions for reference triangle
% http://users.jyu.fi/~oljumali/teaching/TIES594/08/mixed_finite_elements.pdf#page=5
% reference triangle
ref_tria = [0 0; 1 0; 0 1];
N1x = @(x)x; N1y = @(x)x-1; % Bottom edge
N2x = @(x) sqrt(2)*x; N2y = @(x)sqrt(2)*x; %diagonal edge
N3x = @(x)x-1; N3y = @(x)x; % left edge
%%
a = 0.2557211; b = -0.1561293; c = -0.04744832;
a = 0.3615364; b =  -0.2659625; c = -0.006270633;
xyv = [];
for ix = 0:0.05:1
    for iy = 0:0.05:1
        if inpolygon(ix, iy, ref_tria(:,1), ref_tria(:,2))
            xyv = [xyv; ix iy N1x(ix)*a + N2x(ix)*b + N3x(ix)*c ...
                    N1y(iy)*a + N2y(iy)*b + N3y(iy)*c];
        end
    end
end
% plot
clf
plot(ref_tria([1 2 3 1],1), ref_tria([1 2 3 1],2),'b')
hold on
for ii = 1:size(xyv,1)
    px = xyv(ii,1) + xyv(ii,3)*0.05;
    py = xyv(ii,2) + xyv(ii,4)*0.05;
    plot([xyv(ii,1) px],[xyv(ii,2) py],'r')
end
axis equal
%% RT shape functions for reference quad
ref_quad = [0 0; 0 1; 1 1; 1 0];
N1x = @(x)0; N1y = @(x)1-x; % Bottom face
N2x = @(x)x; N2y = @(x)0; % Right face
N3x = @(x)0; N3y = @(x)x; % Top face
N4x = @(x)1-x; N4y = @(x)0; % Left face
%%
a = 1; b = 5; c = 1; d = 1;
a = -19491.29434; b = 15202.0322; c = 0; d = -21030.3737;
a = 0.2659625; b = 0.1760972; c = -0.3749663; d = -0.1258288;
%a = -1+2*rand; b = -1+2*rand; c = -1+2*rand; d = -1+2*rand;
xyv = [];
for ix = 0:0.05:1
    for iy = 0:0.05:1
        xyv = [xyv; ix iy N1x(ix)*a + N2x(ix)*b + N3x(ix)*c + N4x(ix)*d ...
                    N1y(iy)*a + N2y(iy)*b + N3y(iy)*c + N4y(iy)*d];
    end
end
% plot
clf
plot(ref_quad([1 2 3 4 1],1), ref_quad([1 2 3 4 1],2),'b')
hold on
for ii = 1:size(xyv,1)
    px = xyv(ii,1) + xyv(ii,3)*0.15;
    py = xyv(ii,2) + xyv(ii,4)*0.15;
    plot([xyv(ii,1) px],[xyv(ii,2) py],'r')
end
axis equal







