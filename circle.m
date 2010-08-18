function [Z]=circle(center, NOP, radius)
if (nargin <2),
 error('Please see help for INPUT DATA.');
elseif (nargin==3)
    style='b-';
end;
THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
for i=1:NOP
    if i < NOP/8
        Xs(i) = 1;
        Ys(i) = (i-1)/(NOP/8);
    elseif i < 3*NOP/8
        Ys(i) = 1;
        Xs(i) = 1 - (i-(NOP/8))/(NOP/8);
    elseif i < 5*NOP/8
        Xs(i) = -1;
        Ys(i) = 1-(i-3*NOP/8)/(NOP/8);
    elseif i < 7*NOP/8
        Ys(i) = -1;
        Xs(i) = -1 + (i-5*NOP/8)/(NOP/8);
    else 
        Xs(i) = 1;
        Ys(i) = -1 + (i-7*NOP/8)/(NOP/8);
    end
end
X=X+center(1);
X=X';
Y=Y+center(2);
Y=Y';
Z=[X,Y, Xs', Ys'];
% plot(X,Y,Xs,Ys);
% [Xs ;Ys]
axis square;