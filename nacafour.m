% --- The naca4 function creates x and y coordinates of the specified
% --- naca airfoil.

function [x,y] = nacafour(naca4,npanel,x,y)

%
% naca airfoil digits
%

n1 = floor(mod(naca4,10));
n2 = floor(mod(naca4,100)/10);
n3 = floor(mod(naca4,1000)/100);
n4 = floor(mod(naca4,10000)/1000);

%
% maximum camber, thickness and location of maximum camber
%

m = n4 / 100;
p = n3 / 10;
t = (n2*10 +n1) /100;

%
% compute thickness and camber distributions
%

if mod(npanel,2) ~= 0
    sprintf('Please choose an even number of panels');
    sprintf('Exiting...');
    exit;
end

nside = npanel / 2 +1;

%
% bunching parameter
%

an    = 1.5;
anp   = an +1;

%
% camber distribution
%

for i=1:nside
    frac  = (i -1)/(nside -1);
    xx(i) = 1 -anp*frac*(1 -frac)^an -(1 -frac)^anp;
    yt(i) = ( 0.29690*sqrt(xx(i)) -0.12600*xx(i)     ...
             -0.35160*xx(i)^2      +0.28430*xx(i)^3  ... 
             -0.10150*xx(i)^4) * t / 0.20;
    if xx(i) < p
        yc(i) = m/p^2 * (2*p*xx(i) -xx(i)^2);
    else
        yc(i) = m/(1 -p)^2 * ((1 -2*p) + 2*p*xx(i)-xx(i)^2);
    end
end

%
% airfoil shape = camber + thickness
%

for i=1:nside
    x(nside+i-1) = xx(i);
    x(nside-i+1) = xx(i);
    y(nside+i-1) = yc(i) +yt(i);
    y(nside-i+1) = yc(i) -yt(i);

end
x(2:nside*2) = x;
x(1) = 1;
x(nside*2+1) = 1;
y(2:nside*2) = y;
y(1) = 0;
y(nside*2+1) = 0;

return