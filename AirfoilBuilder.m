function profile = AirfoilBuilder(foil, n)
    %%AirfoilBuilder - This function accepts a string of string of text which
    %specifies the series and number of the airfoil in format "key-number"
    %as well as the number of points requested as output for that airfoil.
    %For example calling AirfoilBuilder('NACA-001', 100) would return the
    %profile of the NACA 0012 airfoilmade up of 100 points. 

    if nargin ~= 2
        error('The number of input arguments provided is not consistent with the number required for this function.');
    end
    
    switch foil
        case 'NACA 0010'
            x = linspace(0,1,n);
            y = 1/0.2*.1*(.2969.*sqrt(x)-0.1260.*(x)-0.3516*(x.^2)+0.2843*(x.^3)-0.1015*(x.^4));
            x2 = linspace(1,0,n);
            y2 = -1/0.2*.1*(.2969.*sqrt(x2)-0.1260.*(x2)-0.3516*(x2.^2)+0.2843*(x2.^3)-0.1015*(x2.^4));
    end
    profile = [x x2; y y2];
    profile = profile';
end