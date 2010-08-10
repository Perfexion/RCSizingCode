function RCSizingCode()
    %RCSizingCode - This code is used for dimesnsioning of an RC sized
    %aircraft. This code is designed to work with and accept only the 'IPS'
    %unit standard.

    %% Initialization
    mtrwt = 2; %lbs
    btrywt = 2; %lbs
    %% Begin Main Loop



end

%% Sub Functions
function fk = DragBuildup(type, dim1, dim2, dim3, V)
    %DragBuildup - the purpose of the drag buildup function is to accept
    %   different solid object types (Body, Wing, or Nacelle) with proper dimensions
    %   and return the equivalent flat plate area (dimensionalized) for use in
    %   determining the net drag of a collection of solids which make up the
    %   external geometry of an aircraft. This function is a sub function of
    %   the RCSizingCode function and is designed to work within that function
    %   alone. It must be recognized that in order to keep the required input
    %   argument list minimal, the wing planform area is not accepted. As such
    %   the returned flat plate drag area (f) is dimensional. In order to get
    %   a true drag coefficient one must divide the flat plate area by the
    %   wing planform area. (CDf = f/S;)
    %
    %Four arguments are accepted into this function:
    %
    %argument one - type:
    %   This argument is the type of solid object which drag should be calculated for.
    %   DragBuildup currently supports three types of solid objects. 'Body'
    %   'Wing' or 'Nacelle'
    %
    %argument two - dim1: 
    %   This argument is the first dimension of the solid object. For a
    %   Body or Nacelle this is the length. For a wing this is the mean
    %   aerodynamic chord
    %
    %argument three dim2: 
    %   This argument is the second dimension of the solid object. For a
    %   Body or Naclle this is the diameter of the object. For a wing this
    %   is the thickness ratio of the wing.
    %
    %argument four dim3:
    %   This argument is the third dimension of the solid object. For a
    %   wing this should be the span. Because a body or nacelle is fully
    %   defined geometrically (as far as the drag buildup equations are
    %   concerned) this dimension is unused and should be provided as 0. 
    %
    %argument five V: 
    %   This argument is the velocity of the component (aircraft). This
    %   will likely be the estimated (target) cruise velocity of the
    %   aircraft.
    %
    %It should be known that all atmospheric conditions are based on the
    %standard atmosphere equations provided David G. Hull's "Fundamentals
    %of Airplane Flight Mechanics" for an atmosphere altitude of 2300'
    %(roughly the altitude of Tucson, AZ).
    
    if nargin ~= 5
        error('The number of input arguments is not consisten with the number of arguments required. Please see the function documentation for help on this issue');
    end
    
    %The following three equations come from David G. Hull's Book at page
    %45
    AtmTmp = 518.69 - 3.5662*10^-3 * 2300; %Rankine
    AtmPrs = 1.1376*10^-11 * AtmTmp ^ 5.2560; %lb/ft^2
    AtmDns = 6.6277*10^-15 * AtmTmp ^ 4.2560; %slugs/ft^3
    %The following value was found at:
    %http://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_6
    %01.html
    Mu = 3.90*10^-7; %(lbf-s)/ft^2
    
    %Again, from David G. Hull's Book
    Rek = AtmDns*V*dim1/Mu; %Unitless 
    
    %The folliwing equation is from John D. Anderson's "Modern Compressible
    %Flow" page 76
    Mach = sqrt(1.4*AtmPrs/AtmDns); %Unitless
    
    %All of the equations are from David G. Hull's Book "Fundamentals of
    %Airplane Flight Mechanics" Chapter 3.6.1
    switch type 
        case 'Body'
            Cfk = 0.455/(log10(Rek))^2.58;
            CFk = (1.0 + 0.2*Mach^2)^-0.467;
            FFk = 1.0 + 60/(dim1/dim2)^3 + 0.0025*(dim1/dim2);
            IFk = 1.2;
            Swetk = dim1*pi*dim2;
        
        case 'Nacelle'
            Cfk = 0.455/(log10(Rek))^2.58;
            CFk = (1.0 + 0.2*Mach^2)^-0.467;
            FFk = 1.0 + 60/(dim1/dim2)^3 + 0.0025*(dim1/dim2);
            IFk = 1.3;
            Swetk = dim1*pi*dim2;
        case 'Wing'
            Cfk = 0.455/(log10(Rek))^2.58;
            CFk = (1.0 + 0.2*Mach^2)^-0.467;
            FFk = 1.0 + 1.6/(dim2/dim1) + 100*(dim2/dim1)^4;
            IFk = 1.15;
            Swetk = dim1*dim3*2.04;
        otherwise
        error('The provided solid componet type does not match any solid component types known by this function. Please see the function documentation for help on this issue');
    end
    fk = Cfk*CFk*IFk*FFk*Swetk; 
end