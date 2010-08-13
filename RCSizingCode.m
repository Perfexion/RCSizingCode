function RCSizingCode()
    %RCSizingCode - This code is used for dimesnsioning of an RC sized
    %aircraft. This code is designed to work with and accept only the 'IPS'
    %unit standard.

    %% Initialization
    Continue = true; %This is the flag which control the main while loop. Set it to false to break out of the loop;
    
    
    %% Input Parameters
    %In the future these parameters will be taken as inputs. This will be
    %implimented as soon as a GUI front end has been constructed. For the
    %sake of debugging and running the code without having to manually
    %input parameters, these will exists here. As such this section will
    %also serve as a collection of sorts of inputs which will need to be
    %accepted into the GUI. It should also be noted that the current chosen
    %values for these parameters, although not ridiculous, do not
    %necessarily reflect values which are optimum or nearly optimum for our
    %specific application. 
    
    % -- Weight Guess
    Weight = 10;
    
    % -- Flight Conditions;
    VTarget = 88; %Target Velocity (ft/s) - For scale, 88 corresponds to 60 MPH
    
    % -- Fuselage
    FuselageLength = 4; 
    FuselageDiameter = 8/12;
    
    % -- Main Wing
    WingAR = 9.5; %Aspect Ratio
    WingTPR = .3; %Taper Rati0 
    WingLoading = 1; %Wing loading (lbs/ft^2) - Simply wing area divided by weight. In the RC aircraft world a wing loading of 1 is generally regarded as an intermediate speed aircraft.
    WingAC = 1.5; %Location of aerodynamic center of wing behind nose (in ft)
    WingQCSweep = 15*pi/180; %Quarter Chord Sweep of the Wing
    WingTC = 0.1; %Wing Thickness Ratio
    
    % -- Horizontal Tail
    HTailAR = 4; %Aspect Ratio
    HTailTPR = .5; %Taper Ratio
    HTailVR = 1; %Volume Ratio
    HTailAC = 3.5; %Location of aerodynamic center of horizontal tail behind nose (in ft)
    HTailTC = 0.1; %Horizontal Tail Thickness Ratio
    
    % -- Vertical Tail
    VTailAR = 1.8; %Aspect Ratio
    VTailTPR = .8; %Taper Ratio
    VTailVR = 0.09; %Volume Ratio
    VTailAC = 3.5; %Location of aerodynamic center of vertical tail behind nose (in ft)
    VTailTC = 0.1; %Vertical Tail Thickness Ratio
    
    % -- Electronics
    MotorWeight = 2; % Weight of Motor (lbs)
    BatteryWeight = 4; % Weight of Battery (lbs)
    
    %% Atmosphere

    %The following three equations come from David G. Hull's Book at page
    %45
    AtmTmp = 518.69 - 3.5662*10^-3 * 2300; %Rankine
    AtmPrs = 1.1376*10^-11 * AtmTmp ^ 5.2560; %lb/ft^2
    AtmDns = 6.6277*10^-15 * AtmTmp ^ 4.2560; %slugs/ft^3

   

    
    %% Begin Main Loop
    while Continue == true
       %% Calculate Wing & Tail Dimensions
       
       % -- Main Wing
       WingS = Weight/WingLoading; %Wing planform area (ft^2)
       WingB =  sqrt(WingAR*WingS); %Wing Span (ft)
       WingCR = WingS/(WingB*(1+WingTPR)/2); %Wing Root Chord (ft)
       WingCT = WingTPR*WingCR; %Wing Tip Chord (ft)
       WingCBar = WingS/WingB; %Wing Mean Aerodynamic Chord (ft)
       
       % -- Horizontal Tail
       HTailS = HTailVR*WingB*WingS/(HTailAC-WingAC); %Horizontal tail planform area (ft^2) - assumes that the wing AC is directly over (or at least very near) to the aircraft CG. This may not be a good assumption
       HTailB = sqrt(HTailAR*HTailS); %Horizontal Tail Span (ft)
       HTailCR = HTailS/(HTailB*(1+HTailTPR)/2); %Horizontal Tail Root Chord (ft)
       HTailCT = HTailTPR*HTailCR; %Horizontal Tail Tip Chord (ft)
       HTailCBar = HTailS/HTailB; %Horizontal Tail Mean Aerodynamic Chord (ft)
       
       % -- Horizontal Tail
       VTailS = VTailVR*WingCBar*WingS/(VTailAC-WingAC); %Vertical tail planform area (ft^2) - assumes that the wing AC is directly over (or at least very near) to the aircraft CG. This may not be a good assumption
       VTailB = sqrt(VTailAR*VTailS); %Vertical Tail Span (ft)
       VTailCR = VTailS/(VTailB*(1+VTailTPR)/2); %Horizontal Tail Root Chord (ft)
       VTailCT = VTailTPR*VTailCR; %Horizontal Tail Tip Chord (ft)
       VTailCBar = VTailS/VTailB; %Horizontal Tail Mean Aerodynamic Chord (ft)
       
       %% Drag Buildup
       
       % -- Friction and Pressure Drag
       
       CDZero = 1.1*(DragBuildup('Body', FuselageLength, FuselageDiameter, 0, VTarget)+...
                          DragBuildup('Wing', WingCBar, WingTC, WingB, VTarget)+...
                          DragBuildup('Wing', HTailCBar, HTailTC, HTailB, VTarget)+...
                          DragBuildup('Wing', VTailCBar, VTailTC, VTailB, VTarget))/...
                          WingS; %Zero lift drag coefficient for wing, body, and tail
                      
       % -- Find Planform Efficiency Factor (k) 
       e = (1-0.024*WingAR^0.68)*(1-0.227*WingQCSweep^1.615); %I'm assuming this is a data fit equation. I do not have the history for this equation. I assume it came from somewhere in the MAE344 course. I got it from my sizing code from that class. I will do some research and see what I can find
       eprime = 1/(pi*WingAR*(1/(pi*WingAR*e)+0.38*CDZero)); %Same story here for this equation as for the equation for e. 
       K = 1/(pi*WingAR*eprime); %Planform Efficiency Factor - I think this version comes from David G. Hull's "Fundumentals of Airplane Flight Mechanics"
       
       % -- Find lift Coefficient 
       Cl = Weight*2/(AtmDns*VTarget^2*WingS); %Lift Coefficient
      
       % -- Final Drag Coefficient
       DragCoefficient = CDZero + Cl^2*K; %Drag Coefficient
       
       MinThrust = DragCoefficient*.5*AtmDns*VTarget^2*WingS; %Thrust Required to maintain VTarget; 
       
       %% Structural weight estimates
       
       %This section will estimate the weight of each component of the
       %aircraft based on the 
       
       %% Weight Summation
       
       %NewWeight = sum(); - TODO: Add weight of: Wing, Body, HTail, VTail,
       %Motor, Battery, Payload, Electronics
       
       
       %% Loop Break Test
       
       %{
       if abs(NewWeight-Weight)/Weight > .01
           Continue = false;
       end
       NewWeight = Weight;
       %}
       % TODO: Think about what needs to be done upon convergence;
       
       Continue = false; %Synthetic Loop Break - REMOVE THIS FOR WORKING VERSIONS!
                    
    end
end

%% Sub Functions

%% Drag Buildup Sub Function
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