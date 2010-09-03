function varargout = AirDevil(varargin)
% AIRDEVIL M-file for AirDevil.fig
%      AIRDEVIL, by itself, creates a new AIRDEVIL or raises the existing
%      singleton*.
%
%      H = AIRDEVIL returns the handle to a new AIRDEVIL or the handle to
%      the existing singleton*.
%
%      AIRDEVIL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AIRDEVIL.M with the given input arguments.
%
%      AIRDEVIL('Property','Value',...) creates a new AIRDEVIL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AirDevil_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AirDevil_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AirDevil

% Last Modified by GUIDE v2.5 03-Sep-2010 16:03:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AirDevil_OpeningFcn, ...
                   'gui_OutputFcn',  @AirDevil_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

    % --- Executes on button press in TDViewButton.
function TDViewButton_Callback(hObject, eventdata, handles)
% hObject    handle to TDViewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    figure(1)
    v = get(handles.vTailAirfoilButtonPanel,'SelectedObject');
    h = get(handles.hTailAirfoilButtonPanel,'SelectedObject');
    w = get(handles.wingAirfoilButtonPanel,'SelectedObject');
    if (get(handles.solidWorksCheckBox,'Value') == get(handles.solidWorksCheckBox,'Max'))
       % Checkbox is checked
       solidFlag=1;
    else
       % Checkbox is not checked
       solidFlag=0;
    end
    
    switch get(h,'Tag')   % Get Tag of selected object
    case 'NACA4RadioButton'
        N4=str2double(get(handles.hTailAirfoilTxt, 'String'));
        [hTailx,hTaily] = nacafour(N4,100);
        hTailAirfoil = [hTailx',hTaily'];
        HTailAirfoil = get(handles.hTailAirfoilTxt, 'String');
 
    case 'NACA5RadioButton'
        iaf.designation=get(handles.hTailNACA5AirfoilTxt, 'String');
        iaf.n=100;
        iaf.HalfCosineSpacing=1;
        iaf.wantFile=0;
        iaf.datFilePath='./'; % Current folder
        iaf.is_finiteTE=0;
        af = naca5gen(iaf);
        hTailAirfoil = [af.x,af.z];
        HTailAirfoil = get(handles.hTailNACA5AirfoilTxt, 'String');
 
    case 'customAirfoilRadioButton'
        global hTailAirfoil
        HTailAirfoil = 'custom';

    otherwise
       % Code for when there is no match.
       'otherwise'
 
    end
    
    switch get(v,'Tag')   % Get Tag of selected object
    case 'NACA4RadioButton'
        N4=str2double(get(handles.vTailAirfoilTxt, 'String'));
        [vTailx,vTaily] = nacafour(N4,100);
        vTailAirfoil = [vTailx',vTaily'];
        VTailAirfoil = get(handles.vTailAirfoilTxt, 'String');
 
    case 'NACA5RadioButton'
        iaf.designation=get(handles.hTailNACA5AirfoilTxt, 'String');
        iaf.n=100;
        iaf.HalfCosineSpacing=1;
        iaf.wantFile=0;
        iaf.datFilePath='./'; % Current folder
        iaf.is_finiteTE=0;
        af = naca5gen(iaf);
        vTailAirfoil = [af.x,af.z];
        VTailAirfoil = get(handles.vTailNACA5AirfoilTxt, 'String');
 
    case 'customAirfoilRadioButton'
        global vTailAirfoil
        VTailAirfoil = 'custom';

    otherwise
       % Code for when there is no match.
       'otherwise'
 
    end
    
    switch get(w,'Tag')   % Get Tag of selected object
    case 'NACA4RadioButton'
        N4=str2double(get(handles.wingAirfoilTxt, 'String'));
        [wingx,wingy] = nacafour(N4,100);
        wingAirfoil = [wingx',wingy'];
        WingAirfoil = get(handles.wingAirfoilTxt, 'String');
 
    case 'NACA5RadioButton'
        iaf.designation=get(handles.wingNACA5AirfoilTxt, 'String');
        iaf.n=100;
        iaf.HalfCosineSpacing=1;
        iaf.wantFile=0;
        iaf.datFilePath='./'; % Current folder
        iaf.is_finiteTE=0;
        af = naca5gen(iaf);
        wingAirfoil = [af.x,af.z];
        WingAirfoil = get(handles.wingNACA5AirfoilTxt, 'String');
    
 
    case 'customAirfoilRadioButton'
        global wingAirfoil
        WingAirfoil = 'custom';

    otherwise
       % Code for when there is no match.
       'otherwise'
 
    end
    
    %%Sizing Code 
    Continue = true;
    % -- Weight Guess
    Weight = str2double(get(handles.weightGuessTxt,'String')); %lbs
    ZCg = 0; %ZCg
    
    % -- Flight Conditions;
    VTarget =str2double(get(handles.vTargetTxt,'String')); %Target Velocity (ft/s) - For scale, 88 corresponds to 60 MPH
    
    % -- Fuselage
    FuselageLength = str2double(get(handles.fuselageLengthTxt,'String')); %ft
    FuselageDiameter = str2double(get(handles.fuselageDiameterTxt,'String')); %ft
    
    % -- Main Wing
    WingAR = str2double(get(handles.wingARTxt,'String')); %Aspect Ratio
    WingTPR = str2double(get(handles.wingTPRTxt,'String')); %Taper Rati0 
    WingLoading = str2double(get(handles.wingLoadingTxt,'String')); %Wing loading (lbs/ft^2) - Simply wing area divided by weight. In the RC aircraft world a wing loading of 1 is generally regarded as an intermediate speed aircraft.
    WingAC = str2double(get(handles.wingACTxt,'String')); %Location of aerodynamic center of wing behind nose (in ft)
    WingQCSweep = str2double(get(handles.wingQCSweepTxt,'String')); %Quarter Chord Sweep of the Wing
    WingTC = str2double(get(handles.wingTCTxt,'String')); %Wing Thickness Ratio
    WingDihedral = str2double(get(handles.wingDihedralTxt,'String')); %Wing Dihedral
    WingZ = str2double(get(handles.wingZTxt,'String')); %Wing Tail Z
    WingIncidence = str2double(get(handles.wingIncidenceTxt,'String')); %Wing Incidence Angle
    
    
    % -- Horizontal Tail
    HTailAR = str2double(get(handles.hTailARTxt,'String')); %Aspect Ratio
    HTailTPR = str2double(get(handles.hTailTPRTxt,'String')); %Taper Ratio
    HTailVR = str2double(get(handles.hTailVRTxt,'String')); %Volume Ratio
    HTailAC = str2double(get(handles.hTailACTxt,'String')); %Location of aerodynamic center of horizontal tail behind nose (in ft)
    HTailTC = str2double(get(handles.hTailTCTxt,'String')); %Horizontal Tail Thickness Ratio
    HTailQCSweep = str2double(get(handles.hTailQCSweepTxt,'String')); %Horizontal Tail Quarter Chord Sweep
    HTailDihedral = str2double(get(handles.hTailDihedralTxt,'String')); %Horizontal Tail Dihedral
    HTailZ = str2double(get(handles.hTailZTxt,'String')); %Horizontal Tail Z
    HTailIncidence = str2double(get(handles.hTailIncidenceTxt,'String')); %Horizontal Tail Incidence Angle
    
    
    % -- Vertical Tail
    VTailAR = str2double(get(handles.vTailARTxt,'String')); %Aspect Ratio
    VTailTPR = str2double(get(handles.vTailTPRTxt,'String')); %Taper Ratio
    VTailVR = str2double(get(handles.vTailVRTxt,'String')); %Volume Ratio
    VTailAC = str2double(get(handles.vTailACTxt,'String')); %Location of aerodynamic center of vertical tail behind nose (in ft)
    VTailTC = str2double(get(handles.vTailTCTxt,'String')); %Vertical Tail Thickness Ratio
    VTailQCSweep = str2double(get(handles.vTailQCSweepTxt,'String')); %Vertical Tail Quarter Chord Sweep
    
    
    % -- Electronics
    MotorWeight = str2double(get(handles.motorWeightTxt,'String')); % Weight of Motor (lbs) - Hacker B50-12XL weighs .74 lbs
    MotorCg = str2double(get(handles.motorCgTxt,'String')); % Weight of Motor
    BatteryWeight = str2double(get(handles.batteryWeightTxt,'String')); % Weight of Battery (lbs)
    BatteryCg = str2double(get(handles.batteryCgTxt,'String')); % Weight of Motor 
    PayloadWeight = str2double(get(handles.payloadWeightTxt,'String')); %Weight of Payload (lbs)
    PayloadCg = WingAC; 

 %% Atmosphere

    %The following three equations come from David G. Hull's Book at page
    %45
    AtmAlt = str2double(get(handles.atmAltTxt, 'String')); %Atmospherica Altitude
    AtmTmp = 518.69 - 3.5662*10^-3 * AtmAlt; %Rankine
    AtmPrs = 1.1376*10^-11 * AtmTmp ^ 5.2560; %lbf/ft^2
    AtmDns = 6.6277*10^-15 * AtmTmp ^ 4.2560; %slug/ft^3

    
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
       HTailS = HTailVR*WingCBar*WingS/(HTailAC-WingAC); %Horizontal tail planform area (ft^2) - assumes that the wing AC is directly over (or at least very near) to the aircraft CG. This may not be a good assumption
       HTailB = sqrt(HTailAR*HTailS); %Horizontal Tail Span (ft)
       HTailCR = HTailS/(HTailB*(1+HTailTPR)/2); %Horizontal Tail Root Chord (ft)
       HTailCT = HTailTPR*HTailCR; %Horizontal Tail Tip Chord (ft)
       HTailCBar = HTailS/HTailB; %Horizontal Tail Mean Aerodynamic Chord (ft)
       
       % -- Horizontal Tail
       VTailS = VTailVR*WingB*WingS/(VTailAC-WingAC); %Vertical tail planform area (ft^2) - assumes that the wing AC is directly over (or at least very near) to the aircraft CG. This may not be a good assumption
       VTailB = sqrt(VTailAR*VTailS); %Vertical Tail Span (ft)
       VTailCR = VTailS/(VTailB*(1+VTailTPR)/2); %Vertical Tail Root Chord (ft)
       VTailCT = VTailTPR*VTailCR; %Vertical Tail Tip Chord (ft)
       VTailCBar = VTailS/VTailB; %Vertical Tail Mean Aerodynamic Chord (ft)

       
       %% Drag Buildup
       
       % -- Friction and Pressure Drag
       FuselageDrag = DragBuildup('Body', FuselageLength, FuselageDiameter, 0, VTarget)/WingS; %Drag Coefficient due to Fuselage
       WingDrag = DragBuildup('Wing', WingCBar, WingTC, WingB, VTarget)/WingS; %Drag Coefficient due to Wing
       HTailDrag = DragBuildup('Wing', HTailCBar, HTailTC, HTailB, VTarget)/WingS; %Drag Coefficient due to HTail
       VTailDrag = DragBuildup('Wing', VTailCBar, VTailTC, VTailB, VTarget)/WingS; %Drag Coefficient due to VTail
       
       CdZero = 1.1*(FuselageDrag + WingDrag + HTailDrag + VTailDrag); %Zero lift drag coefficient for wing, body, and tail
                      
       % -- Find Planform Efficiency Factor (k) 
       e = (1-0.024*WingAR^0.68)*(1-0.227*WingQCSweep^1.615); %I'm assuming this is a data fit equation. I do not have the history for this equation. I assume it came from somewhere in the MAE344 course. I got it from my sizing code from that class. I will do some research and see what I can find
       eprime = 1/(pi*WingAR*(1/(pi*WingAR*e)+0.38*CdZero)); %Same story here for this equation as for the equation for e. 
       K = 1/(pi*WingAR*eprime); %Planform Efficiency Factor - I think this version comes from David G. Hull's "Fundumentals of Airplane Flight Mechanics"
       
       % -- Find lift Coefficient 
       Cl = Weight*2/(AtmDns*VTarget^2*WingS); %Lift Coefficient
      
       % -- Final Drag Coefficient
       Cd = CdZero + Cl^2*K; %Drag Coefficient
       
       MinThrust = Cd*.5*AtmDns*VTarget^2*WingS; %Thrust Required to maintain VTarget; 
       
       %% Structural weight estimates
       
       %This section will estimate the weight of each component of the
       %aircraft based on the component size. 
       BalsaDensity = .00576; %lb/in^3
       CarbonFiberDensity = .064667; %lb/in^3
       WingSparCSA = .0706; %in^2 based on .5in OD and .4 in ID
       FuselageSparCSA = .149226; %in^2 based in 1in OD and .9in ID
       WingBCSAE = 1/3*1/24; %Dimensionless
       FuselageBCSAE = ((FuselageDiameter*.8)^2/4*pi-(FuselageDiameter*.65)^2/4*pi)/...
           (FuselageDiameter*.8)^2/4*pi*1/24; %Dimensionless
       
       WingWeight = WingB*CarbonFiberDensity*WingSparCSA*12 + ...
           WingBCSAE*BalsaDensity*WingTC*WingS*12^3;
       HTailWeight = HTailB*CarbonFiberDensity*WingSparCSA*12 + ...
           WingBCSAE*BalsaDensity*HTailTC*HTailS*12^3;
       VTailWeight = VTailB*CarbonFiberDensity*WingSparCSA*12 + ...
           WingBCSAE*BalsaDensity*VTailTC*VTailS*12^3;
       FuselageWeight = FuselageLength*CarbonFiberDensity*FuselageSparCSA*12 + ...
           FuselageBCSAE*BalsaDensity*FuselageLength*FuselageDiameter^2/4*pi*12^3;
       
       %% Weight Summation
       
       NewWeight = WingWeight + HTailWeight + VTailWeight + FuselageWeight...
           + BatteryWeight + MotorWeight + PayloadWeight; %Sum all weights (lbm)
       
       
       %% Loop Break Test
       
       if abs(NewWeight-Weight)/Weight < .01
           Continue = false;
           %% Weight Balance Calcs
            XCg = (WingWeight*WingAC+HTailWeight*HTailAC+VTailWeight*VTailAC...
                +FuselageWeight*1/3*FuselageLength+MotorWeight*MotorCg+...
                BatteryWeight*BatteryCg)/...
                (WingWeight+HTailWeight+VTailWeight+FuselageWeight+...
                MotorWeight+BatteryWeight);
               %% Organize Function Output
               Aircraft.FuselageLength = double(FuselageLength);
               Aircraft.FuselageDiameter = FuselageDiameter;
               Aircraft.FuselageWeight = FuselageWeight;
               Aircraft.FuselageDrag = FuselageDrag;
               Aircraft.WingAC = WingAC;
               Aircraft.WingS = WingS;
               Aircraft.WingB = WingB;
               Aircraft.WingCR = WingCR;
               Aircraft.WingCT = WingCT;
               Aircraft.WingCBar = WingCBar;
               Aircraft.WingAR = WingAR;
               Aircraft.WingTPR = WingTPR;
               Aircraft.WingQCSweep = WingQCSweep;
               Aircraft.WingDihedral = WingDihedral;
               Aircraft.WingTC = WingTC;
               Aircraft.WingAirfoil = WingAirfoil;
               Aircraft.WingWeight = WingWeight;
               Aircraft.WingDrag = WingDrag;
               Aircraft.WingLoading = WingLoading;
               Aircraft.WingK = K;
               Aircraft.HTailAC = HTailAC;
               Aircraft.HTailS = HTailS;
               Aircraft.HTailB = HTailB;
               Aircraft.HTailCR = HTailCR;
               Aircraft.HTailCT = HTailCT;
               Aircraft.HTailCBar = HTailCBar;
               Aircraft.HTailAR = HTailAR;
               Aircraft.HTailTPR = HTailTPR;
               Aircraft.HTailQCSweep = HTailQCSweep;
               Aircraft.HTailDihedral = HTailDihedral;
               Aircraft.HTailTC = HTailTC;
               Aircraft.HTailAirfoil = HTailAirfoil;
               Aircraft.HTailWeight = HTailWeight;
               Aircraft.HTailDrag = HTailDrag;
               Aircraft.HTailVR = HTailVR;
               Aircraft.VTailAC = VTailAC;
               Aircraft.VTailS = VTailS;
               Aircraft.VTailB = VTailB;
               Aircraft.VTailCR = VTailCR;
               Aircraft.VTailCT = VTailCT;
               Aircraft.VTailCBar = VTailCBar;
               Aircraft.VTailAR = VTailAR;
               Aircraft.VTailTPR = VTailTPR;
               Aircraft.VTailQCSweep = VTailQCSweep;
               Aircraft.VTailTC = VTailTC;
               Aircraft.VTailAirfoil = VTailAirfoil;
               Aircraft.VTailWeight = VTailWeight;
               Aircraft.VTailDrag = VTailDrag;
               Aircraft.VTailVR = VTailVR;
               Aircraft.BatteryWeight = BatteryWeight;
               Aircraft.BatteryCg = BatteryCg;
               Aircraft.MotorWeight = MotorWeight;
               Aircraft.MotorCg = MotorCg;
               Aircraft.PayloadWeight = PayloadWeight;
               Aircraft.PayloadCg = PayloadCg;
               Aircraft.CdZerp = CdZero;
               Aircraft.MinThrust = MinThrust;
               Aircraft.AtmAlt = AtmAlt;
               Aircraft.AtmTmp = AtmTmp;
               Aircraft.AtmDns = AtmDns; 
               Aircraft.FlightWeight = NewWeight;
               Aircraft.WingIncidence = WingIncidence;
               Aircraft.WingZ = WingZ;
               Aircraft.HTailIncidence = HTailIncidence;
               Aircraft.HTailZ = HTailZ;
               Aircraft.XCg = XCg;
               Aircraft.ZCg = ZCg; 

       end
       Weight = NewWeight;      
    end
    assignin('base', 'Aircraft', Aircraft);
    WriteDatcom
    %%Optimization
    if (get(handles.optimizeCheckBox,'Value') == get(handles.optimizeCheckBox,'Max'))
        %this code will hopefully be used to optimize our plane, right now
        %it gets the values of ideal lift coefficient from the datcom output.

        %read datcom output
        text = fileread('test.out');
        alpha = regexp(text, '\W(?<=IDEAL ANGLE OF ATTACK =    .)\d\d\d\d\d', 'match');

        alpha_ideal_wing=str2double(alpha(1));
        
        set(handles.wingIncidenceTxt, 'String', alpha_ideal_wing*pi/180);
        Aircraft.WingIncidence = alpha_ideal_wing*pi/180;
        assignin('base', 'Aircraft', Aircraft);
        WriteDatcom
        
        text = fileread('test.out');
        epsilon = regexp(text, '\W(?<= .0      .....      ).....', 'match');
        alpha = regexp(text, '\W(?<=IDEAL ANGLE OF ATTACK =    .)\d\d\d\d\d', 'match');
        
        epsilon = str2double(epsilon);
        alpha_ideal_htail=str2double(alpha(2));
        
        set(handles.hTailIncidenceTxt, 'String', (epsilon + alpha_ideal_htail)*pi/180);
        Aircraft.HTailIncidence = (epsilon+alpha_ideal_htail)*pi/180;
        
        assignin('base', 'Aircraft', Aircraft);
        WriteDatcom
        
        stopcg = false;
        count = 0;
        while stopcg == false
        text = fileread('test.out');
        Cm = regexp(text, '\W(?<= .0    .....    .....    )......', 'match');
        Cm = str2double(Cm)
        
            if Cm > .003
                XCg = XCg * .95;
            elseif Cm < -.003
                XCg = XCg * 1.05;
            else
               stopcg = true;
            end
            Aircraft.XCg = XCg;
            count = count + 1
            assignin('base', 'Aircraft', Aircraft);
            WriteDatcom
        end
    else
       %who knows
    end
    
    BatteryCg = (XCg*(WingWeight+HTailWeight+VTailWeight+FuselageWeight+...
    MotorWeight+BatteryWeight)-(WingWeight*WingAC+HTailWeight*HTailAC...
    +VTailWeight*VTailAC+FuselageWeight*1/3*FuselageLength+...
    MotorWeight*MotorCg))/BatteryWeight;
    Aircraft.BatteryCg = BatteryCg;
    set(handles.batteryCgTxt, 'String', BatteryCg);
    assignin('base', 'Aircraft', Aircraft);

        
        
        
        
        
        %% Run TDVIew
        TDView(solidFlag, wingAirfoil, vTailAirfoil, hTailAirfoil)
    
    
    

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
    AtmPrs = 1.1376*10^-11 * AtmTmp ^ 5.2560; %lbf/ft^2
    AtmDns = 6.6277*10^-15 * AtmTmp ^ 4.2560 * (32.1740486); %lbm/ft^3
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

    % --- Executes just before AirDevil is made visible.
function AirDevil_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AirDevil (see VARARGIN)

% Choose default command line output for AirDevil
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes AirDevil wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AirDevil_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fuselageLengthTxt_Callback(hObject, eventdata, handles)
% hObject    handle to fuselageLengthTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fuselageLengthTxt as text
%        str2double(get(hObject,'String')) returns contents of fuselageLengthTxt as a double


% --- Executes during object creation, after setting all properties.
function fuselageLengthTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fuselageLengthTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fuselageDiameterTxt_Callback(hObject, eventdata, handles)
% hObject    handle to fuselageDiameterTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fuselageDiameterTxt as text
%        str2double(get(hObject,'String')) returns contents of fuselageDiameterTxt as a double


% --- Executes during object creation, after setting all properties.
function fuselageDiameterTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fuselageDiameterTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wingARTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingARTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingARTxt as text
%        str2double(get(hObject,'String')) returns contents of wingARTxt as a double


% --- Executes during object creation, after setting all properties.
function wingARTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingARTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wingTPRTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingTPRTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingTPRTxt as text
%        str2double(get(hObject,'String')) returns contents of wingTPRTxt as a double


% --- Executes during object creation, after setting all properties.
function wingTPRTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingTPRTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wingLoadingTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingLoadingTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingLoadingTxt as text
%        str2double(get(hObject,'String')) returns contents of wingLoadingTxt as a double


% --- Executes during object creation, after setting all properties.
function wingLoadingTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingLoadingTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wingACTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingACTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingACTxt as text
%        str2double(get(hObject,'String')) returns contents of wingACTxt as a double


% --- Executes during object creation, after setting all properties.
function wingACTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingACTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wingQCSweepTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingQCSweepTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingQCSweepTxt as text
%        str2double(get(hObject,'String')) returns contents of wingQCSweepTxt as a double


% --- Executes during object creation, after setting all properties.
function wingQCSweepTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingQCSweepTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wingTCTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingTCTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingTCTxt as text
%        str2double(get(hObject,'String')) returns contents of wingTCTxt as a double


% --- Executes during object creation, after setting all properties.
function wingTCTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingTCTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hTailARTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailARTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailARTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailARTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailARTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailARTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hTailTPRTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailTPRTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailTPRTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailTPRTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailTPRTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailTPRTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hTailVRTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailVRTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailVRTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailVRTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailVRTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailVRTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hTailACTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailACTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailACTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailACTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailACTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailACTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hTailTCTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailTCTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailTCTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailTCTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailTCTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailTCTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vTargetTxt_Callback(hObject, eventdata, handles)
% hObject    handle to vTargetTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vTargetTxt as text
%        str2double(get(hObject,'String')) returns contents of vTargetTxt as a double


% --- Executes during object creation, after setting all properties.
function vTargetTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vTargetTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vTailARTxt_Callback(hObject, eventdata, handles)
% hObject    handle to vTailARTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vTailARTxt as text
%        str2double(get(hObject,'String')) returns contents of vTailARTxt as a double


% --- Executes during object creation, after setting all properties.
function vTailARTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vTailARTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vTailTPRTxt_Callback(hObject, eventdata, handles)
% hObject    handle to vTailTPRTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vTailTPRTxt as text
%        str2double(get(hObject,'String')) returns contents of vTailTPRTxt as a double


% --- Executes during object creation, after setting all properties.
function vTailTPRTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vTailTPRTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vTailVRTxt_Callback(hObject, eventdata, handles)
% hObject    handle to vTailVRTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vTailVRTxt as text
%        str2double(get(hObject,'String')) returns contents of vTailVRTxt as a double


% --- Executes during object creation, after setting all properties.
function vTailVRTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vTailVRTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vTailACTxt_Callback(hObject, eventdata, handles)
% hObject    handle to vTailACTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vTailACTxt as text
%        str2double(get(hObject,'String')) returns contents of vTailACTxt as a double


% --- Executes during object creation, after setting all properties.
function vTailACTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vTailACTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vTailTCTxt_Callback(hObject, eventdata, handles)
% hObject    handle to vTailTCTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vTailTCTxt as text
%        str2double(get(hObject,'String')) returns contents of vTailTCTxt as a double


% --- Executes during object creation, after setting all properties.
function vTailTCTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vTailTCTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function weightGuessTxt_Callback(hObject, eventdata, handles)
% hObject    handle to weightGuessTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of weightGuessTxt as text
%        str2double(get(hObject,'String')) returns contents of weightGuessTxt as a double


% --- Executes during object creation, after setting all properties.
function weightGuessTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to weightGuessTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function motorWeightTxt_Callback(hObject, eventdata, handles)
% hObject    handle to motorWeightTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of motorWeightTxt as text
%        str2double(get(hObject,'String')) returns contents of motorWeightTxt as a double


% --- Executes during object creation, after setting all properties.
function motorWeightTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to motorWeightTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function batteryWeightTxt_Callback(hObject, eventdata, handles)
% hObject    handle to batteryWeightTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of batteryWeightTxt as text
%        str2double(get(hObject,'String')) returns contents of batteryWeightTxt as a double


% --- Executes during object creation, after setting all properties.
function batteryWeightTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to batteryWeightTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function payloadWeightTxt_Callback(hObject, eventdata, handles)
% hObject    handle to payloadWeightTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of payloadWeightTxt as text
%        str2double(get(hObject,'String')) returns contents of payloadWeightTxt as a double


% --- Executes during object creation, after setting all properties.
function payloadWeightTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to payloadWeightTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vTailQCSweepTxt_Callback(hObject, eventdata, handles)
% hObject    handle to vTailQCSweepTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vTailQCSweepTxt as text
%        str2double(get(hObject,'String')) returns contents of vTailQCSweepTxt as a double


% --- Executes during object creation, after setting all properties.
function vTailQCSweepTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vTailQCSweepTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hTailQCSweepTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailQCSweepTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailQCSweepTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailQCSweepTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailQCSweepTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailQCSweepTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hTailDihedralTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailDihedralTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailDihedralTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailDihedralTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailDihedralTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailDihedralTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wingDihedralTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingDihedralTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingDihedralTxt as text
%        str2double(get(hObject,'String')) returns contents of wingDihedralTxt as a double


% --- Executes during object creation, after setting all properties.
function wingDihedralTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingDihedralTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColorpushbutton8','white');
end



function atmAltTxt_Callback(hObject, eventdata, handles)
% hObject    handle to atmAltTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of atmAltTxt as text
%        str2double(get(hObject,'String')) returns contents of atmAltTxt as a double


% --- Executes during object creation, after setting all properties.
function atmAltTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atmAltTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function vTailAirfoilTxt_Callback(hObject, eventdata, handles)
% hObject    handle to vTailAirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vTailAirfoilTxt as text
%        str2double(get(hObject,'String')) returns contents of vTailAirfoilTxt as a double


% --- Executes during object creation, after setting all properties.
function vTailAirfoilTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vTailAirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in vTailImportAirfoilButton.
function vTailImportAirfoilButton_Callback(hObject, eventdata, handles)
% hObject    handle to vTailImportAirfoilButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Airfoil
uiwait
global airfoil vTailAirfoil
vTailAirfoil = airfoil;

% --- Executes on button press in hTailImportAirfoilButton.
function hTailImportAirfoilButton_Callback(hObject, eventdata, handles)
% hObject    handle to hTailImportAirfoilButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Airfoil
uiwait
global airfoil hTailAirfoil
hTailAirfoil=airfoil;

% --- Executes on button press in wingImportAirfoilAirfoilButton.
function wingImportAirfoilAirfoilButton_Callback(hObject, eventdata, handles)
% hObject    handle to wingImportAirfoilAirfoilButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Airfoil
uiwait
global airfoil wingAirfoil
wingAirfoil=airfoil;

function wingAirfoilTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingAirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingAirfoilTxt as text
%        str2double(get(hObject,'String')) returns contents of wingAirfoilTxt as a double


% --- Executes during object creation, after setting all properties.
function wingAirfoilTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingAirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hTailAirfoilTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailAirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailAirfoilTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailAirfoilTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailAirfoilTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailAirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in solidWorksCheckBox.
function solidWorksCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to solidWorksCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of solidWorksCheckBox



function hTailIncidenceTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailIncidenceTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailIncidenceTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailIncidenceTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailIncidenceTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailIncidenceTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wingIncidenceTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingIncidenceTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingIncidenceTxt as text
%        str2double(get(hObject,'String')) returns contents of wingIncidenceTxt as a double


% --- Executes during object creation, after setting all properties.
function wingIncidenceTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingIncidenceTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xCgTxt_Callback(hObject, eventdata, handles)
% hObject    handle to xCgTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xCgTxt as text
%        str2double(get(hObject,'String')) returns contents of xCgTxt as a double


% --- Executes during object creation, after setting all properties.
function xCgTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xCgTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zCgTxt_Callback(hObject, eventdata, handles)
% hObject    handle to zCgTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zCgTxt as text
%        str2double(get(hObject,'String')) returns contents of zCgTxt as a double


% --- Executes during object creation, after setting all properties.
function zCgTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zCgTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hTailZTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailZTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailZTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailZTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailZTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailZTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wingZTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingZTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingZTxt as text
%        str2double(get(hObject,'String')) returns contents of wingZTxt as a double


% --- Executes during object creation, after setting all properties.
function wingZTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingZTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function wingNACA5AirfoilTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingNACA5AirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingNACA5AirfoilTxt as text
%        str2double(get(hObject,'String')) returns contents of wingNACA5AirfoilTxt as a double


% --- Executes during object creation, after setting all properties.
function wingNACA5AirfoilTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingNACA5AirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hTailNACA5AirfoilTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailNACA5AirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailNACA5AirfoilTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailNACA5AirfoilTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailNACA5AirfoilTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailNACA5AirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vTailNACA5AirfoilTxt_Callback(hObject, eventdata, handles)
% hObject    handle to vTailNACA5AirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vTailNACA5AirfoilTxt as text
%        str2double(get(hObject,'String')) returns contents of vTailNACA5AirfoilTxt as a double


% --- Executes during object creation, after setting all properties.
function vTailNACA5AirfoilTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vTailNACA5AirfoilTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in vTailAirfoilButtonPanel.
function vTailAirfoilButtonPanel_SelectionChangeFcn(hObject, eventdata, handles)
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject);    
switch get(hObject,'Tag')   % Get Tag of selected object
    case 'NACA4RadioButton'
        set(handles.vTailImportAirfoilButton, 'Visible', 'Off');
 
    case 'NACA5RadioButton'
        set(handles.vTailImportAirfoilButton, 'Visible', 'Off');
 
    case 'customAirfoilRadioButton'
        set(handles.vTailImportAirfoilButton, 'Visible', 'On');

    otherwise
       % Code for when there is no match.
       'otherwise'
 
end
%updates the handles structure
guidata(hObject, handles);


% --- Executes when selected object is changed in hTailAirfoilButtonPanel.
function hTailAirfoilButtonPanel_SelectionChangeFcn(hObject, eventdata, handles)
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject);    
switch get(hObject,'Tag')   % Get Tag of selected object
    case 'NACA4RadioButton'
        set(handles.hTailImportAirfoilButton, 'Visible', 'Off');
 
    case 'NACA5RadioButton'
        set(handles.hTailImportAirfoilButton, 'Visible', 'Off');
 
    case 'customAirfoilRadioButton'
        set(handles.hTailImportAirfoilButton, 'Visible', 'On');

    otherwise
       % Code for when there is no match.
       'otherwise'
 
end
%updates the handles structure
guidata(hObject, handles);


% --- Executes when selected object is changed in wingAirfoilButtonPanel.
function wingAirfoilButtonPanel_SelectionChangeFcn(hObject, eventdata, handles)
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject);    
switch get(hObject,'Tag')   % Get Tag of selected object
    case 'NACA4RadioButton'
        set(handles.wingImportAirfoilAirfoilButton, 'Visible', 'Off');
 
    case 'NACA5RadioButton'
        set(handles.wingImportAirfoilAirfoilButton, 'Visible', 'Off');
 
    case 'customAirfoilRadioButton'
        set(handles.wingImportAirfoilAirfoilButton, 'Visible', 'On');

    otherwise
       % Code for when there is no match.
       'otherwise'
 
end
%updates the handles structure
guidata(hObject, handles);


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in optimizeCheckBox.
function optimizeCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to optimizeCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optimizeCheckBox



function batteryCgTxt_Callback(hObject, eventdata, handles)
% hObject    handle to batteryCgTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of batteryCgTxt as text
%        str2double(get(hObject,'String')) returns contents of batteryCgTxt as a double


% --- Executes during object creation, after setting all properties.
function batteryCgTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to batteryCgTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
