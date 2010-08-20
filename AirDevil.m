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

% Last Modified by GUIDE v2.5 19-Aug-2010 23:21:57

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

% --- Executes on button press in runButton.
function runButton_Callback(hObject, eventdata, handles)
% hObject    handle to runButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    

    Continue = true;
    % -- Weight Guess
    Weight = str2double(get(handles.weightGuessTxt,'String')); %lbs
    
    % -- Flight Conditions;
    VTarget =str2double(get(handles.vTargetTxt,'String')); %Target Velocity (ft/s) - For scale, 88 corresponds to 60 MPH
    
    % -- Fuselage
    FuselageLength = str2double(get(handles.fuselageLengthTxt,'String')); %ft
    FuselageDiameter = 8/12; %ft
    
    % -- Main Wing
    WingAR = str2double(get(handles.wingARTxt,'String')); %Aspect Ratio
    WingTPR = str2double(get(handles.wingTPRTxt,'String')); %Taper Rati0 
    WingLoading = str2double(get(handles.wingLoadingTxt,'String')); %Wing loading (lbs/ft^2) - Simply wing area divided by weight. In the RC aircraft world a wing loading of 1 is generally regarded as an intermediate speed aircraft.
    WingAC = str2double(get(handles.wingACTxt,'String')); %Location of aerodynamic center of wing behind nose (in ft)
    WingQCSweep = str2double(get(handles.wingQCSweepTxt,'String')); %Quarter Chord Sweep of the Wing
    WingTC = str2double(get(handles.wingTCTxt,'String')); %Wing Thickness Ratio
    WingDihadrial = str2double(get(handles.wingDihadrialTxt,'String')); %Wing Dihadrial
    
    % -- Horizontal Tail
    HTailAR = str2double(get(handles.hTailARTxt,'String')); %Aspect Ratio
    HTailTPR = str2double(get(handles.hTailTPRTxt,'String')); %Taper Ratio
    HTailVR = str2double(get(handles.hTailVRTxt,'String')); %Volume Ratio
    HTailAC = str2double(get(handles.hTailACTxt,'String')); %Location of aerodynamic center of horizontal tail behind nose (in ft)
    HTailTC = str2double(get(handles.hTailTCTxt,'String')); %Horizontal Tail Thickness Ratio
    HTailQCSweep = str2double(get(handles.hTailQCSweepTxt,'String')); %Horizontal Tail Quarter Chord Sweep
    HTailDihadrial = str2double(get(handles.hTailDihadrialTxt,'String')); %Horizontal Tail Dihadrial
    
    % -- Vertical Tail
    VTailAR = str2double(get(handles.vTailARTxt,'String')); %Aspect Ratio
    VTailTPR = str2double(get(handles.vTailTPRTxt,'String')); %Taper Ratio
    VTailVR = str2double(get(handles.vTailVRTxt,'String')); %Volume Ratio
    VTailAC = str2double(get(handles.vTailACTxt,'String')); %Location of aerodynamic center of vertical tail behind nose (in ft)
    VTailTC = str2double(get(handles.vTailTCTxt,'String')); %Vertical Tail Thickness Ratio
    VTailQCSweep = str2double(get(handles.vTailQCSweepTxt,'String')); %Vertical Tail Quarter Chord Sweep
    
    % -- Electronics
    MotorWeight = str2double(get(handles.motorWeightTxt,'String')); % Weight of Motor (lbs) - Hacker B50-12XL weighs .74 lbs
    BatteryWeight = str2double(get(handles.batteryWeightTxt,'String')); % Weight of Battery (lbs)
    PayloadWeight = str2double(get(handles.payloadWeightTxt,'String')); %Weight of Payload (lbs)

 %% Atmosphere

    %The following three equations come from David G. Hull's Book at page
    %45
    AtmAlt = str2double(get(handles.atmAltTxt, 'String')); %Atmospherica Altitude
    AtmTmp = 518.69 - 3.5662*10^-3 * AtmAlt; %Rankine
    AtmPrs = 1.1376*10^-11 * AtmTmp ^ 5.2560; %lbf/ft^2
    AtmDns = 6.6277*10^-15 * AtmTmp ^ 4.2560 * (32.1740486); %lbm/ft^3

    
    %% Begin Main Loop
    while Continue == true
       %% Calculate Wing & Tail Dimensions
       
       % -- Main Wing
       WingS = Weight/WingLoading; %Wing planform area (ft^2)
       WingB =  sqrt(WingAR*WingS); %Wing Span (ft)
       WingCR = WingS/(WingB*(1+WingTPR)/2); %Wing Root Chord (ft)
       WingCT = WingTPR*WingCR; %Wing Tip Chord (ft)
       WingCBar = WingS/WingB; %Wing Mean Aerodynamic Chord (ft)
       WingAirfoil = 'NACA 0010'; %Wing Airfoil 
       
       % -- Horizontal Tail
       HTailS = HTailVR*WingCBar*WingS/(HTailAC-WingAC); %Horizontal tail planform area (ft^2) - assumes that the wing AC is directly over (or at least very near) to the aircraft CG. This may not be a good assumption
       HTailB = sqrt(HTailAR*HTailS); %Horizontal Tail Span (ft)
       HTailCR = HTailS/(HTailB*(1+HTailTPR)/2); %Horizontal Tail Root Chord (ft)
       HTailCT = HTailTPR*HTailCR; %Horizontal Tail Tip Chord (ft)
       HTailCBar = HTailS/HTailB; %Horizontal Tail Mean Aerodynamic Chord (ft)
       HTailAirfoil = 'NACA 0010'; %Horizontal Tail Airfoil
       
       % -- Horizontal Tail
       VTailS = VTailVR*WingB*WingS/(VTailAC-WingAC); %Vertical tail planform area (ft^2) - assumes that the wing AC is directly over (or at least very near) to the aircraft CG. This may not be a good assumption
       VTailB = sqrt(VTailAR*VTailS); %Vertical Tail Span (ft)
       VTailCR = VTailS/(VTailB*(1+VTailTPR)/2); %Vertical Tail Root Chord (ft)
       VTailCT = VTailTPR*VTailCR; %Vertical Tail Tip Chord (ft)
       VTailCBar = VTailS/VTailB; %Vertical Tail Mean Aerodynamic Chord (ft)
       VTailAirfoil = 'NACA 0010'; %Vertical Tail Airfoil
       
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
               %% Organize Function Output
               AirDevilsOut = {'Fuselage Length' double(FuselageLength);
            'Fuselage Diameter' FuselageDiameter;
            'Fuselage Weight Estimate' FuselageWeight;
            'Fuselage Drag' FuselageDrag;
            'Wing AC' WingAC;
            'Wing Area' WingS;
            'Wing Span' WingB;
            'Wing Root Chord' WingCR;
            'Wing Tip Chord' WingCT;
            'Wing Mean Aerodynamic Chord' WingCBar;
            'Wing Aspect Ratio' WingAR;
            'Wing Taper Ratio' WingTPR;
            'Wing Quarter Chord Sweep' WingQCSweep;
            'Wing Dihadrial' WingDihadrial;
            'Wing Thickness Ratio' WingTC;
            'Wing Airfoil' WingAirfoil;
            'Wing Weight Estimate' WingWeight;
            'Wing Drag' WingDrag;
            'Wing Loading' WingLoading;
            'Wing Efficiency Factor' K;
            'HTail AC' HTailAC;
            'Htail Area' HTailS;
            'HTail Span' HTailB;
            'HTail Root Chord' HTailCR;
            'HTail Tip Chord' HTailCT;
            'HTail Mean Aerodynamic Chord' HTailCBar;
            'HTail Aspect Ratio' HTailAR;
            'HTail Taper Ratio' HTailTPR;
            'HTail Quarter Chord Sweep' HTailQCSweep;
            'HTail Dihadrial' HTailDihadrial;
            'HTail Thickness Ratio' HTailTC;
            'HTail Airfoil' HTailAirfoil;
            'HTail Weight Estimate' HTailWeight;
            'HTail Drag' HTailDrag;
            'HTail Volume Ratio' HTailVR;
            'VTail AC' VTailAC;
            'VTail Area' VTailS;
            'VTail Span' VTailB;
            'VTail Root Chord' VTailCR;
            'VTail Tip Chord' VTailCT;
            'VTail Mean Aerodynamic Chord' VTailCBar;
            'VTail Aspect Ratio' VTailAR;
            'VTail Taper Ratio' VTailTPR;
            'VTail Quarter Chord Sweep' VTailQCSweep;
            'VTail Thickness Ratio' VTailTC;
            'VTail Airfoil' VTailAirfoil;
            'VTail Weight Estimate' VTailWeight;
            'VTail Drag' VTailDrag;
            'VTail Volume Ratio' VTailVR;
            'Battery Weight' BatteryWeight;
            'Motor Weight' MotorWeight;
            'Payload Weight' PayloadWeight;
            'Level Flight Lift Coefficient' Cl;
            'Level Flight Drag Coefficient' Cd;
            'Zero Lift Drag Coefficient' CdZero;
            'Thrust at Target Flight Conditions' MinThrust;
            'Atmospheric Altitude' AtmAlt;
            'Atmospheric Temperature' AtmTmp;
            'Atmospheric Pressure' AtmTmp;
            'Atmospheric Density' AtmDns;
            };%Collection of Outs
        
            %%Write Output Data to a file
            fid=fopen('AirDevilsOut.csv','wt');
            [rows,cols]=size(AirDevilsOut);
            for i=1:rows
                fprintf(fid,'%s,',AirDevilsOut{i,1:end-1});
                if isnumeric(AirDevilsOut{i,end})
                    fprintf(fid,'%d\n',AirDevilsOut{i,end});
                else
                    fprintf(fid,'%s\n',AirDevilsOut{i,end});
                end
            end
            fclose(fid);
            
            assignin('base', 'AirDevilsOut', AirDevilsOut);
            
       
       end
       Weight = NewWeight;      
    end

    % --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    figure(1)
    TDView()
    
    

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



function hTailDihadrialTxt_Callback(hObject, eventdata, handles)
% hObject    handle to hTailDihadrialTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hTailDihadrialTxt as text
%        str2double(get(hObject,'String')) returns contents of hTailDihadrialTxt as a double


% --- Executes during object creation, after setting all properties.
function hTailDihadrialTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hTailDihadrialTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wingDihadrialTxt_Callback(hObject, eventdata, handles)
% hObject    handle to wingDihadrialTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wingDihadrialTxt as text
%        str2double(get(hObject,'String')) returns contents of wingDihadrialTxt as a double


% --- Executes during object creation, after setting all properties.
function wingDihadrialTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wingDihadrialTxt (see GCBO)
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
