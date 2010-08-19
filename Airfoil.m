function varargout = Airfoil(varargin)
% AIRFOIL M-file for Airfoil.fig
%      AIRFOIL, by itself, creates a new AIRFOIL or raises the existing
%      singleton*.
%
%      H = AIRFOIL returns the handle to a new AIRFOIL or the handle to
%      the existing singleton*.
%
%      AIRFOIL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AIRFOIL.M with the given input arguments.
%
%      AIRFOIL('Property','Value',...) creates a new AIRFOIL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Airfoil_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Airfoil_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Airfoil

% Last Modified by GUIDE v2.5 18-Aug-2010 17:16:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Airfoil_OpeningFcn, ...
                   'gui_OutputFcn',  @Airfoil_OutputFcn, ...
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


% --- Executes just before Airfoil is made visible.
function Airfoil_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Airfoil (see VARARGIN)

% Choose default command line output for Airfoil
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Airfoil wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Airfoil_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in doneButton.
function doneButton_Callback(hObject, eventdata, handles)
% hObject    handle to doneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global airfoil
airfoil=get(handles.airfoilXYTable, 'Data');
close Airfoil

function rowsTxt_Callback(hObject, eventdata, handles)
% hObject    handle to rowsTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rowsTxt as text
%        str2double(get(hObject,'String')) returns contents of rowsTxt as a double


% --- Executes during object creation, after setting all properties.
function rowsTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rowsTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in updateTableButton.
function updateTableButton_Callback(hObject, eventdata, handles)
% hObject    handle to updateTableButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rows = str2double(get(handles.rowsTxt,'string'));
num_elem = cell(rows,2);
num_elem(:,:) = {''};
set(handles.airfoilXYTable,'Data',num_elem);
set(handles.airfoilXYTable,'ColumnEditable',true(1,2));
set(handles.airfoilXYTable);
