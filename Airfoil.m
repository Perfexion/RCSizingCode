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

% Last Modified by GUIDE v2.5 19-Aug-2010 15:52:30

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
vars = evalin('base', 'who');
set(handles.listbox1, 'String', vars)

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
variablename = get_var_name(handles);
airfoil = evalin('base', [variablename ';']);
close Airfoil


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [var] = get_var_name(handles)
list_entries = get(handles.listbox1,'String');
index_selected = get(handles.listbox1,'Value');
var = list_entries{index_selected(1)};


% --- Executes on button press in updateListboxButton.
function updateListboxButton_Callback(hObject, eventdata, handles)
% hObject    handle to updateListboxButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vars = evalin('base','who');
set(handles.listbox1,'String',vars)