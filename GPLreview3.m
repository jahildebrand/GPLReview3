function varargout = GPLreview3(varargin)
% GPLreview3 MATLAB code for GPLreview3.fig
%      GPLreview3, by itself, creates a new GPLreview3 or raises the existing
%      singleton*.
%
%      H = GPLreview3 returns the handle to a new GPLreview3 or the handle to
%      the existing singleton*.
%
%      GPLreview3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GPLreview3.M with the given input arguments.
%
%      GPLreview3('Property','Value',...) creates a new GPLreview3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GPLreview3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GPLreview3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GPLreview3

% Last Modified by GUIDE v2.5 08-May-2025 05:46:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GPLreview3_OpeningFcn, ...
    'gui_OutputFcn',  @GPLreview3_OutputFcn, ...
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


% --- Executes just before GPLreview3 is made visible.
function GPLreview3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GPLreview3 (see

handles.j=1;
handles.marker_count=0;
handles.reverse_vector=0;
handles.reverse_counter=0;
handles.dim_coords=0;
handles.filter=0;
handles.brightness=0.3;
handles.NextFile=0;


% Choose default command line output for GPLreview3
handles.output = hObject;

set(handles.figure1,'KeyPressFcn',@myFunction);

% Resize the figure to ensure "% completed" is visible
figPos = get(handles.figure1, 'Position'); % Get current window size
figPos(4) = figPos(4) + 30;  % Increase height to fit "% completed"
set(handles.figure1, 'Position', figPos); % Apply new size

% Force a UI refresh
drawnow;

% Save changes to handles structure
guidata(hObject, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GPLreview3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GPLreview3_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function FFTL_Callback(hObject, eventdata, handles)
% hObject    handle to FFTL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FFTL as text
%        str2double(get(hObject,'String')) returns contents of FFTL as a double

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function FFTL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FFTL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sample_freq_Callback(hObject, eventdata, handles)
% hObject    handle to sample_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sample_freq as text
%        str2double(get(hObject,'String')) returns contents of sample_freq as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sample_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start_freq_Callback(hObject, eventdata, handles)
% hObject    handle to start_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_freq as text
%        str2double(get(hObject,'String')) returns contents of start_freq as a double

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function start_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function end_freq_Callback(hObject, eventdata, handles)
% hObject    handle to end_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_freq as text
%        str2double(get(hObject,'String')) returns contents of end_freq as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function end_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_length_Callback(hObject, eventdata, handles)
% hObject    handle to plot_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_length as text
%        str2double(get(hObject,'String')) returns contents of plot_length as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function overlap_Callback(hObject, eventdata, handles)
% hObject    handle to overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of overlap as text
%        str2double(get(hObject,'String')) returns contents of overlap as a double
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function overlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_spec.
function plot_spec_Callback(hObject, eventdata, handles)
% hObject    handle to plot_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


guidata(hObject, handles);
% bt_struct = table2struct(handles.bt);  % Convert only for function call
draw_spectogram3(handles.data, ...
    str2double(get(handles.FFTL, 'String')), ...
    str2double(get(handles.overlap, 'String')), ...
    str2double(get(handles.plot_length, 'String')), ...
    str2double(get(handles.sample_freq, 'String')), ...
    str2double(get(handles.start_freq, 'String')), ...
    str2double(get(handles.end_freq, 'String')), ...
    handles.brightness, handles.yes_Bird_labels, handles.yes_A_labels, handles.markers, ...
    handles.j_included,handles.bt,handles.whiten, handles.dim_coords,handles.A_labels);

%
set(handles.start_time, 'String', ...
    ['#',num2str(handles.j_included(1)),' ',datestr(handles.bt.julian_start_time(handles.j_included(1)), 0)]); % Use the correct field name
set(handles.end_time, 'String', ...
    ['#',num2str(handles.j_included(end)),' ',datestr(handles.bt.julian_end_time(handles.j_included(end)), 0)]); % Use field name

% index = 1 + sum(handles.reverse_vector(1:end-1)); % Compute index safely
% set(handles.start_time, 'String', datestr(handles.bt.julian_start_time(index), 0)); % Use the correct field name
% index = sum(handles.reverse_vector(1:end)); % Compute index
% if index > height(handles.bt)
%     index = height(handles.bt);
% end
% set(handles.end_time, 'String', datestr(handles.bt.julian_end_time(index), 0)); % Use field name


% --- Executes on button press in pushbutton2. FORWARD
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


data=[]; handles.data=[]; handles.markers=0; length_index=0;

bt=handles.bt;
try
    handles.bt = struct2table(handles.bt);
catch
end
% handles.previous_count=handles.reverse_counter;
handles.reverse_counter=0;
handles.j_included = []; % Store values of j that contribute to the plot

while (handles.j <= height(bt))
    % test for existence of file
    if ~ismember(handles.bt.fname(handles.j), string({handles.time_file_list.name}))
        disp(['File not found: ',handles.bt.fname(handles.j)]);
        return
    else
        handles.wavefile = fullfile(handles.time_file_path,handles.bt.fname(handles.j));
        handles.wavefile = string(handles.wavefile);  % Convert to string
    end

    info = audioinfo(handles.wavefile);

    while  strcmp(fullfile(handles.time_file_path,bt.fname(handles.j)),...
            handles.wavefile)

        handles.j_included = [handles.j_included, handles.j];
        handles.reverse_counter=handles.reverse_counter+1;

        buffer=str2double(get(handles.buffer,'String'))*info.SampleRate;

        data = audioread(handles.wavefile,...
            [bt.start_time(handles.j) - buffer, ...
            bt.end_time(handles.j) + buffer]);

        % Extract only the first channel (left)
        data = data(:,1); % Use sub_data(:,2) for the right channel

        handles.data=[handles.data;data];
        handles.j=handles.j+1;
        length_index=length_index+length(data);
        handles.markers=[handles.markers,length_index];
        if (length(handles.data) > handles.hyd.detection.parm.SampleFreq*...
                str2double(handles.plot_length.String)) || ...
                handles.j > length(bt.fname)
            break
        end

    end
    if (length(handles.data) > handles.hyd.detection.parm.SampleFreq*...
            str2double(handles.plot_length.String)) || ...
            handles.j > length(bt.fname)
        if handles.j > length(bt.fname)
            handles.j = handles.j -1;
        end
        break
    end
end
set(handles.percent_completed,'String',handles.j/height(handles.bt)*100)

% handles.marker_count
% handles.marker_count=handles.marker_count+length(handles.markers);
%handles.reverse_vector
handles.reverse_vector=[handles.reverse_vector,handles.reverse_counter];

handles.plot_length_prev = get(handles.plot_length,'string');
if(handles.j==height(bt))
    set(handles.plot_length,'String',length(handles.data)/str2double(get(handles.sample_freq,'String')))
end


guidata(hObject, handles);
plot_spec_Callback(hObject, eventdata, handles)

set(handles.plot_length,'String',handles.plot_length_prev)

set(handles.wave_filename,'String',handles.wavefile);
set(handles.detection_filename,'String',handles.det_file);

pcent_false_Callback(hObject, eventdata, handles)
pcen_A_false_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
handles.data=[];

%  pushbutton2_Callback(hObject, eventdata, handles)





% --- Executes on button press in pushbutton3. REVERSE <
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Ensure reverse_vector has at least two elements before accessing end-1
if length(handles.reverse_vector) >= 2
    handles.j = handles.j - handles.reverse_vector(end) - handles.reverse_vector(end-1);
elseif length(handles.reverse_vector) == 1
    handles.j = handles.j - handles.reverse_vector(end);  % Use only the last element
else
    warning('reverse_vector is empty, decrementing j by 1 as fallback.');
    handles.j = handles.j - 1;  % Default to decrementing by 1
end
% Ensure handles.j remains a valid positive integer
handles.j = max(1, round(handles.j));

% handles.j=handles.j-handles.reverse_vector(end)-handles.reverse_vector(end-1);
% handles.marker_count=handles.marker_count-handles.reverse_vector(end);
handles.reverse_vector=handles.reverse_vector(1:end-2);


guidata(hObject, handles);
pushbutton2_Callback(hObject, eventdata, handles)

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.brightness=get(hObject,'Value');
guidata(hObject,handles);

plot_spec_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Bird_labels.
function Bird_labels_Callback(hObject, eventdata, handles)
% hObject    handle to Bird_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Bird_labels


handles.yes_Bird_labels=get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in mark. MARK SOME
function mark_Callback(hObject, eventdata, handles)
% hObject    handle to mark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles  structure with handles and user data (see GUIDATA)
coordinates=ginput(2);

[~, k2]=find(handles.markers/str2double(get(handles.overlap,'String')) > coordinates(1,1));
[~, k4]=find(handles.markers/str2double(get(handles.overlap,'String')) <= coordinates(2,1));

k5=intersect(k2-1,k4);

% indices = sum(handles.reverse_vector(1:end-1)) + k5;  % Compute indices
indices = handles.j_included(1) + k5 -1;  % Compute indices

% Assigns yes to each selected element
handles.bt.b_cross_flag(indices) = 1; % Vectorized update
Labels = num2cell(handles.bt.b_cross_flag(:));  % Convert to column vector cell array
file_path = strcat(handles.det_file_path, handles.det_file);
save(file_path, 'Labels','-append','-v7');

guidata(hObject,handles);

plot_spec_Callback(hObject, eventdata, handles)
pcent_false_Callback(hObject, eventdata, handles)
pcen_A_false_Callback(hObject, eventdata, handles)

% --- Executes on button press in whiten.
function whiten_Callback(hObject, eventdata, handles)
% hObject    handle to whiten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.whiten=get(hObject,'Value');

guidata(hObject,handles);


% --- Executes on button press in play_audio.
function play_audio_Callback(hObject, eventdata, handles)
% hObject    handle to play_audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

coordinates=ginput(2);
round(coordinates);
% Keep the crosshairs visible by adding vertical red lines
axes(gca); % Use the current active axes
hold on;
% Extract x-coordinates (time markers)
x1 = coordinates(1,1);
x2 = coordinates(2,1);
y1 = coordinates(1,2);
y2 = coordinates(2,2);

% plot([x1 x1], ylim, 'r-', 'LineWidth', 2);  % First vertical line
% plot([x2 x2], ylim, 'r-', 'LineWidth', 2);  % Second vertical line
% hold off;

[~,k2]=find(handles.markers/str2double(get(handles.overlap,'String')) > coordinates(1,1));
[~,k4]=find(handles.markers/str2double(get(handles.overlap,'String')) <= coordinates(2,1));

%added so that Amanda can see exact time of audio that's playing:
k5=intersect(k2-1,k4);
if k5(1) > 0
    % index1 = sum(handles.reverse_vector(1:end-1)) + k5(1);
    % index2 = sum(handles.reverse_vector(1:end-1)) + k5(end);
    index1 = handles.j_included(1) + k5(1) -1;
    index2 = handles.j_included(1) + k5(end) - 1;
end

% Ensure indices are within bounds
if index1 > height(handles.bt) || index2 > height(handles.bt)
    error('Index exceeds the number of elements in bt');
end

% Convert to date strings
disp(['Sound end:',datestr(handles.bt.julian_start_time(index1), 0)])
disp(['Sound start:',datestr(handles.bt.julian_start_time(index2), 0)])

if(handles.filter==1)

    plot([x1 x2], [y1 y1], 'r-', 'LineWidth', 2);  % First vertical line
    plot([x1 x2],  [y2 y2], 'r-', 'LineWidth', 2);  % Second vertical line
     plot([x1 x1], [y1 y2], 'r-', 'LineWidth', 2);  % First vertical line
    plot([x2 x2],  [y1 y2], 'r-', 'LineWidth', 2);  % Second vertical line
    hold off;
    lower_freq=coordinates(1,2)/str2double(get(handles.FFTL,'String'))*str2double(get(handles.sample_freq,'String'));
    upper_freq=coordinates(2,2)/str2double(get(handles.FFTL,'String'))*str2double(get(handles.sample_freq,'String'));

    fs=str2double(get(handles.sample_freq,'String'));          % sampling rate
    % fs
    F=[lower_freq-50 lower_freq upper_freq upper_freq+50];  % band limits
    A=[0 1 0];                % band type: 0='stop', 1='pass'
    dev=[0.0001 10^(0.1/20)-1 0.0001]; % ripple/attenuation spec
    [M,Wn,beta,typ]= kaiserord(F,A,dev,fs);  % window parameters
    b=fir1(M,Wn,typ,kaiser(M+1,beta),'noscale'); % filter design
    DATA = filter(b,1,handles.data(coordinates(1,1)*str2double(get(handles.overlap,'String')):coordinates(2,1)*str2double(get(handles.overlap,'String'))));

    handles.dim_coords = [floor(coordinates(1,2)),floor(coordinates(2,2))];

    % plot_spec_Callback(hObject, eventdata, handles)

    if(handles.speedup==1)
        DATA_resampled = resample(DATA, 1, 2); %2x speedup
        soundsc(DATA_resampled, str2double(get(handles.sample_freq,'String')));
    else
        soundsc(DATA,str2double(get(handles.sample_freq,'String')));
    end

    handles.dim_coords=0;
    % plot_spec_Callback(hObject, eventdata, handles)

else

    plot([x1 x1], ylim, 'r-', 'LineWidth', 2);  % First vertical line
    plot([x2 x2], ylim, 'r-', 'LineWidth', 2);  % Second vertical line
    hold off;
    start_idx = coordinates(1,1) * str2double(get(handles.overlap,'String'));
    end_idx   = coordinates(2,1) * str2double(get(handles.overlap,'String'));

    DATA = handles.data(start_idx:end_idx);

    if(handles.speedup==1)
        DATA_resampled = resample(DATA, 1, 2); %2x speedup
        soundsc(DATA_resampled, str2double(get(handles.sample_freq,'String')));
    else
        soundsc(DATA, str2double(get(handles.sample_freq,'String')));
    end
end

guidata(hObject,handles);



% --- Executes on button press in filter.
function filter_Callback(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filter
handles.filter=get(hObject,'Value');


guidata(hObject,handles);

% --- Executes on button press in mark_yes. MARK ALL YES
function mark_yes_Callback(hObject, eventdata, handles)
% hObject    handle to mark_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Compute the indices
% indices = sum(handles.reverse_vector(1:end-1)) + (1:length(handles.markers)-1);
indices = handles.j_included ;

% Assign value to each element
handles.bt.b_cross_flag(indices) = 1; % Vectorized update
Labels = num2cell(handles.bt.b_cross_flag(:));  % Convert to column vector cell array
file_path = strcat(handles.det_file_path, handles.det_file);
save(file_path, 'Labels','-append','-v7');

pcent_false_Callback(hObject, eventdata, handles)
pcen_A_false_Callback(hObject, eventdata, handles)
pushbutton2_Callback(hObject, eventdata, handles)


% --- Executes on button press in mark_no. MARK ALL NO
function mark_no_Callback(hObject, eventdata, handles)
% hObject    handle to mark_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% indices = sum(handles.reverse_vector(1:end-1)) + (1:length(handles.markers)-1);
indices = handles.j_included;

% Assign value 0 to the field 'b_cross_flag' for all selected elements
handles.bt.b_cross_flag(indices) = 0; % Vectorized update

Labels = num2cell(handles.bt.b_cross_flag(:));  % Convert to column vector cell array
file_path = strcat(handles.det_file_path, handles.det_file);
save(file_path, 'Labels','-append','-v7');

guidata(hObject,handles);

pcent_false_Callback(hObject, eventdata, handles)
pcen_A_false_Callback(hObject, eventdata, handles)
pushbutton2_Callback(hObject, eventdata, handles)

% --- Executes on key press with focus on filter and none of its controls.
function filter_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


%myFunction(src,evnt)
%function myFunction(src,evnt, handles)
function myFunction(src, eventdata, handles, hObject)

%this function takes in two inputs by default

%src is the gui figure
%evnt is the keypress information

%this line brings the handles structures into the local workspace
%now we can use handles.cats in this subfunction!

handles = guidata(src);
hObject = handles.output;
%switch evnt.Key
switch eventdata.Key
    case 'leftarrow'
        pushbutton3_Callback(hObject, eventdata, handles)
    case 'rightarrow'
        pushbutton2_Callback(hObject, eventdata, handles)
    case 'y'
        mark_yes_Callback(hObject, eventdata, handles)
    case 'n'
        mark_no_Callback(hObject, eventdata, handles)
    case 'a'
        pushbutton10_Callback(hObject, eventdata, handles)
    case 's'
        pushbutton12_Callback(hObject, eventdata, handles)
    case 'm'
        mark_Callback(hObject, eventdata, handles)
    case '0'
        mark_subset_no_Callback(hObject, eventdata, handles)

    case 'escape'
end



% --------------------------------------------------------------------
function uipushtool3_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName,FilterIndex] = uigetfile('mat');

addpath(PathName);
handles.det_file_list=dir(PathName);
handles.det_file=FileName;
handles.det_file_path=PathName;

% load([PathName,handles.det_file]);
% Load the selected file
loadedData = load([PathName, handles.det_file]);
% Determine the format of the loaded data
% The first time v3 data are read, they are in 'hyd"
% after editing they are in the v1 "bt" format.
if isfield(loadedData, 'bt')
    disp('Loaded data is in v1 format.');
    handles.format_version = 'v1';
end
if isfield(loadedData, 'hyd')
    disp('Loaded data is in v3 format.');
    handles.hyd = loadedData.hyd;
    handles.format_version = 'v3';
end
if isfield(loadedData, 'Times')
    disp('Loaded data is in v4 Alba format.');
    handles.hyd.detection.calls = loadedData.Times;
    labels = cell2mat(loadedData.Labels);
    numCalls = numel(handles.hyd.detection.calls);
    numLabels = numel(labels);
    if numCalls == numLabels
        for i = 1:numCalls
            handles.hyd.detection.calls(i).boundary_cross_flag = labels(i);
        end
    else
        error('Mismatch: calls has %d elements, but Labels has %d elements', numCalls, numLabels);
    end

    handles.hyd.detection.parm = loadedData.parm;
    handles.format_version = 'v4';
end
if ~isfield(loadedData, 'hyd') && ~isfield(loadedData, 'bt')  ...
        && ~isfield(loadedData, 'Times')
    disp('Loaded data unknown format');
    return
end

if ~isfield(loadedData, 'bt')
    disp('GPL_v3 data format input');
    handles.sample_freq.String = num2str(handles.hyd.detection.parm.SampleFreq);
    handles.start_freq.String = num2str(handles.hyd.detection.parm.freq_lo);
    handles.end_freq.String = num2str(handles.hyd.detection.parm.freq_hi);
    handles.FFTL.String = num2str(handles.hyd.detection.parm.fftl);
    handles.overlap.String = num2str(round(handles.hyd.detection.parm.fftOverlap));

    % Extract start and end times from struct array correctly
    start_times = arrayfun(@(x) x.start_time, handles.hyd.detection.calls);
    end_times = arrayfun(@(x) x.end_time, handles.hyd.detection.calls);
    b_cross_flags = arrayfun(@(x) x.boundary_cross_flag, handles.hyd.detection.calls);
    julian_start_times = arrayfun(@(x) x.julian_start_time, handles.hyd.detection.calls);
    julian_end_times = arrayfun(@(x) x.julian_end_time, handles.hyd.detection.calls);
    fnames = string(arrayfun(@(x) x.fname, handles.hyd.detection.calls, 'UniformOutput', false));

    % Ensure all numeric variables are column vectors
    start_times = round(start_times(:));
    end_times = round(end_times(:));
    b_cross_flags = b_cross_flags(:);
    julian_start_times = julian_start_times(:);
    julian_end_times = julian_end_times(:);
    fnames = fnames(:); % Ensure fnames is also a column vector

    % Create a table for handles.bt
    handles.bt = struct2table(struct(...
        'start_time', num2cell(start_times), ...
        'end_time', num2cell(end_times), ...
        'b_cross_flag', num2cell(b_cross_flags), ...
        'julian_start_time', num2cell(julian_start_times), ...
        'julian_end_time', num2cell(julian_end_times), ...
        'fname', cellstr(fnames))); % Convert struct to table
    %
else
    disp('GPL_v1 data format');
    handles.bt = struct2table(loadedData.bt); %
    if isfield(loadedData, 'hyd') % using v3 but after initial edit
        handles.sample_freq.String = num2str(handles.hyd.detection.parm.SampleFreq);
        handles.start_freq.String = num2str(handles.hyd.detection.parm.freq_lo);
        handles.end_freq.String = num2str(handles.hyd.detection.parm.freq_hi);
        handles.FFTL.String = num2str(handles.hyd.detection.parm.fftl);
        handles.overlap.String = num2str(round(handles.hyd.detection.parm.fftOverlap));
    else
        handles.sample_freq.String = num2str(loadedData.parm.sample_freq);
        handles.start_freq.String = num2str(loadedData.parm.freq_lo);
        handles.end_freq.String = num2str(loadedData.parm.freq_hi);
        handles.FFTL.String = num2str(loadedData.parm.fftl);
        handles.overlap.String = num2str(loadedData.parm.skip);
    end

end

%JAH new code to split handles.bt into enteries of 1-2 sec length
% Determine sampling rate
SampleFreq = str2double(handles.sample_freq.String);
one_sec_samples = SampleFreq; % 44100 samples per second

% Test total duration before splitting
duration = handles.bt.end_time - handles.bt.start_time;
tot_time = sum(duration) / SampleFreq;
% disp(['Total time Initial (sec): ', num2str(tot_time)])
% disp(['Initial # detections): ', num2str(height(handles.bt))])

% Save original for checking (optional)
handles.bt_original = handles.bt;

% Preallocate rough guess
estimated_total_rows = ceil(height(handles.bt) * 2);

start_time_array = zeros(estimated_total_rows,1);
end_time_array = zeros(estimated_total_rows,1);
b_cross_flag_array = zeros(estimated_total_rows,1);
julian_start_time_array = zeros(estimated_total_rows,1);
julian_end_time_array = zeros(estimated_total_rows,1);
fname_array = strings(estimated_total_rows,1);

row_idx = 1; % Initialize

for i = 1:height(handles.bt)
    % === Compute basic info ===
    start_sample = handles.bt.start_time(i);
    end_sample = handles.bt.end_time(i);
    total_samples = end_sample - start_sample + 1;

    num_full_chunks = floor(total_samples / one_sec_samples);
    leftover_samples = total_samples - num_full_chunks * one_sec_samples;

    current_start = start_sample;

    % === Create full 1-second chunks ===
    for j = 1:num_full_chunks
        new_start = current_start;
        new_end = new_start + one_sec_samples - 1;

        % Interpolate Julian times
        frac_start = (new_start - start_sample) / (end_sample - start_sample);
        frac_end = (new_end - start_sample) / (end_sample - start_sample);

        new_julian_start = handles.bt.julian_start_time(i) + frac_start * ...
            (handles.bt.julian_end_time(i) - handles.bt.julian_start_time(i));
        new_julian_end = handles.bt.julian_start_time(i) + frac_end * ...
            (handles.bt.julian_end_time(i) - handles.bt.julian_start_time(i));

        % Save this chunk
        start_time_array(row_idx) = new_start;
        end_time_array(row_idx) = new_end;
        b_cross_flag_array(row_idx) = handles.bt.b_cross_flag(i);
        julian_start_time_array(row_idx) = new_julian_start;
        julian_end_time_array(row_idx) = new_julian_end;
        fname_array(row_idx) = handles.bt.fname(i);

        current_start = new_end + 1; % move to next
        row_idx = row_idx + 1;
    end

    % === Handle leftover ===
    if leftover_samples > 0
        if (leftover_samples / SampleFreq) < 0.1 && row_idx > 1
            % Merge leftover into previous chunk
            end_time_array(row_idx-1) = end_sample;

            % Correct Julian end time
            julian_end_time_array(row_idx-1) = handles.bt.julian_end_time(i);
        else
            % Create a final small chunk
            new_start = current_start;
            new_end = end_sample;

            frac_start = (new_start - start_sample) / (end_sample - start_sample);
            frac_end = (new_end - start_sample) / (end_sample - start_sample);

            new_julian_start = handles.bt.julian_start_time(i) + frac_start * ...
                (handles.bt.julian_end_time(i) - handles.bt.julian_start_time(i));
            new_julian_end = handles.bt.julian_start_time(i) + frac_end * ...
                (handles.bt.julian_end_time(i) - handles.bt.julian_start_time(i));

            start_time_array(row_idx) = new_start;
            end_time_array(row_idx) = new_end;
            b_cross_flag_array(row_idx) = handles.bt.b_cross_flag(i);
            julian_start_time_array(row_idx) = new_julian_start;
            julian_end_time_array(row_idx) = new_julian_end;
            fname_array(row_idx) = handles.bt.fname(i);

            row_idx = row_idx + 1;
        end
    end
end

% === Trim arrays to final size ===
start_time_array = start_time_array(1:row_idx-1);
end_time_array = end_time_array(1:row_idx-1);
b_cross_flag_array = b_cross_flag_array(1:row_idx-1);
julian_start_time_array = julian_start_time_array(1:row_idx-1);
julian_end_time_array = julian_end_time_array(1:row_idx-1);
fname_array = fname_array(1:row_idx-1);

% === Build final table ===
handles.bt = table(start_time_array, end_time_array, ...
    b_cross_flag_array, julian_start_time_array, julian_end_time_array, fname_array, ...
    'VariableNames', {'start_time', 'end_time', 'b_cross_flag', 'julian_start_time', 'julian_end_time', 'fname'});

% === Test total duration after splitting ===
duration = handles.bt.end_time - handles.bt.start_time;
tot_time = sum(duration) / (SampleFreq);
% disp(['Total time Final (sec): ', num2str(tot_time)])
% disp(['Final # detections): ', num2str(height(handles.bt))])

handles.total_length.String = num2str(round(tot_time/60,3));
% === Extra Validation ===
% initial_samples = sum(handles.bt_original.end_time - handles.bt_original.start_time + 1);
% final_samples = sum(handles.bt.end_time - handles.bt.start_time + 1);
% disp(['Initial total samples: ', num2str(initial_samples)])
% disp(['Final total samples:   ', num2str(final_samples)])

% Store file
Times = table2struct(handles.bt);
Labels = num2cell(handles.bt.b_cross_flag(:));  % Convert to column vector cell array
file_path = strcat(handles.det_file_path, handles.det_file);
save(file_path, 'Times','Labels','-append');

handles.FileName = FileName;
handles.yes_Bird_labels = 0;
handles.yes_A_labels = 0;
handles.A_labels = [];
% Update GUI data
guidata(hObject, handles);


% --------------------------------------------------------------------
function uipushtool5_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName, FilterIndex] = uigetfile( '*.*', 'Select a Wav File');

addpath(PathName);
handles.time_file_path = PathName;
handles.time_file_list = dir(fullfile(PathName, '*.wav'));
handles.wavefile=FileName;

info = audioinfo(fullfile(handles.time_file_path,handles.wavefile));
if handles.sample_freq.String ~= num2str(info.SampleRate)
    disp(['handles.sample_freq.String',handles.sample_freq.String, but ',...' ...
        'infor.SampleRate',num2str(info.SampleRate)])
    return
end
guidata(hObject,handles);



% --- Executes on button press in pushbutton10. previous AIRPLANE
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Extract important variables
% Extract important variables
A_labels = handles.A_labels;
yes_A_labels = handles.yes_A_labels;
j_included = handles.j_included;

% Proceed only if airplane labels are active
if yes_A_labels == 1

    % Find previous plane detection index before current j_included(1)
    previousPlaneRel = find(A_labels(1:j_included(1)-1) > 0, 1, 'last');

    if ~isempty(previousPlaneRel)
        % Now walk backward to find where the plane block started
        idx = previousPlaneRel;

        while idx > 1 && A_labels(idx-1) > 0
            idx = idx - 1;
        end

        % Also update handles.j to the starting index
        handles.j = idx;

        % Save changes
        guidata(hObject, handles);

        % Refresh the plot
        pushbutton2_Callback(hObject, eventdata, handles);

        % Adjust axes if needed (optional: shift view here)

    else
        disp('No previous airplane detection found.');
    end

end



% --- Executes on button press in 2x speedup.
function speedup_Callback(hObject, eventdata, handles)
% hObject    handle to speedup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.speedup=get(hObject,'Value');
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of speedup

function marker_number_Callback(hObject, eventdata, handles)
% hObject    handle to marker_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of marker_number as text
handles.marker_number = str2double(get(hObject,'String'));%
if handles.marker_number > height(handles.bt)
    handles.marker_number =  height(handles.bt)  - 30;
end
jump = handles.j - handles.marker_number;
handles.j = handles.j - jump;

% handles.marker_count = handles.marker_count - jump;

guidata(hObject,handles);
pushbutton2_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function marker_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to marker_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.marker_number=get(hObject,'Value');

guidata(hObject,handles);



function buffer_Callback(hObject, eventdata, handles)
% hObject    handle to buffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of buffer as text
%        str2double(get(hObject,'String')) returns contents of buffer as a double


% --- Executes during object creation, after setting all properties.
function buffer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to buffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mark_subset_no. MARK SOME NO
function mark_subset_no_Callback(hObject, eventdata, handles)
% hObject    handle to mark_subset_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

coordinates=ginput(2);

[~, k2]=find(handles.markers/str2double(get(handles.overlap,'String')) > coordinates(1,1));
[~, k4]=find(handles.markers/str2double(get(handles.overlap,'String')) <= coordinates(2,1));


k5=intersect(k2-1,k4);

% Compute the indices
% indices = sum(handles.reverse_vector(1:end-1)) + k5;
indices = handles.j_included(1) + k5 -1;

% Assign the value '0' to the 'b_cross_flag' field for all selected elements
handles.bt.b_cross_flag(indices) = 0; % Vectorized update
Labels = num2cell(handles.bt.b_cross_flag(:));  % Convert to column vector cell array
file_path = strcat(handles.det_file_path, handles.det_file);
save(file_path, 'Labels','-append','-v7');

guidata(hObject,handles);

plot_spec_Callback(hObject, eventdata, handles)
pcent_false_Callback(hObject, eventdata, handles)
pcen_A_false_Callback(hObject, eventdata, handles)

% --- Executes on button press in pushbutton12, next Airplane >
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Extract important variables
A_labels = handles.A_labels;
yes_A_labels = handles.yes_A_labels;
j_included = handles.j_included;

% Proceed only if airplane labels are active
if yes_A_labels == 1

    % Find next plane detection index greater than j_included
    nextPlaneRel = find(A_labels(j_included(end)+1:end) > 0, 1, 'first');

    if ~isempty(nextPlaneRel)
        nextPlane = nextPlaneRel + j_included(end); % Correct for relative indexing
    else
        nextPlane = []; % No next plane found
    end
    if ~isempty(nextPlane)
        % Update handles.j to the next plane detection
        handles.j = nextPlane;

        % Update handles structure
        guidata(hObject, handles);

        % Call pushbutton2 callback to update the GUI
        pushbutton2_Callback(hObject, eventdata, handles);
    else
        disp('No next airplane detection found.');
    end

end


% --- Executes on button press in Airplane_labels.
function Airplane_labels_Callback(hObject, eventdata, handles)
% hObject    handle to Airplane_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.yes_A_labels = get(hObject,'Value'); %returns toggle state of Airplane_labels
guidata(hObject,handles);

if handles.yes_A_labels ==  1

    %test if variable A_labels exists
    if exist('handles.A_labels', 'var')
        disp('Variable "A_labels" already exists in workspace.');
    else
        [AFileName, APathName] = uigetfile('*.mat', 'Select a MAT file');

        % Save file info to handles
        handles.a_file = AFileName;
        handles.a_file_path = APathName;

        % Load the selected MAT file
        loadedAData = load(fullfile(APathName, AFileName));
    end
    % Bird detection start and end times
    bird_start = handles.bt.julian_start_time;
    bird_end = handles.bt.julian_end_time;

    % Airplane detection start and end times
    plane_start = loadedAData.shipTimes(:,1);
    plane_end = loadedAData.shipTimes(:,2);

    % Map ship labels to numeric plane types
    % 'jet' -> 1, 'prop' -> 2
    num_planes = length(loadedAData.shipLabels);
    plane_type = zeros(num_planes,1); % preallocate

    for j = 1:num_planes
        if strcmpi(loadedAData.shipLabels{j}, 'jet')
            plane_type(j) = 1;
        elseif strcmpi(loadedAData.shipLabels{j}, 'prop')
            plane_type(j) = 2;
        else
            plane_type(j) = NaN; % unknown type
        end
    end

    % Initialize labels array
    handles.A_labels = zeros(size(bird_start)); % now it stores 0, 1, or 2

    % Check each bird detection
    for i = 1:length(bird_start)
        % Find which planes overlap
        overlapping_planes = find( (bird_start(i) <= plane_end) & (bird_end(i) >= plane_start) );

        if ~isempty(overlapping_planes)
            % If multiple planes overlap, you can decide what to do.
            % For now, take the first overlapping plane's type.
            handles.A_labels(i) = plane_type(overlapping_planes(1));
        else
            handles.A_labels(i) = 0; % no plane overlap
        end
    end

    duration = plane_end - plane_start;
    tot_time = sum(duration) *24*60*60;
    % disp(['Total time Plane (sec): ', num2str(tot_time)])
    handles.A_length.String = num2str(round(tot_time,3));
    guidata(hObject, handles);

end



function total_length_Callback(hObject, eventdata, handles)
% hObject    handle to total_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of total_length as text
%        str2double(get(hObject,'String')) returns contents of total_length as a double


% --- Executes during object creation, after setting all properties.
function total_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to total_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pcent_false_Callback(hObject, eventdata, handles)
% hObject    handle to pcent_false (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SampleFreq = str2double(handles.sample_freq.String); % Get sample frequency

% Find indices of detections that are True
itrue = find(handles.bt.b_cross_flag);

% Calculate total time for True detections
duration_true = handles.bt.end_time(itrue) - handles.bt.start_time(itrue);
tot_true_sec = sum(duration_true) / SampleFreq;
% disp(['Total time True (sec): ', num2str(tot_true_sec)])

% Calculate total time overall
duration_total = handles.bt.end_time - handles.bt.start_time;
tot_time_sec = sum(duration_total) / SampleFreq;

% Numeric calculation
tot_false_sec = tot_time_sec - tot_true_sec;
pcent_false = (tot_false_sec / tot_time_sec) * 100;
pcent_false = round(pcent_false, 3);

% Save the number
handles = guidata(hObject);  % Refresh

% Update GUI box (the edit text)
set(handles.pcent_false, 'String', sprintf('%.3f', pcent_false));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pcent_false_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcent_false (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A_length_Callback(hObject, eventdata, handles)
% hObject    handle to A_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A_length as text
%        str2double(get(hObject,'String')) returns contents of A_length as a double


% --- Executes during object creation, after setting all properties.
function A_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pcen_A_false_Callback(hObject, eventdata, handles)
% hObject    handle to pcen_A_false (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pcen_A_false as text
%        str2double(get(hObject,'String')) returns contents of pcen_A_false as a double

SampleFreq = str2double(handles.sample_freq.String); % Get sample frequency

% Find indices of detections that are True
if handles.yes_A_labels == 1
    iair = find(handles.A_labels > 0);
    itrue = find(handles.bt.b_cross_flag & handles.A_labels > 0);

    % Calculate total time for True detections with airplanes
    duration_air = handles.bt.end_time(iair) - handles.bt.start_time(iair);
    tot_air_sec = sum(duration_air) / SampleFreq;
    % disp(['Total Airplane (sec): ', num2str(tot_air_sec)])

    % Calculate total bird and airplane true
    duration_true = handles.bt.end_time(itrue) - handles.bt.start_time(itrue);
    tot_true_sec = sum(duration_true) / SampleFreq;

    % Numeric calculation
    tot_false_sec = tot_air_sec - tot_true_sec;
    pcen_A_false = (tot_false_sec / tot_air_sec) * 100;
    pcen_A_false = round(pcen_A_false, 3);

    % Update GUI box (the edit text)
    set(handles.pcen_A_false, 'String', sprintf('%.3f', pcen_A_false));
    % Save the number
    handles = guidata(hObject);  % Refresh
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pcen_A_false_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcen_A_false (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reset_true. RESET ALL BIRD DET to be TRUE
function reset_true_Callback(hObject, eventdata, handles)
% hObject    handle to reset_true (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the state of the checkbox
handles.yes_reset = get(hObject,'Value'); % returns toggle state
guidata(hObject,handles);

if handles.yes_reset == 1
    % Ask user for confirmation
    choice = questdlg('Do you want to Reset Bird Detections to True and save?', ...
        'Confirm Update', ...
        'Yes','No','No');  % Default = 'No'

    % Handle response
    switch choice
        case 'Yes'
            % === Execute the reset ===
            handles.bt.b_cross_flag(:) = 1;  % Set all flags to true
            Labels = num2cell(handles.bt.b_cross_flag(:));  % Convert to cell array
            file_path = strcat(handles.det_file_path, handles.det_file);
            save(file_path, 'Labels', '-append', '-v7');

            disp('Reset and save completed.');
            guidata(hObject, handles);

            % Refresh the GUI
            pcent_false_Callback(hObject, eventdata, handles)
            pcen_A_false_Callback(hObject, eventdata, handles)
            pushbutton2_Callback(hObject, eventdata, handles);

        case 'No'
            disp('Reset canceled by user.');
    end

    % Uncheck the checkbox after operation
    set(hObject, 'Value', 0);
    handles.yes_reset = 0;
    guidata(hObject, handles);
end
