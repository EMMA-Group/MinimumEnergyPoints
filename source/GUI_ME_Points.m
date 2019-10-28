function varargout = GUI_ME_Points(varargin)
%GUI_ME_POINTS MATLAB code file for GUI_ME_Points.fig
%      GUI_ME_POINTS, by itself, creates a new GUI_ME_POINTS or raises the existing
%      singleton*.
%
%      H = GUI_ME_POINTS returns the handle to a new GUI_ME_POINTS or the handle to
%      the existing singleton*.
%
%      GUI_ME_POINTS('Property','Value',...) creates a new GUI_ME_POINTS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to GUI_ME_Points_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUI_ME_POINTS('CALLBACK') and GUI_ME_POINTS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUI_ME_POINTS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES, GENERATEPOINTS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% GUI_ME_Points.m
% Copyright (C) 2018, Felix Fritzen and Oliver Kunc
% All rights reserved.
%
% This source code is licensed under the BSD 3-Clause License found in the
% LICENSE file in the root directory of this source tree.
%
% This program employs a modified version of the softwares
%
% 1) Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox.
%    Release 1.10 2005-06-26
%
%    written by Paul Leopardi for the University of New South Wales.
%
%    See COPYING in the subfolder eq_sphere_partitions for
%    licensing information regarding this software.
%
%    See CHANGELOG in the subfolder eq_sphere_partitions for
%    a concise list of changes that were made to the original code.
%
% 2) UIGETVAR_PUTVAR
%    Release 1.0 2010-02-19
%
%    written by John D'Errico
%
%    See uigetvar.m in the subfolder uigetvar_putvar for more information.
%
%    See license.txt in the subfolder uigetvar_putvar for licensing
%    information regarding this software.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This software package is related to the research article
%
% Oliver Kunc and Felix Fritzen: ''
% JOURNAL NAME, Number/Volume, p. XX-YY, 2019
% DOI   ...
% URL   dx.doi.org/...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ME_Points_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ME_Points_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before GUI_ME_Points is made visible.
function GUI_ME_Points_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

addpath(genpath('..'))

% Choose default command line output for GUI_ME_Points
handles.output = hObject;

% Application specific data ...
handles.execution_directory = pwd;         % store the working directory of the GUI. make sure this is still the working directory whenever a button is pushed.
% ... input parameters
handles.parameters = struct;
handles.parameters.D = [];
handles.parameters.N = [];
handles.parameters.sym_flag = [];
handles.parameters.energy_index = [];
handles.parameters.n_cycle_energy = [];
handles.parameters.n_it_energy = [];
handles.parameters.n_cycle_gradient = [];
handles.parameters.n_it_gradient = [];
handles.parameters.Xstart = [];            % matrix containing the initial point coordinates as columns
handles.parameters.Xstart_origin = [];
handles.parameters.Xstart_name = [];
handles.parameters.creation_date = [];
% ... results
handles.results = struct;
handles.results.X = [];                    % matrix containing the point coordinates as columns, final or temporary
handles.results.NN_MinMeanMax = [];        % nearest neighbor: min, mean, and max value
handles.results.meshnorm = [];             % the mesh norm of X
handles.results.meshratio = [];            % the mesh ratio of X, = mesh norm * 2 / min(NN)
handles.results.largestgap = [];           % the coordinate of the largest gap
handles.results.gamma_POU_value = [];      % value of gamma, optimized to fit constant one function ('partition of unity')
handles.results.gamma_POU_RMSE = [];       % rooted mean square error of SBF interpolation to fit constant one function

% Disable buttons that should not be pushed at first ...
set(handles.pushbutton_plotEDF,'enable','off');
set(handles.pushbutton_export_workspace,'enable','off');
set(handles.pushbutton_export_path,'enable','off');
set(handles.pushbutton_export_files,'enable','off');
set(handles.pushbutton_plotpoints3d,'enable','off');
set(handles.pushbutton_gamma_POU,'enable','off');
% ... and also edit fields
set(handles.edit_export,'enable','off');
set(handles.edit_gamma_min,'enable','off');
set(handles.edit_gamma_max,'enable','off');
set(handles.edit_gamma_steps,'enable','off');
set(handles.edit_gamma_Nvali,'enable','off');

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_ME_Points_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_generate_points.
function pushbutton_generate_points_Callback(hObject, eventdata, handles)

    % Disable all buttons
    SetEnableAllPushbuttons('off', handles, gcbf);

    % Store current working directory as a string. Change to execution
    % directory.
    current_directory = pwd;
    cd(handles.execution_directory);

    % Copy user parameters to parameters struct: part I
    handles.parameters.D = CheckInteger(str2double(get(handles.edit_D,'string')),2);
    handles.parameters.N = CheckInteger(str2double(get(handles.edit_N,'string')),2);
    handles.parameters.sym_flag = get(handles.checkbox_sym,'value');
        raw_val = get(handles.popupmenu_energy,'value');
        if raw_val==3
            idx = 0.5;
        else
            idx = raw_val-3;
        end
    handles.parameters.energy_index = idx;
    handles.parameters.n_cycle_energy = CheckInteger(str2double(get(handles.edit_energy_cycles,'string')),0);
    handles.parameters.n_it_energy = CheckInteger(str2double(get(handles.edit_energy_iterations,'string')),0);
    handles.parameters.n_cycle_gradient = CheckInteger(str2double(get(handles.edit_gradient_cycles,'string')),0);
    handles.parameters.n_it_gradient = CheckInteger(str2double(get(handles.edit_gradient_iterations,'string')),0);

    % Copy user parameters to easy-to-read local variables
    D = handles.parameters.D;
    N = handles.parameters.N;
    sym_flag         = handles.parameters.sym_flag;
    energy_index     = handles.parameters.energy_index;
    n_cycle_energy   = handles.parameters.n_cycle_energy;
    n_it_energy      = handles.parameters.n_it_energy;
    n_cycle_gradient = handles.parameters.n_cycle_gradient;
    n_it_gradient    = handles.parameters.n_it_gradient;

    % Copy user parameters to parameter struct: part II
    % handles.parameters.Xstart = [];               % matrix containing the initial point coordinates as columns
    if get(handles.radiobutton_init_equalarea,'value') == 1
        handles.parameters.Xstart_origin = 'Equal Area';
        handles.parameters.Xstart_name = 'Equal Area';
        if( sym_flag == 1 && mod(N,2)==0 )
            fprintf('# new initialization\n');
            Xstart = eq_point_set(D-1,2*N);
            handles.parameters.Xstart = RenormalizeColumns( 0.002*RandomPoints(D,N) +  Xstart(:, 1:N) );
%             Xstart = eq_point_set(D-1,N+1);
%             handles.parameters.Xstart = RenormalizeColumns( 0.001*RandomPoints(D,N) +  Xstart(:, 1:N) );
        else
            handles.parameters.Xstart = RenormalizeColumns( 0.001*RandomPoints(D,N) + eq_point_set(D-1,N) );
        end
    elseif get(handles.radiobutton_init_random,'value') == 1
        handles.parameters.Xstart_origin = 'random';
        handles.parameters.Xstart_name = 'random';
        handles.parameters.Xstart = RandomPoints(D,N);
    elseif get(handles.radiobutton_init_workspace,'value') == 1
        handles.parameters.Xstart_origin = 'workspace';
        % handles.parameters.Xstart_name was set by corresponding pushbutton_export_workspace callback
    elseif get(handles.radiobutton_init_matfile,'value') == 1
        handles.parameters.Xstart_origin = 'matfile';
        % handles.parameters.Xstart_name was set by corresponding pushbutton_export_workspace callback
    elseif get(handles.radiobutton_init_txtfile,'value') == 1
        handles.parameters.Xstart_origin = 'txtfile';
        % handles.parameters.Xstart_name was set by corresponding pushbutton_export_workspace callback
    end
    Xstart = handles.parameters.Xstart;
    handles.parameters.creation_date = datestr(now);

    % TODO introduce fixpoint functionality
    disp('------------------------------------')
    disp('Use fminunc to minimize energy I ...')
    tic
    X = MinimizeEnergy( D, N, n_it_energy,  n_cycle_energy, Xstart, energy_index, sym_flag );
    toc

    disp('------------------------------------')
    disp('Use lsqnonlin to minimize gradient g ...')
    tic
    X = MinimizeGradient( D, N, n_it_gradient, n_cycle_gradient, X, energy_index, sym_flag );
    toc

%     % Plot empirical distribution functions
%     disp('Plotting empirical distribution functions')
%     tic
%     cla(handles.axes1)
%     PlotEDF(X, sym_flag, handles.axes1);
%     toc

    % Nearest neighbor (NN) and mesh statistics
    [NN_xi_MinMeanMax, meshnorm, meshratio, largestgap] = NN_and_mesh_statistics(X,sym_flag);

    % Display mesh and NN statistics
    set(handles.text_meshstat_norm,'string',        sprintf('%7.5f',meshnorm));
    set(handles.text_meshstat_separation,'string',  sprintf('%7.5f',0.5*NN_xi_MinMeanMax(1)));
    set(handles.text_meshstat_ratio,'string',       sprintf('%7.5f',meshratio));
    set(handles.text_NN_min,'string',               sprintf('%7.5f',NN_xi_MinMeanMax(1)));
    set(handles.text_NN_mean,'string',              sprintf('%7.5f',NN_xi_MinMeanMax(2)));
    set(handles.text_NN_max,'string',               sprintf('%7.5f',NN_xi_MinMeanMax(3)));

    % Store results in handles structure
    handles.results.X = X;
    handles.results.NN_MinMeanMax = NN_xi_MinMeanMax;
    handles.results.meshnorm = meshnorm;
    handles.results.meshratio = meshratio;
    handles.results.largestgap = largestgap;

    % Enable all edit fields
    set(handles.edit_export,'enable','on');
    set(handles.edit_gamma_min,'enable','on');
    set(handles.edit_gamma_max,'enable','on');
    set(handles.edit_gamma_steps,'enable','on');
    set(handles.edit_gamma_Nvali,'enable','on');

    % Update handles structure with data
    guidata(gcbf, handles);

    % Enable all pushbuttons
    SetEnableAllPushbuttons('on', handles, gcbf);

    disp(['Finished gegnerating ', num2str(N), ' points in ', ...
        num2str(D), ' dimensions.']);

    % Change back to directory that was active on click time
    cd(current_directory);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_sym.
function checkbox_sym_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sym (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_sym



function edit_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_N as text
%        str2double(get(hObject,'String')) returns contents of edit_N as a double


% --- Executes during object creation, after setting all properties.
function edit_N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D_Callback(hObject, eventdata, handles)
% hObject    handle to edit_D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_D as text
%        str2double(get(hObject,'String')) returns contents of edit_D as a double


% --- Executes during object creation, after setting all properties.
function edit_D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_energy.
function popupmenu_energy_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_energy contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_energy


% --- Executes during object creation, after setting all properties.
function popupmenu_energy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in radiobutton_init_equalarea.
function radiobutton_init_equalarea_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_init_equalarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_init_type, 'string', 'type: Equal Area');
guidata(gcbf,handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton_init_equalarea


% --- Executes on button press in radiobutton_init_random.
function radiobutton_init_random_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_init_random (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_init_type, 'string', 'type: random');
guidata(gcbf,handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton_init_random


% --- Executes on button press in radiobutton_init_workspace.
function radiobutton_init_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_init_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text_init_type, 'string', 'type: from workspace');
guidata(gcbf,handles);
% Call corresponding pushbutton callback function
success_flag = pushbutton_load_workspace_Callback(handles.pushbutton_load_workspace, eventdata, handles);
% If successful, everything is fine. If not, then print warning, because
% behavior is undefined.
if success_flag==0
    msgbox('WARNING: undefined behavior because loading process has been canceled! Make sure initial points are loaded properly.')
end

% Hint: get(hObject,'Value') returns toggle state of radiobutton_init_workspace


% --- Executes on button press in radiobutton_init_matfile.
function radiobutton_init_matfile_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_init_matfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text_init_type, 'string', 'type: from .mat file');
guidata(gcbf,handles);
% Call corresponding pushbutton callback function
success_flag = pushbutton_load_mat_file_Callback(handles.pushbutton_load_mat_file, eventdata, handles);
% If successful, everything is fine. If not, then print warning, because
% behavior is undefined.
if success_flag==0
    msgbox('WARNING: undefined behavior because loading process has been canceled! Make sure initial points are loaded properly.')
end

% Hint: get(hObject,'Value') returns toggle state of radiobutton_init_matfile


% --- Executes on button press in radiobutton_init_txtfile.
function radiobutton_init_txtfile_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_init_txtfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text_init_type, 'string', 'type: from .txt file');
guidata(gcbf,handles);
% Call corresponding pushbutton callback function
success_flag = pushbutton_load_txt_file_Callback(handles.pushbutton_load_txt_file, eventdata, handles);
% If successful, everything is fine. If not, then print warning, because
% behavior is undefined.
if success_flag==0
    msgbox('WARNING: undefined behavior because loading process has been canceled! Make sure initial points are loaded properly.')
end

% Hint: get(hObject,'Value') returns toggle state of radiobutton_init_txtfile


% --- Executes on button press in pushbutton_load_workspace.
function [success_flag] = pushbutton_load_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Xstart, Xstart_name] = uigetvar('double');
if ~isempty(Xstart) && ~isempty(Xstart_name)
    % Set compatible parameters
    handles.parameters.Xstart = Xstart;
    handles.parameters.Xstart_name = Xstart_name;
    D = size(Xstart,1);
    N = size(Xstart,2);
    handles.parameters.D = D;
    handles.parameters.N = N;
    % Set compatible GUI object properties
    set(handles.edit_D,'string',num2str(D));
    set(handles.edit_N,'string',num2str(N));
    handles.uibuttongroup_initialpoints.SelectedObject = handles.radiobutton_init_workspace;
    set(handles.text_init_D,'string',['D = ', num2str(D)]);
    set(handles.text_init_N,'string',['N = ', num2str(N)]);

    % Output and return value
    disp(['Loaded workspace variable ', Xstart_name, ' as Xstart.']);
    success_flag = 1;

    % Update handles structure
    guidata(gcbf, handles);
else
%     msgbox('CANCEL loading initial points from workspace. If algorithm is started now, the previously loaded Xstart will be used. If none was loaded before, it will not start.')
    success_flag = 0;
end


% --- Executes on button press in pushbutton_load_mat_file.
function [success_flag] = pushbutton_load_mat_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_mat_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[matfile_name, matfile_path] = uigetfile('*.mat','Pick a .mat file');
if matfile_name == 0
    success_flag = 0;
    return
else
    matfile_totalname = [matfile_path, matfile_name];
    load(matfile_totalname);
end
if ~exist('X','var')
    msgbox(['WARNING: no variable X contained in file ', matfile_totalname, ...
        '. Undefined behavior! Make sure initial points are loaded properly.']);
    success_flag = 0;
else
    % Set compatible parameters
    handles.parameters.Xstart = X;
    handles.parameters.Xstart_name = matfile_totalname;
    D = size(X,1);
    N = size(X,2);
    handles.parameters.D = D;
    handles.parameters.N = N;
    % Set compatible GUI object properties
    set(handles.edit_D,'string',num2str(D));
    set(handles.edit_N,'string',num2str(N));
    handles.uibuttongroup_initialpoints.SelectedObject = handles.radiobutton_init_matfile;
    set(handles.text_init_D,'string',['D = ', num2str(D)]);
    set(handles.text_init_N,'string',['N = ', num2str(N)]);

    success_flag = 1;
    disp(['Loaded .mat file ', matfile_totalname, ' containing X = Xstart.']);

    % Update handles structure
    guidata(gcbf, handles);
end


% --- Executes on button press in pushbutton_load_txt_file.
function [success_flag] = pushbutton_load_txt_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_txt_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[txtfile_name, txtfile_path] = uigetfile('*.txt','Pick a .txt file');
if txtfile_name == 0
    success_flag = 0;
    return
end
txtfile_totalname = [txtfile_path, txtfile_name];
Xstart = dlmread(txtfile_totalname);
Xstart = Xstart';  % ATTENTION: transposition!
if isempty(Xstart)
    msgbox(['WARNING: failed reading textfile ', txtfile_totalname]);
    success_flag = 0;
else
    % Set compatible parameters
    handles.parameters.Xstart = Xstart;
    handles.parameters.Xstart_name = txtfile_totalname;
    D = size(Xstart,1);
    N = size(Xstart,2);
    handles.parameters.D = D;
    handles.parameters.N = N;
    % Set compatible GUI object properties
    set(handles.edit_D,'string',num2str(D));
    set(handles.edit_N,'string',num2str(N));
    handles.uibuttongroup_initialpoints.SelectedObject = handles.radiobutton_init_txtfile;
    set(handles.text_init_D,'string',['D = ', num2str(D)]);
    set(handles.text_init_N,'string',['N = ', num2str(N)]);

    success_flag = 1;
    disp(['Loaded .txt file ', txtfile_totalname, ' containing X = Xstart.']);

    % Update handles structure
    guidata(gcbf, handles);
end


% --- Executes on button press in pushbutton_gamma_POU.
function pushbutton_gamma_POU_Callback(hObject, eventdata, handles)
    % Get contents of edit text fields
    gamma_min = str2double(get(handles.edit_gamma_min,'string'));
    gamma_max = str2double(get(handles.edit_gamma_max,'string'));
    gamma_steps = CheckInteger(str2double(get(handles.edit_gamma_steps,'string')),2);
    gamma_Nvali = CheckInteger(str2double(get(handles.edit_gamma_Nvali,'string')),1);
    % Basic consistency checks
    if gamma_min < 0.1
        msgbox('gamma_min must be at least 0.1')
    elseif isinf(gamma_min)
        msgbox('gamma_min must be < infinity')
    elseif gamma_max < 0.1
        msgbox('gamma_max must be at least 0.1')
    elseif isinf(gamma_max)
        msgbox('gamma_max must be < infinity')
    else
        % Begin the optimization
        disp('Begin kernel parameter optimization to fit constant function...')
        tic
        [ gamma_best, gamma_RMSE ] = SearchGamma( handles.results.X, gamma_Nvali, gamma_min, gamma_max, gamma_steps );
        toc

        % Display results
        set(handles.text_gamma_POU,'string',['result: gamma =  ', num2str(gamma_best)]);
        set(handles.text_gamma_POU_RMSE,'string',['RMSE = ', num2str(gamma_RMSE)]);

        % Store results
        handles.results.gamma_POU_value = gamma_best;
        handles.results.gamma_POU_RMSE = gamma_RMSE;

        % Update handles structure
        guidata(gcbf, handles);
    end




function edit_gamma_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gamma_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gamma_min as text
%        str2double(get(hObject,'String')) returns contents of edit_gamma_min as a double


% --- Executes during object creation, after setting all properties.
function edit_gamma_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gamma_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gamma_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gamma_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gamma_max as text
%        str2double(get(hObject,'String')) returns contents of edit_gamma_max as a double


% --- Executes during object creation, after setting all properties.
function edit_gamma_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gamma_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gamma_steps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gamma_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gamma_steps as text
%        str2double(get(hObject,'String')) returns contents of edit_gamma_steps as a double


% --- Executes during object creation, after setting all properties.
function edit_gamma_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gamma_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gamma_Nvali_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gamma_Nvali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gamma_Nvali as text
%        str2double(get(hObject,'String')) returns contents of edit_gamma_Nvali as a double


% --- Executes during object creation, after setting all properties.
function edit_gamma_Nvali_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gamma_Nvali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gradient_cycles_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gradient_cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gradient_cycles as text
%        str2double(get(hObject,'String')) returns contents of edit_gradient_cycles as a double


% --- Executes during object creation, after setting all properties.
function edit_gradient_cycles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gradient_cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gradient_iterations_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gradient_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gradient_iterations as text
%        str2double(get(hObject,'String')) returns contents of edit_gradient_iterations as a double


% --- Executes during object creation, after setting all properties.
function edit_gradient_iterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gradient_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_energy_cycles_Callback(hObject, eventdata, handles)
% hObject    handle to edit_energy_cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_energy_cycles as text
%        str2double(get(hObject,'String')) returns contents of edit_energy_cycles as a double


% --- Executes during object creation, after setting all properties.
function edit_energy_cycles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_energy_cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_energy_iterations_Callback(hObject, eventdata, handles)
% hObject    handle to edit_energy_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_energy_iterations as text
%        str2double(get(hObject,'String')) returns contents of edit_energy_iterations as a double


% --- Executes during object creation, after setting all properties.
function edit_energy_iterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_energy_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_export_workspace.
function pushbutton_export_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_export_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prefix            = get(handles.edit_export,'string');
prefix_properties = ( size(prefix)==[1,0] );
if strcmp(prefix, '') || strcmp(prefix, '.') || min(prefix_properties) == 1
    msgbox('ERROR: prefix must be chosen first')
else
    prefix_parameters = [prefix,'_parameters'];
    prefix_results    = [prefix,'_results'];
    abort_flag = 0;
    if ismember(prefix_parameters,evalin('base','who'))
        msgbox(['Struct ', prefix_parameters, ' already exists in base workspace. Not storing anything!'])
        abort_flag = 1;
    end
    if ismember(prefix_results,evalin('base','who'))
        msgbox(['Struct ', prefix_results, ' already exists in base workspace. Not storing anything!'])
        abort_flag = 1;
    end
    if abort_flag == 0
        assignin('base', prefix_parameters, handles.parameters);
        assignin('base', prefix_results, handles.results);
        disp(['Exported results to base workspace structures ', ...
            prefix_parameters, ' and ', prefix_results])
    end
end



% --- Executes on button press in pushbutton_export_files.
function pushbutton_export_files_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_export_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prefix            = get(handles.edit_export,'string');
path_export       = get(handles.text_export_path,'string');

prefix_properties = ( size(prefix)==[1,0] );
if strcmp(prefix, '') || strcmp(prefix, '.') || min(prefix_properties) == 1
    msgbox('ERROR: prefix must be chosen first')
else
    filename_mat = [path_export, '/', prefix, '.mat'];
    filename_txt = [path_export, '/', prefix, '.txt'];
    abort_flag = 0;
    if exist(filename_mat,'file') ~= 0
        msgbox(['File ', filename_mat, ' already exists. Not storing anything!'])
        abort_flag = 1;
    end
    if exist(filename_txt,'file') ~= 0
        msgbox(['File ', filename_txt, ' already exists. Not storing anything!'])
        abort_flag = 1;
    end
    if abort_flag == 0
        parameters = handles.parameters;
        results = handles.results;
        save(filename_mat,'parameters','results');
        dlmwrite(filename_txt,handles.results.X','delimiter','\t','precision','%30.23e'); % remember: transpose
    end
end

% --- Executes on button press in pushbutton_plotpoints3d.
function pushbutton_plotpoints3d_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plotpoints3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.parameters.D ~= 3
    msgbox(['Plot only possible for D = 3, but currently D = ', ...
        num2str(handles.parameters.D), '.']);
else
    PlotPointsOnSphere(handles.results.X, handles.results.largestgap, handles.parameters.sym_flag);
end



function edit_export_Callback(hObject, eventdata, handles)
% hObject    handle to edit_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_export as text
%        str2double(get(hObject,'String')) returns contents of edit_export as a double


% --- Executes during object creation, after setting all properties.
function edit_export_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_export_path.
function pushbutton_export_path_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_export_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname = uigetdir('.','Pick directory for export'); % returns 0 upon user caused cancel
if pathname ~= 0
    set(handles.text_export_path,'string',pathname);
else
    set(handles.text_export_path,'string','.');
end


% auxiliary function: set 'enable' status of all pushbuttons
function SetEnableAllPushbuttons(state, handles, current_fig)
if strcmp(state,'off')
    set(handles.pushbutton_generate_points,'enable','off');
    set(handles.pushbutton_plotEDF,'enable','off');
    set(handles.pushbutton_load_workspace,'enable','off');
    set(handles.pushbutton_load_mat_file,'enable','off');
    set(handles.pushbutton_load_txt_file,'enable','off');
    set(handles.pushbutton_export_workspace,'enable','off');
    set(handles.pushbutton_export_path,'enable','off');
    set(handles.pushbutton_export_files,'enable','off');
    set(handles.pushbutton_plotpoints3d,'enable','off');
    set(handles.pushbutton_gamma_POU,'enable','off');
    guidata(current_fig, handles);
elseif strcmp(state,'on')
    set(handles.pushbutton_generate_points,'enable','on');
    set(handles.pushbutton_plotEDF,'enable','on');
    set(handles.pushbutton_load_workspace,'enable','on');
    set(handles.pushbutton_load_mat_file,'enable','on');
    set(handles.pushbutton_load_txt_file,'enable','on');
    set(handles.pushbutton_export_workspace,'enable','on');
    set(handles.pushbutton_export_path,'enable','on');
    set(handles.pushbutton_export_files,'enable','on');
    set(handles.pushbutton_plotpoints3d,'enable','on');
    set(handles.pushbutton_gamma_POU,'enable','on');
    guidata(current_fig, handles);
else
    error('SetEnableAllPushbuttons needs state argument either ''on'' or ''off''');
end


% --- Executes on button press in pushbutton_plotEDF.
function pushbutton_plotEDF_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plotEDF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Plot empirical distribution functions
if( isempty(handles.results.X) )
    return;
end
disp('Plotting empirical distribution functions')
tic
cla(handles.axes1)
PlotEDF(handles.results.X, handles.parameters.sym_flag, handles.axes1);
toc


function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
handles.axes1 = hObject;
guidata(gcbf, handles);
axes1_OpeningFcn(hObject, eventdata, hObject)

function axes1_OpeningFcn(hObject, eventdata, ax, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
PlotEDF([],0,ax)
