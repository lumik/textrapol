function varargout = textrapol(varargin)
% TEXTRAPOL MATLAB code for textrapol.fig
%      TEXTRAPOL, by itself, creates a new TEXTRAPOL or raises the existing
%      singleton*.
%
%      H = TEXTRAPOL returns the handle to a new TEXTRAPOL or the handle to
%      the existing singleton*.
%
%      TEXTRAPOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEXTRAPOL.M with the given input arguments.
%
%      TEXTRAPOL('Property','Value',...) creates a new TEXTRAPOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before textrapol_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to textrapol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help textrapol

% Last Modified by GUIDE v2.5 07-Jan-2015 10:28:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @textrapol_OpeningFcn, ...
                   'gui_OutputFcn',  @textrapol_OutputFcn, ...
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


% --- Executes just before textrapol is made visible.
function textrapol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to textrapol (see VARARGIN)

% Choose default command line output for textrapol
handles.output = hObject;

% functional handles
h=textrapol_functions;
[handles.textrapol_load, handles.faktorka,...
    handles.func, handles.residuum, handles.strfromvec,...
    handles.lmmin]=h{:};

set(handles.main_figure,'toolbar','figure');

set(handles.select_spec_panel,'Title','Select spectrum','Units',...
    'normalized','Position',[.05,.88,.4,.1],'Visible','off');
set(handles.chosen_spec_text,'Style','pushbutton','Units','normalized',...
    'Position',[.47 .15 .5 .8],'String','Chosen spectrum:','Parent',...
    handles.select_spec_panel,'BackgroundColor','green','Enable',...
    'inactive');
set(handles.choose_spec_pushbutton,'Style','pushbutton','Units','normalized','Position',...
    [.02 .15 .13 .8],'String','#','Parent',handles.select_spec_panel,...
    'TooltipString','Will chose one spectrum','Visible','on',...
    'Enable','on','FontSize',10);
set(handles.down_pushbutton,'Style','pushbutton','Units','normalized','Position',...
    [.17 .15 .13 .8],'String','<','Parent',handles.select_spec_panel,...
    'TooltipString','Will chose one spectrum','FontSize',10);
set(handles.up_pushbutton,'Style','pushbutton','Units','normalized','Position',...
    [.32 .15 .13 .8],'String','>','Parent',handles.select_spec_panel,...
    'TooltipString','Will chose one spectrum','FontSize',10);
set(handles.invert_pushbutton,'String','invert','Units',...
    'normalized','Position',[.47,.89,.06,.06],'Visible','off');

set(handles.U_axes,'Units','normalized');
set(handles.V_axes,'Units','normalized');
set(handles.main_figure,'Resize','on');


handles.sumfunc = @(vector) sum(vector.^2);

handles.current_directory = pwd; % nastaveni aktualniho adresare
handles.filepath = handles.current_directory;
handles.filename = '';
  
handles.spectra_orig = zeros(2,2);
handles.N_spectra = 1;
handles.factor_num = 1:2;
handles.invert_K = 1;

handles.polydeg = 0;
handles.extrapolation_interval = [3,20];
% handles.extrapolation_time = handles.extrapolation_interval(1);
handles.extrapolation_time = 0.5;
handles.plotstep = 1;

handles.data_loaded = false;

handles.H = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes textrapol wait for user response (see UIRESUME)
% uiwait(handles.main_figure);


% --- Outputs from this function are returned to the command line.
function varargout = textrapol_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function load_menuitem_Callback(hObject, eventdata, handles)
% hObject    handle to load_menuitem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filter_spec={'*.*','All Files (*.*)';
 '*.mat','MAT-files (*.mat)';
 '*.txt;*.dat','text files (*.txt,*.dat)'}; 
[status,soubor,cesta,data_orig]=handles.textrapol_load(0,handles.current_directory,...
  filter_spec,'Load Data','Off');
if status==0 % Pokud dojde k chybe pri nacitani dat (stisknuti cancel),
  % zobrazi se chybova hlaska.   
  h_errordlg=errordlg('No data load or invalid format of data !',...
      'loading data');
  waitfor(h_errordlg);
else % V pripade nacteni dat se provadi vse dalsi v Callbacku (prepisou se
  % predchozi data a smazou predchozi fity).
  
  handles.data_loaded = true;
  handles.x_scale = data_orig(:,1);
  handles.spectra_orig = data_orig(:,2:end);
  handles.N_spectra = size(handles.spectra_orig,2);
  handles.invert_K = ones(handles.N_spectra,1);
  
  handles.chosen_subspectrum = 1;
  
  set(handles.down_pushbutton,'Enable','off');
  if(handles.N_spectra == 1)
      set(handles.up_pushbutton,'Enable','off');
  else
      set(handles.up_pushbutton,'Enable','on');
  end
  handles.filepath = cesta;
  handles.current_directory = cesta;
  handles.filename = soubor;
  
  guidata(hObject, handles)
  
  treat_data(hObject)
  
  plot_function(hObject);
  
  handles = guidata(hObject);
  
  set(handles.invert_pushbutton, 'Visible', 'on', 'Enable', 'on')
  set(handles.select_spec_panel,'Visible','on');
  set(handles.chosen_spec_text, 'String', ...
    sprintf('Chosen spectrum: %d', handles.chosen_subspectrum), ...
    'BackgroundColor', 'green');

  guidata(hObject,handles)
 
end

function treat_data(hObject)

handles = guidata(hObject);

[handles.U,handles.W,handles.V,handles.E]=handles.faktorka(handles.spectra_orig);
handles.x_V = (1:size(handles.V,1))';
  
  
  
%   par0 = zeros(7 + (length(handles.factor_num) - 1) * 4, 1);
%   par0(1) = -3e-3;
%   par0(2) = -5;
%   par0(3) = 0.5;
%   par0(4) = -3e-3;
%   par0(5) = -1.5;
%   par0(6) = -5e-5;
%   par0(7) = -.157;
%   par0(8) = 1e0;
%   par0(9) = -5e-1;
%   par0(10) = 0;
%   par0(11) = 0;
% 
%   [handles.par,sigma,iter] = handles.lmmin(@(x,p) handles.residuum(...
%           handles.func,x,handles.V,handles.W,handles.factor_num,p),...
%           [],handles.sumfunc,handles.x_V,par0,1e2);
  
handles.polypar = zeros(handles.polydeg + 1, length(handles.factor_num));
fitint = handles.extrapolation_interval(1):handles.extrapolation_interval(2);
handles.V_extrapol = zeros(1,length(handles.factor_num));
for ii = 1:length(handles.factor_num)
    handles.polypar(:,ii) = polyfit(handles.x_V(fitint),handles.V(fitint,handles.factor_num(ii)),handles.polydeg);
    handles.V_extrapol(ii) = polyval(handles.polypar(:,ii),handles.extrapolation_time);
end

handles.extrapol_spectrum = handles.U(:,handles.factor_num)*diag(handles.W(handles.factor_num))*handles.V_extrapol';

guidata(hObject, handles);

% --- Executes on button press in choose_spec_pushbutton.
function choose_spec_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to choose_spec_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
prompt = {sprintf('Select one subspectrum\n(integer from 1 to %d)',...
    handles.N_spectra)}; 
dlg_title_matrix = 'Selection of subspectrum';
num_lines = 1;
options.Resize='on';
odpoved_num_of_spec=inputdlg(prompt,dlg_title_matrix,num_lines,{'',''},...
    options);
%--------------------------------------------------------------------------
% Testovani, zda nebylo stisknuto cancel
%--------------------------------------------------------------------------
if ~isequal(odpoved_num_of_spec,{})
    [num_chos_spec,status]=str2num(odpoved_num_of_spec{1});
    if num_chos_spec < 1 || num_chos_spec > handles.N_spectra || ...
            status == 0 || ...
            ~isequal(floor(num_chos_spec), num_chos_spec) || ...
            size(num_chos_spec,2) ~= 1
        % Nepripustne hodnoty
        set(handles.chosen_spec_text,'String', 'Incorrect number', ...
        'BackgroundColor','red'); % Chybna hodnota a jeji cervena indikace
        %errordlg('You selected incorret number of spectra','Incorret input','on'); 
    else % Pripustne hodnoty mezi matice se ulozi do globalnich promennych
        if num_chos_spec > 1
            set(handles.down_pushbutton, 'Enable', 'on');
        else
           set(handles.down_pushbutton, 'Enable', 'off');
        end
        if num_chos_spec<handles.N_spectra
            set(handles.up_pushbutton, 'Enable', 'on');
        else
            set(handles.up_pushbutton, 'Enable', 'off');
        end
        set(handles.chosen_spec_text, 'String', ...
            sprintf('Chosen spectrum: %d', num_chos_spec), ...
            'BackgroundColor', 'green');
        handles.chosen_subspectrum = num_chos_spec;
        guidata(hObject, handles);
        plot_function(hObject);
    end
end 


% --- Executes on button press in down_pushbutton.
function down_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to down_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.chosen_subspectrum <= 2
    handles.chosen_subspectrum = 1;
    set(handles.down_pushbutton,'Enable','off');
else
    handles.chosen_subspectrum = handles.chosen_subspectrum - 1;
end
if handles.chosen_subspectrum < handles.N_spectra;
    set(handles.up_pushbutton,'Enable','on');
end

set(handles.chosen_spec_text, 'String', ...
    sprintf('Chosen spectrum: %d', handles.chosen_subspectrum), ...
    'BackgroundColor', 'green');
guidata(hObject,handles);
plot_function(hObject);

% --- Executes on button press in up_pushbutton.
function up_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to up_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.chosen_subspectrum >= handles.N_spectra - 1
    handles.chosen_subspectrum = handles.N_spectra;
    set(handles.up_pushbutton,'Enable','off');
else
    handles.chosen_subspectrum = handles.chosen_subspectrum + 1;
end
if handles.chosen_subspectrum > 1
    set(handles.down_pushbutton,'Enable','on');
end

set(handles.chosen_spec_text, 'String', ...
    sprintf('Chosen spectrum: %d', handles.chosen_subspectrum), ...
    'BackgroundColor', 'green');

guidata(hObject,handles);
plot_function(hObject);

function plot_function(hObject)

handles = guidata(hObject);

gray1 = [0 0 0] + .6;
gray2 = [0 0 0] + .9;

cla(handles.U_axes);
cla(handles.V_axes);  

ii = handles.chosen_subspectrum;

plot(handles.U_axes,handles.x_scale,handles.invert_K(ii) * handles.U(:,ii));
aLims = axis(handles.U_axes);
xmin = min(handles.x_scale);
xmax = max(handles.x_scale);
dx = (xmax - xmin) / 50;
axis(handles.U_axes,[xmin - dx, xmax + dx, aLims(3:4)])

if handles.extrapolation_interval(1) > 1
    int1 = 1:handles.extrapolation_interval(1);
    int2 = 1:handles.extrapolation_interval(1) - 1;
    plot(handles.V_axes,handles.x_V(int1),handles.invert_K(ii) * handles.V(int1,ii),'-','color',gray1)
    hold(handles.V_axes,'on');    
    plot(handles.V_axes,handles.x_V(int2),handles.invert_K(ii) * handles.V(int2,ii),'o','MarkerEdgeColor',gray1,'MarkerFaceColor',gray2)
end

if handles.extrapolation_interval(1) == 1
    hold(handles.V_axes,'on');    
end

if handles.extrapolation_interval(2) < handles.N_spectra
    int1 = handles.extrapolation_interval(2) : handles.N_spectra;
    int2 = handles.extrapolation_interval(2) + 1 : handles.N_spectra;
    plot(handles.V_axes,handles.x_V(int1),handles.invert_K(ii) * handles.V(int1,ii),'-','color',gray1)
    if handles.extrapolation_interval(1) > 1
        hold(handles.V_axes,'on');
    end
    plot(handles.V_axes,handles.x_V(int2),handles.invert_K(ii) * handles.V(int2,ii),'o','MarkerEdgeColor',gray1,'MarkerFaceColor',gray2)
end

int1 = handles.extrapolation_interval(1):handles.extrapolation_interval(2);
H1 = plot(handles.V_axes,handles.x_V(int1),handles.invert_K(ii) *...
    handles.V(int1,ii),'-ok','MarkerFaceColor','g');

if ismember(ii,handles.factor_num)
    jj = find(ii == handles.factor_num);
    hold(handles.V_axes,'all');
    x_V_fit = handles.extrapolation_interval(1):handles.plotstep:handles.extrapolation_interval(2);
    H2 = plot(handles.V_axes, x_V_fit,handles.invert_K(ii) * polyval(handles.polypar(:,jj),x_V_fit),'-');
    H3 = plot(handles.V_axes,...
        handles.extrapolation_time,handles.invert_K(ii) *...
        polyval(handles.polypar(:,jj),handles.extrapolation_time),...
        'xr','MarkerSize',14, 'LineWidth',2);
    legend_text = {'V','fit'};
    legend([H1,H2],legend_text,'Location','Best')
else
    legend_text = {'V'};
    legend(H1,legend_text,'Location','Best')
end
% plot(handles.orig_axes,handles.x_scale,handles.spectra_orig(:,ii)/sqrt(sum(handles.spectra_orig(handles.fitintind,ii).^2)))
% handles.par(1) = -3e-3;
% handles.par(2) = -10;
% handles.par(3) = 0.5;
% handles.par(4) = -3e-3;
% handles.par(5) = -0.5;
% handles.par(6) = -5e-5;
% handles.par(7) = -.157;
% handles.par(8) = 1e0;
% handles.par(9) = -5e-1;
% handles.par(10) = 0;
% handles.par(11) = 0;
% par = handles.par
% if ismember(ii, handles.factor_num)
%     legend_text{end+1} = 'fit';
%     pind = 1:7;
%     I = find(ii == handles.factor_num);
%     if I > 1
%         pind(1) = 8 + (I-2) * 4;
%         pind(4) = 9 + (I-2) * 4;
%         pind(6) = 10 + (I-2) * 4;
%         pind(7) = 11 + (I-2) * 4;
%     end
%     plot(handles.V_axes,handles.x_V,handles.invert_K(ii) * handles.func(handles.x_V,handles.par(pind)),'-');
% end

% legend(handles.orig_axes,{'orig spectrum', 'subtracted spectrum','cacodylate','water','glass'})
guidata(hObject,handles);


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in invert_pushbutton.
function invert_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to invert_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.invert_K(handles.chosen_subspectrum) = -handles.invert_K(handles.chosen_subspectrum);
% handles.U(:,handles.chosen_subspectrum) = -handles.U(:,handles.chosen_subspectrum);
% handles.V(:,handles.chosen_subspectrum) = -handles.V(:,handles.chosen_subspectrum);


guidata(hObject, handles)
plot_function(hObject);


% --------------------------------------------------------------------
function save_menuitem_Callback(hObject, eventdata, handles)
% hObject    handle to save_menuitem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
savedata = [handles.x_scale,handles.extrapol_spectrum];
I = find(handles.filename == '.', 1, 'last');
filename = [handles.filepath,handles.filename(1:I-1),'_e.txt'];
fprintf('Saving file:\n%s\n',filename);
save(filename,'savedata','-ascii');
fprintf('Done.\n');


% --------------------------------------------------------------------
function options_menu_Callback(hObject, eventdata, handles)
% hObject    handle to options_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function settings_menuitem_Callback(hObject, eventdata, handles)
% hObject    handle to settings_menuitem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.settings_figure = figure('Position',[560,528,700,200],...
    'Toolbar','none','Resize','off','Name','Settings',...
    'Color',get(handles.main_figure,'Color'));

extrapol_interval_edit_string = mat2str(handles.extrapolation_interval);
handles.settings_extrapol_panel = uipanel('Title','Extrapolation',...
    'Units','pixels','Position',[20,150,660,40],'Visible','on');
handles.settings_extrapol_interval_text = uicontrol('Style','text',...
    'Parent',handles.settings_extrapol_panel,...
    'Units','pixels','Position',[10,6,50,20],'String','interval:',...
    'HorizontalAlignment','left');
handles.settings_extrapol_interval_edit = uicontrol('Style','edit',...
    'Parent',handles.settings_extrapol_panel,...
    'Units','pixels','Position',[60,6,100,20],...
    'String',extrapol_interval_edit_string,'HorizontalAlignment','left',...
    'BackgroundColor',[1,1,1]);

handles.settings_extrapol_time_text = uicontrol('Style','text',...
    'Parent',handles.settings_extrapol_panel,...
    'Units','pixels','Position',[170,6,50,20],'String','time:',...
    'HorizontalAlignment','left');
handles.settings_extrapol_time_edit = uicontrol('Style','edit',...
    'Parent',handles.settings_extrapol_panel,...
    'Units','pixels','Position',[220,6,100,20],...
    'String',num2str(handles.extrapolation_time),'HorizontalAlignment','left',...
    'BackgroundColor',[1,1,1]);

handles.settings_polydeg_panel = uipanel('Title','Set degree of polynom',...
    'Units','pixels','Position',[20,100,660,40],'Visible','on');
handles.settings_polydeg_text = uicontrol('Style','text',...
    'Parent',handles.settings_polydeg_panel,...
    'Units','pixels','Position',[10,6,50,20],'String','Degree:',...
    'HorizontalAlignment','left');
handles.settings_polydeg_edit = uicontrol('Style','edit',...
    'Parent',handles.settings_polydeg_panel,...
    'Units','pixels','Position',[60,6,590,20],...
    'String',num2str(handles.polydeg),'HorizontalAlignment','left',...
    'BackgroundColor',[1,1,1]);


handles.settings_important_factors_panel = uipanel('Title','Important factors',...
    'Units','pixels','Position',[20,50,660,40],'Visible','on');
handles.settings_important_factors_text = uicontrol('Style','text',...
    'Parent',handles.settings_important_factors_panel,...
    'Units','pixels','Position',[10,6,50,20],'String','list:',...
    'HorizontalAlignment','left');
handles.settings_important_factors_edit = uicontrol('Style','edit',...
    'Parent',handles.settings_important_factors_panel,...
    'Units','pixels','Position',[60,6,590,20],...
    'String',handles.strfromvec(handles.factor_num),'HorizontalAlignment','left',...
    'BackgroundColor',[1,1,1]);

handles.settings_ok_pushbutton = uicontrol('Style','pushbutton',...
    'Units','normalized','Position',[.8,.05,.08,.1],'String','OK',...
    'Callback',{@settings_ok_pushbutton_Callback,handles});
handles.settings_cancel_pushbutton = uicontrol('Style','pushbutton',...
    'Units','normalized','Position',[.9,.05,.08,.1],'String','Cancel',...
    'Callback',{@settings_cancel_pushbutton_Callback,handles});
guidata(hObject, handles);

% --------------------------------------------------------------------
function settings_ok_pushbutton_Callback(hObject, eventdata, handles)

handles = guidata(handles.main_figure);

handles.extrapolation_interval = eval(get(handles.settings_extrapol_interval_edit,'String'));
handles.extrapolation_time = str2double(get(handles.settings_extrapol_time_edit,'String'));
handles.polydeg = str2double(get(handles.settings_polydeg_edit,'String'));
handles.factor_num = eval(get(handles.settings_important_factors_edit,'String'));
close(handles.settings_figure);

guidata(handles.main_figure, handles);

if handles.data_loaded
    treat_data(handles.main_figure);
    plot_function(handles.main_figure);
end


function settings_cancel_pushbutton_Callback(hObject, eventdata, handles)

handles = guidata(handles.main_figure);

close(handles.settings_figure);

guidata(handles.main_figure, handles);


% --------------------------------------------------------------------
function view_menu_Callback(hObject, eventdata, handles)
% hObject    handle to view_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function spectrum_menuitem_Callback(hObject, eventdata, handles)
% hObject    handle to spectrum_menuitem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.H = figure('Toolbar','figure');
plot(handles.x_scale,handles.extrapol_spectrum,'-');

guidata(hObject, handles);
