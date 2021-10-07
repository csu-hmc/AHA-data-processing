function dflow_inspector
% loads a D-Flow data file and plots the channels

	scrsz = get(0,'ScreenSize');
	width = scrsz(3);
	height = scrsz(4);
	% position it at top of screen, whole width
	f = figure('Position',[10 height-370 width-15 320]);

	hfile=[];
	hslidera=[];
	hsliderb=[];
	hzoom=[];
	hunzoom=[];
	hrefresh=[];
    hquit = [];
	initializewindow;
	% initialize some variables to make sure they have this scope
	f1 = 0;
	f2 = 0;
	nsamples = 0;
	channela = 0;
	channelb = 0;
	t = [];
	data = [];
	channelnames={};
	nchannels = 0;	
	pathname = pwd;
	
	return
 
	%  Callbacks for simple_gui. These callbacks automatically
	%  have access to component handles and initialized data 
	%  because they are nested at a lower level.
   
	% Push button callbacks. Each callback plots current_data in
	% the specified plot type.
 
	function filebutton_Callback(source,eventdata) 
		% let user select file and load data
		[filename,pathname] = uigetfile('*.txt','Select a D-Flow data file', pathname);
		filename = [pathname filename];
        d = importdata(filename);
        channelnames = d.colheaders;
        t = d.data(:,1)-d.data(1,1);  % start at 0
        data = d.data;
		nsamples = size(data,1);
		nchannels = size(data,2);

		set(hslidera,'SliderStep',[1/(nchannels-1) 20/(nchannels-1)]);
		set(hsliderb,'SliderStep',[1/(nchannels-1) 20/(nchannels-1)]);
		set(f,'Name',['D-Flow file ' filename]);

		f1 = 1;
		f2 = nsamples;
		channela = 1;
		channelb = 2;
		showdata;
	end
	
	function slidera_Callback(source,eventdata) 
		% Change the channel number
		if (nchannels == 0)
			return
		end
		channela = round((nchannels-1)*get(source,'Value')+1);
		showdata;
	end
 
	function sliderb_Callback(source,eventdata) 
		% Change the channel number
		if (nchannels == 0)
			return
		end
		channelb = round((nchannels-1)*get(source,'Value')+1);
		showdata;
	end
 
	function zoombutton_Callback(source,eventdata) 
		% Wait for user to click on plot, and then zoom while centered there
		if (nchannels == 0)
			return
		end
		pos = ginput(1);
		isample = max(find(pos(1)>t));
		f1 = round(isample-(f2-f1)/4);
		f2 = round(isample+(f2-f1)/4);
		showdata
	end
 
	function unzoombutton_Callback(source,eventdata) 
		% View full timeline
		if (nchannels == 0)
			return
		end
		f1 = 1;
		f2 = nsamples;
		showdata
	end
 
	function refreshbutton_Callback(source,eventdata) 
		% Redraw the screen
		initializewindow;
		showdata;
    end

	function quitbutton_Callback(source,eventdata) 
        close(gcf);
	end
	
	function initializewindow
		clf
		%  Construct the components.
		hfile = uicontrol('Style','pushbutton','String','FILE',...
			'Units','Pixels','Position',[width-100,280,70,20],...
			'Callback',{@filebutton_Callback});
		hzoom = uicontrol('Style','pushbutton','String','ZOOM',...
			'Units','Pixels','Position',[width-100,250,70,20],...
			'Callback',{@zoombutton_Callback});
		hunzoom = uicontrol('Style','pushbutton','String','UNZOOM',...
			'Units','Pixels','Position',[width-100,220,70,20],...
			'Callback',{@unzoombutton_Callback});
		hrefresh = uicontrol('Style','pushbutton','String','REFRESH',...
			'Units','Pixels','Position',[width-100,190,70,20],...
			'Callback',{@refreshbutton_Callback});
		hquit = uicontrol('Style','pushbutton','String','QUIT',...
			'Units','Pixels','Position',[width-100,160,70,20],...
			'Callback',{@quitbutton_Callback});
		hslidera = uicontrol('Style','slider','BackgroundColor','b',...
	               'Units','Pixels','Position',[width-175 20 20 280],...
	               'Callback',@slidera_Callback);
		hsliderb = uicontrol('Style','slider','BackgroundColor','r',...
	               'Units','Pixels','Position',[width-145 20 20 280],...
	               'Callback',@sliderb_Callback);
		ha = axes('Units','Pixels','Position',[50,40,width-250,250]); 
	   
		% Initialize the GUI.
		% Change units to normalized so components resize 
		% automatically.
		set([f,hfile,hslidera,hsliderb,hzoom,hunzoom,hrefresh],'Units','normalized');
	end
 
	function showdata
		% plots two channels of data
		if (nchannels == 0)
			return
		end
		[ax,h1,h2] = plotyy(t(f1:f2), data(f1:f2,channela),t(f1:f2), data(f1:f2,channelb), 'plot');
		namea = channelnames{channela};
		namea = strrep(namea,'_','\_');		% escape the underline to prevent Latex formatting
		nameb = channelnames{channelb};
		nameb = strrep(nameb,'_','\_');		% escape the underline to prevent Latex formatting
		legend([h1 h2],[num2str(channela) ': ' namea],[num2str(channelb) ': ' nameb]);
		xlabel('time (s)')
	end
 
end 
