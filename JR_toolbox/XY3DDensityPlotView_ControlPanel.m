classdef XY3DDensityPlotView_ControlPanel < handle
% Private class.  Creates the control panel.  Dependents (XY3DDensityPlotView.m)
% are informed of changes in values through events SmoothingChanged and 
% GridChanged. Default panel settings, as indicated in the properties below, 
% may be modified.

    properties

        % interface components and associated variables
        panel;
       
        gridText;
        gridEditBox;
        gridValue = 500;   % default may be modified
        
        smoothingText;
        smoothingEditBox;
        smoothingValue = 20;   % default may be modified
        
        smoothingSlider;
        sliderMin = 1;   % defaults may be modified
        sliderMax = 100;
        sliderStep = [0.02 0.1];
        
        rightPanel; 
        radioGroup;
        filterRadio;
        eilerRadio;  
        method = 'Eiler';  % default may be modified
        closeButton;

        panelMargin = 7;
        textBoxHeight = 21;
    end
    
    events
        SmoothingChanged
        GridChanged
    end
    
    methods
        
        function obj=XY3DDensityPlotView_ControlPanel(parent)
            % Create the panel interface in the supplied parent container.
            
            % the panel
            obj.panel = uipanel('parent', parent, 'units', 'pixels',...
                'BorderType', 'beveledout');
       
            % components to remain right aligned, placed in a separate
            % child panel
            obj.rightPanel = uipanel('parent', obj.panel, 'units',...
                'pixels','BorderType', 'none');
            
            % radio buttons
            obj.radioGroup = uibuttongroup('parent', obj.rightPanel, 'BorderType', 'none');
            pos = [0 0 100 obj.textBoxHeight];
            obj.eilerRadio = uicontrol('Style','Radio','String','Eiler',...
                'pos',pos,'parent',obj.radioGroup,'HandleVisibility','off');            
            gap = 7;            
            ext = get(obj.eilerRadio, 'extent');                         
            radioWidthKludge=10;
            pos = [ ext(3)+gap+radioWidthKludge 0 100 obj.textBoxHeight];
            obj.filterRadio = uicontrol('Style','Radio','String','Filter',...
                'pos',pos,'parent',obj.radioGroup,'HandleVisibility','off');
            set(obj.radioGroup,'SelectionChangeFcn', @obj.radioSelectionCallback); 
            if strcmp(obj.method, 'Eiler')
                set(obj.radioGroup,'SelectedObject',obj.eilerRadio);
            else
                 set(obj.radioGroup,'SelectedObject',obj.filterRadio);
            end
            
            % close button             
            pos = [ 110 0 80 obj.textBoxHeight ];
            obj.closeButton = uicontrol('style', 'pushbutton', 'string', 'Close',...
                'callback', @obj.closeButtonCallback, 'position', pos, 'parent', obj.rightPanel);
            
            % position right aligned panel
            panelWidth = get(obj.panel, 'position');
            panelWidth = panelWidth(3);
            rightPanelWidth = pos(1)+pos(3);          
            pos = [ panelWidth - rightPanelWidth - obj.panelMargin, obj.panelMargin, rightPanelWidth, obj.textBoxHeight];            
            set(obj.rightPanel, 'position', pos);
                        
            % the rest of the components
              
            % grid text
            kludgeToLowerText=3;
            pos = [obj.panelMargin obj.panelMargin-kludgeToLowerText 100 obj.textBoxHeight];         
            obj.gridText = uicontrol('style', 'text', 'Parent', obj.panel,...
               'HorizontalAlignment', 'left',... 
               'string', 'Grid:', 'position', pos);
           
            % grid edit box
            gap = 3;
            ext = get(obj.gridText, 'extent');
            posleft = pos(1) + ext(3) + gap;
            pos = [posleft obj.panelMargin 50 obj.textBoxHeight];
            gridEditBox = uicontrol('style', 'edit', 'parent', obj.panel,...
               'HorizontalAlignment', 'left', 'value', 500,...
               'callback', @obj.gridEditBoxCallback,...
               'position', pos, 'string', obj.gridValue);
           
            % smoothing text
            gap = 20;
            posleft = pos(1) + pos(3) + gap;
            pos = [posleft obj.panelMargin-kludgeToLowerText 100 obj.textBoxHeight];
            obj.smoothingText = uicontrol('style', 'text', 'Parent', obj.panel,...
               'HorizontalAlignment', 'left',...             
               'string', 'Smoothing:', 'position', pos);
            
            % smoothing edit box
            gap = 3;
            ext = get(obj.smoothingText, 'extent');
            posleft = pos(1) + ext(3) + gap;
            pos = [posleft obj.panelMargin 50 obj.textBoxHeight];
            obj.smoothingEditBox = uicontrol('style', 'edit', 'parent', obj.panel,...
               'HorizontalAlignment', 'left', 'value', obj.smoothingValue,...
               'callback', @obj.smoothingEditBoxCallback,...
               'position', pos, 'string', obj.smoothingValue);
           
            % smoothing slider
            gap = 7;
            posleft = pos(1) + pos(3) + gap;           
            pos=get(obj.rightPanel, 'position');
            width = pos(1) - gap - posleft;
            pos = [posleft obj.panelMargin width obj.textBoxHeight];
            obj.smoothingSlider =  uicontrol('Style','slider', 'Parent',obj.panel,...        
               'min', obj.sliderMin, 'max', obj.sliderMax,...
               'sliderStep', obj.sliderStep, 'value', obj.smoothingValue,...
               'Callback', @obj.sliderCallback, 'position', pos);
             
        end
        
        function resizeCallback(obj, eventdata, handles)
            % Called by container on resize to resize panel.
            
            % get container and its 'dimensions'
            parent = get(obj.panel, 'parent');
            pos = get(parent, 'position');
            
            % set new control panel extent 
            set(obj.panel, 'position', [0 0 pos(3) pos(4)] );
            
            % fix right aligned panel
            panelWidth = get(obj.panel, 'position');
            panelWidth = panelWidth(3);
            rightPanelWidth = get(obj.rightPanel, 'position');
            rightPanelWidth = rightPanelWidth(3);            
            pos = [ panelWidth - rightPanelWidth - obj.panelMargin, obj.panelMargin, rightPanelWidth, obj.textBoxHeight];            
            set(obj.rightPanel, 'position', pos);
                        
            % fix slider
            posRightPanel=get(obj.rightPanel, 'position');
            posSlider = get(obj.smoothingSlider, 'position');
            gap = 7;
            width = posRightPanel(1) - gap - posSlider(1);
            posSlider(3) = max(width,1);
            set(obj.smoothingSlider, 'position', posSlider);
            
        end
    
        function applySettings(obj,gridRes, method, param)
            % Programmatically change settings.            
            
            obj.gridValue = gridRes;
            set(obj.gridEditBox, 'string', num2str(gridRes));            
            
            obj.method = method;
            if strcmp(obj.method, 'Eiler')
                set(obj.radioGroup,'SelectedObject',obj.eilerRadio);
            else
                 set(obj.radioGroup,'SelectedObject',obj.filterRadio);
            end
            
            obj.smoothingValue = param;
            set(obj.smoothingEditBox, 'string', num2str(param));                       
            
            notify(obj, 'GridChanged'); 
            % no smoothing changed notification to avoid multiple plot redraws
        end
        
    end
        
    methods(Access='private')
        
        function gridEditBoxCallback(obj, eventdata, handles)
            % Grid value changed
            obj.gridValue = str2num(get(eventdata, 'string'));
            notify(obj, 'GridChanged');
        end
        
        function radioSelectionCallback(obj, eventdata, handles)
            % Radio button selection changed
            uic=get(eventdata, 'SelectedObject');
            obj.method = get(uic, 'string');
            notify(obj, 'SmoothingChanged');
        end
       
        function sliderCallback(obj, eventdata, handles)
            % Slider changed and associated smoothing value.
            v =get(eventdata, 'value');
            obj.smoothingValue = v;
            set(obj.smoothingEditBox, 'string', num2str(v));
            notify(obj, 'SmoothingChanged');
        end
                 
        function smoothingEditBoxCallback(obj, eventdata, handles)
            % Smoothing value changed.
            v = str2num(get(eventdata, 'string'));        
            set(obj.smoothingSlider, 'value', v);
            obj.smoothingValue = v;
            notify(obj, 'SmoothingChanged');
        end

        function closeButtonCallback(obj,eventdata, handles)
            % Close button clicked.
            close;
        end
        
    end
    
end