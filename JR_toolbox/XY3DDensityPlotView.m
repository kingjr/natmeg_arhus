classdef XY3DDensityPlotView < handle
    
% XY3DDensityPlotView is used to visualize the density of points described 
% by two variables X and Y, as a surface plot.  This may be useful when
% points are densely packed and overlapping, as in this case a scatterplot 
% may be limited in revealing their actual distribution.
%
% The XY data must be labelled, each data point belonging to one of two 
% classes (e.g. positive and negative).  A surface is produced on the 
% density of points in each class, and the two surfaces are shown together 
% on the same axes.
%
% Interface controls are provided to interactively refine smoothing
% parameters.  Two smoothing methods are provided, Eiler and Filter, which
% give a slightly different effect in terms of the surface produced.
%
% This class works in conjuction with two private classes. These are 
% XY3DensityPlotData, which computes density information for plotting, and 
% XY3DDensityPlotView_ControlPanel, which displays the interface components 
% which appear in the bar along the bottom of the figure window through 
% which parameters can be interactively modified.
%
% Usage:
%
% XY3DDensityPlotView(x, y, removeOutliers, xLabels)
%
% x: a 2 column array of numbers containing data points 'XY' ('Y'
%     being distinct from the class labels y, beneath).
%
% y: a 1 column array of class labels for each data point.  Class
%     may be either 1 or -1.  The former may be regarded as 'positive',
%     and the corresponding density surface plotted is coloured red; the 
%     latter blue.
%
% removeOutliers: a boolean indicating whether or not outlying points
%     should be removed from x before plotting.  The effect of this, if
%     there are extreme outliers, is to shorten the range of the axes and 
%     focus the plot on the main mass of points, rather than extremities.  
%     This may reduce whitespace.  A point is deemed an outlier based on 
%     the Mahalanobis distance, if it is greater than n standard deviations 
%     from the mean distance (n is defined in property "outlierSDThresh" 
%     in XY3DDensityPlotData.m).
% 
% xLabels: a cell array containing two strings which are the names of the
%     two variables in x.  These are used to label the axes.
%
% * Note that all parameters are required. 
% * Grid and smooth parameter settings may be applied programmically using 
%   the function "applySettings".
% * It has been tested only on Matlab version 2010b.
% * Data points are assumed to exist for both classes.
%
% Example:
%
% This example uses the demo data provided in the source files.
%
% load('demo data.mat');
% app = XY3DDensityPlotView(data.x, data.y, true, {'Var1', 'Var2'});
% app.applySettings(500, 'Filter', 20);
%
% Acknowledgements:
%
% This work is heavily based on the smoothhist2D code of Peter Perkins of
% Mathworks from the Matlab File Exchange at:
% http://www.mathworks.com/matlabcentral/fileexchange/13352. A copyright
% notice related to the use of this material is provided in license.txt.
% 
% The source file usercolormap.m, also from the File Exchange, is included 
% with kind permission of its author, Yo Fukushima. It can be found at:
% http://www.mathworks.com/matlabcentral/fileexchange/7144-usercolormap
%
% Author: 
%
% Jaspar Cahill, jjcahill@gmail.com. 19 Jan 2011.
% 

    properties
       data; % density data calculated for plotting

       % interface
       fig;   
       topPanel;
       bottomPanel;
       controlPanel;
       ax;

       labels;  % for labelling axes
    end
    
    methods
        
        function obj=XY3DDensityPlotView(x,y,removeOutliers, labels)
            % The main function, with arguments as described in the opening
            % comments.
            
            obj.labels = labels;
            
            % set up the figure interface
            obj.fig=figure('name', 'XY 3D Density Plot', 'resizefcn', @obj.resize);            
            obj.topPanel=uipanel('parent', obj.fig, 'units', 'pixels', 'BorderType', 'none');
            obj.ax = axes('parent', obj.topPanel);
            obj.bottomPanel=uipanel('parent', obj.fig, 'units', 'pixels', 'BorderType', 'none');
            obj.controlPanel=XY3DDensityPlotView_ControlPanel(obj.bottomPanel);
            
            % listen to changes in the control bar
            addlistener(obj.controlPanel, 'GridChanged', @obj.gridChanged);
            addlistener(obj.controlPanel, 'SmoothingChanged', @obj.smoothingChanged);

            % compute density data
            obj.data=XY3DDensityPlotData(x,y, removeOutliers);
            obj.data.computeRawDensityGrid([obj.controlPanel.gridValue obj.controlPanel.gridValue]);
            obj.data.smooth(obj.controlPanel.method, obj.controlPanel.smoothingValue);
                       
            % display plot
            obj.displayPlot();
            
        end
         
        function applySettings(obj, gridRes, method, param)
            % Programmatically change settings.
            % E.g. app.applySettings(500, 'Filter', 30);

            obj.controlPanel.applySettings(gridRes, method, param);            
        end
        
    end
    
    methods(Access='private')
        
        function displayPlot(obj)
            % Displays plot based on computed data in "obj.data".
   
            % get smoothed density data
            data = obj.data;
            dposs = data.dposs;
            dnegs = data.dnegs;

            % erase very small densities
            cutOff = 0.00000001;
            dposs(dposs<cutOff)=NaN;
            dnegs(dnegs<cutOff)=NaN;
            set(obj.ax, 'nextplot', 'replace');

            % display surfaces
            sNeg = surf(obj.ax,data.ctrs1,data.ctrs2,dnegs,'edgealpha',0, 'facealpha', 0.5);
            set(obj.ax, 'nextplot', 'add');
            sPos=surf(obj.ax,data.ctrs1,data.ctrs2,dposs,'edgealpha',0 );
 
            % colour surfaces          
            blueMap = usercolormap([0 0 0.4], [0 0 1],1);
            redMap = usercolormap([0.4 0 0], [1 0 0],1);
            set(obj.fig, 'colormap', [blueMap;redMap]);
            cmin = min(dnegs(:));
            cmax = max(dnegs(:));            
            m=256; % colors in each colormap
            C1 = min(m,round((m-1)*(dnegs-cmin)/(cmax-cmin))+1); 
            cmin = min(dposs(:));
            cmax = max(dposs(:));
            C2 = min(m,round((m-1)*(dposs-cmin)/(cmax-cmin))+1); 
            C2 = 256 + C2;
            set(sNeg,'CData',C1);
            set(sPos,'CData',C2);
            
            % label plot
            xlabel(obj.ax, obj.labels{1}, 'fontweight', 'bold');
            ylabel(obj.ax, obj.labels{2}, 'fontweight', 'bold');
            zlabel(obj.ax, 'Density');
            title(obj.ax, ['XY Density plot of ' obj.labels{1} ' and ' obj.labels{2}], 'fontweight', 'bold','interpreter', 'none');
        end
        
        function gridChanged(obj, eventdata, handles)            
            % Response to change in grid value event from control panel.           
            v=obj.controlPanel.gridValue;
            obj.data.computeRawDensityGrid([v v]);
            obj.data.smooth(obj.controlPanel.method, obj.controlPanel.smoothingValue);
            obj.displayPlot();
        end
        
        function smoothingChanged(obj, eventdata, handles)
            % Response to change in smoothing value event from control panel.
            obj.data.smooth(obj.controlPanel.method, obj.controlPanel.smoothingValue);          
            obj.displayPlot();
        end
        
        function resize(obj, eventdata, handles)
            % Callback on figure resize.
            bottomPanelHeight =  obj.controlPanel.panelMargin*2 + obj.controlPanel.textBoxHeight;
            fpos = get(obj.fig, 'position');            
            set(obj.topPanel, 'position', [0 bottomPanelHeight fpos(3) fpos(4)-bottomPanelHeight]);
            set(obj.bottomPanel, 'position', [0 0 fpos(3) bottomPanelHeight]);
            obj.controlPanel.resizeCallback();
        end
        
    end
end