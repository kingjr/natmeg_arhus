classdef XY3DDensityPlotData < handle
% Private class. It calculates density values (from supplied user data) on 
% which the surface plot is based.
%
% Calcuation of raw densities and smoothing is based on the Matlab File
% Exhange project "smoothist2D" (see acknowledgements in XY3DDensityPlotView.m).
    
    properties
        % user data
        x = [];
        y =[];
        
        % calculated data
        ctrs1 = [];  % grid centers
        ctrs2 = [];
        dpos = [];  % raw densities
        dneg = [];
        dposs = []; % smoothed densities
        dnegs = [];
        
        % outlier removal
        outlierSDThresh = 3;  % default may be modified
    end
    
    methods
       
        function obj = XY3DDensityPlotData(x,y, removeOutliers)
            % Constructor in which user data is assigned
            obj.x = x;
            obj.y = y;
            if nargin <=2 || removeOutliers
                obj.removeOutliers();
            end
        end
    
        function computeRawDensityGrid(obj, nbins)  
            % Compute raw densities from grid.
            
            % define grid over extent of data
            minx = min(obj.x,[],1);            
            maxx = max(obj.x,[],1);
            edges1 = linspace(minx(1), maxx(1), nbins(1)+1);
            obj.ctrs1 = edges1(1:end-1) + .5*diff(edges1);
            edges1 = [-Inf edges1(2:end-1) Inf];
            edges2 = linspace(minx(2), maxx(2), nbins(2)+1);
            obj.ctrs2 = edges2(1:end-1) + .5*diff(edges2);
            edges2 = [-Inf edges2(2:end-1) Inf];
            
            % raw pos densities
            pos = find( obj.y== 1 );
            xpos = obj.x(pos,[1 2]);
            [n,p] = size(xpos);
            bin = zeros(n,2);
            [na,bin(:,2)] = histc(xpos(:,1),edges1);
            [na,bin(:,1)] = histc(xpos(:,2),edges2);
            obj.dpos = accumarray(bin,1,nbins([2 1])) ./ length(obj.x);
            obj.dpos = accumarray(bin,1,nbins([1 2])) ./ length(obj.x);
            
            % raw neg densities
            neg = find(obj.y==-1);
            xneg = obj.x(neg,[1 2]);
            [n,p] = size(xneg);
            bin = zeros(n,2);
            [n,bin(:,2)] = histc(xneg(:,1),edges1);
            [n,bin(:,1)] = histc(xneg(:,2),edges2);            
            obj.dneg = accumarray(bin,1,nbins([1 2])) ./ length(obj.x);
        end        
        
        function smooth(obj, method, param)
            % Apply smoothing to the raw densities, using specified method
            % and associated smoothing parameter.  Methods are Eiler and
            % Filter.  The Eiler parameter is lambda. The Filter parameter
            % is the kernel matrix size.
            
            if strcmp(method, 'Eiler')
                temp = XY3DDensityPlotData.smooth1D(obj.dpos,param);
                obj.dposs = XY3DDensityPlotData.smooth1D(temp',param)';                               
                temp = XY3DDensityPlotData.smooth1D(obj.dneg,param);
                obj.dnegs = XY3DDensityPlotData.smooth1D(temp',param)';                
            else
                obj.dposs = XY3DDensityPlotData.filter2D(obj.dpos,param);
                obj.dnegs = XY3DDensityPlotData.filter2D(obj.dneg,param);
            end
        end
        
    end
    
    methods(Access='private')
        
        function removeOutliers(obj)
            % Removes outliers in x based on Mahalanobis distance.
            
            mscores = mahal(obj.x, obj.x);    
            m = mean(mscores);
            thresh = m + (obj.outlierSDThresh* std(mscores));             
            obj.x = obj.x(mscores < thresh,:);
            obj.y = obj.y(mscores < thresh,:);
        end

    end
    
    methods(Static)
        
         function Z = smooth1D(Y,lambda)
            % Eiler smoothing
            [m,n] = size(Y);
            E = eye(m);
            D1 = diff(E,1);
            D2 = diff(D1,1);
            P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
            Z = (E + P) \ Y;
        end
        
        function Z = filter2D(Y,bw)
            % Filter smoothing
            z = -1:(1/bw):1;
            k = .75 * (1 - z.^2); % epanechnikov-like weights
            k = k ./ sum(k);
            Z = filter2(k'*k,Y);
            assignin('base', 'k2', k'*k);
        end

    end
    
end