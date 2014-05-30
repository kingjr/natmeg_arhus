function data = build_tree(data)
% data.value
% data.width
% data.color
% data.children
% (C) JeanRémi King, jeanremi.king+matlab@gmail.com
colors = colormap('jet');
branch_line = bezierInterp(0, 0, 1, 1);
startx = linspace(0,1,101)';
starty = zeros(101,1);
if ~isfield(data,'color'),data.color = [0 0 0]; end
data.handle = plot_branch(data,0,0,1,1,data.width,data.color,branch_line,1,0);
data = add_child(1,0,branch_line,colors,data,1);
return

function data = add_child(startx,starty,branch_line,colors,data,gener )
brothers = length(data.children);
plot_on = true;
gener = gener+1;% change generation
add_y = 0;
for brother = 1:brothers
    % define new branch
    eval(['data.' data.children{brother} '=plot_branch(data.' data.children{brother} ',startx,starty,brother,brothers,data.width,colors,branch_line,gener,add_y);']);
    %add_y = add_y+ eval(['data.' data.children{brother} '.width']);
    % if next node has children
    if ismember('children',...
            fieldnames(eval(['data.' data.children{brother}])))
        data_child = add_child(...
            startx+1,add_y.*0.1+starty+reduce(1,gener,brother,brothers),...
            branch_line,colors,...
            eval(['data.' data.children{brother}]),gener);
        eval(['data.' data.children{brother} '=data_child;']);
    end
    add_y = add_y+ eval(['data.' data.children{brother} '.width']);
end
return

function branch_line= reduce(branch_line,gener,brother,brothers)
% reduce branch according to its generation
branch_line = branch_line'/factorial(gener)*((brother-1)/(brothers-1)-.5);
if isnan(branch_line)
    branch_line = 0;
end
return

function data = plot_branch(data,startx,starty,brother,brothers,width,colors,branch_line,gener,add_y)
% calculate branch from data structure and generations infos and plot them
if width > 0
    smooth = 3;
    smooth_type = 'sigmoid';
    smooth_type = 'straight';
    base_width = .1;
    % define new branch
    data.branch         = [];
    %-- x
    data.branch(1,:)    = [...
        startx+linspace(0,1,length(branch_line)), ...
        startx+linspace(1,0,length(branch_line))];
    %-- y
    if strcmp(smooth_type,'straight')
        starty = starty + (add_y+data.width)*base_width;
    end
    if brothers == 1 % handle calculation error for no brothers
        data.branch(2,:) = [...
            starty+0*branch_line', ...
            starty+0*branch_line(end:-1:1)'];
    else
        data.branch(2,:) = [...
            starty+reduce(branch_line,gener,brother,brothers), ...
            starty+reduce(branch_line(end:-1:1),gener,brother,brothers)];
    end
    % evaluate color
    if ~isfield(data, 'color'), data.color = [0 0 0]; end
    if length(data.color) == 1
        data.color = colors(round((length(colors) - 1) * data.color)+1,:);
    end
    %data.handle = plot(data.branch(1,:),data.branch(2,:),...
    % 'branch_linewidth', base_width * width,...
    % 'color', color);

    switch smooth_type
        case 'sigmoid'
            sig =  [1 1-(1-data.width) ./(1 + exp(-(linspace(-1,1,length(data.branch)/2-2)*smooth))) data.width];
            top     = data.branch(2,1:floor(end/2))+sig*base_width*width;
            bottom  = data.branch(2,(floor(end/2)+1):end)-sig(end:-1:1)*base_width*width;
        case 'branch_linear'
            sig = linspace(1,data.width,length(data.branch)/2);
            top     = data.branch(2,1:floor(end/2))+sig*base_width*width;
            bottom  = data.branch(2,(floor(end/2)+1):end)-sig(end:-1:1)*base_width*width;        
        case 'straight'
            top     = data.branch(2,1:floor(end/2));
            bottom  = -data.width*base_width+data.branch(2,(floor(end/2)+1):end);
            
    end
    data.handle = fill(...
        data.branch(1,:),...
        [top bottom],...
        data.color,'EdgeColor','none');hold on;
%     set(data.handle,'FaceAlpha',0);
end
return
