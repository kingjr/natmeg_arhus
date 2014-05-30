function angle = embed_pie(data,level,angle,total,colors)
% angle = embed_pie(data,level,angle,total,colors)
%data = {{{10, 20, 30},{10, 20, 30}},{{10, 20, 30},{10, 20, 30}},{{10, 20, 30},{10, 20, 30,}}}
x = 0;
y = 0;
r = 1;
border = .1;
border2 = .01;
linewdith = 1;

if nargin == 1;
    levels=0;new_level = data;
    while not(isempty(new_level));
        try new_level = new_level{1}; levels = levels+1;
        catch e,break;end
    end
    level = 0;
    angle = 0;
    total = sum(embed_pie_sum(data,1));
    %colors = {{'g', [1 .8 0], [1 0 0]} {'g', [1 0 0]} {'c', [.5 .5 .5], 'k'}};
    colors = {{[.8 .8 1], [.6 .6 .8], [.4 .4 .6]} {'g',[1 0 0]} {'c', [.5 .5 .5], 'k'}};
end


if iscell(data(1))
    gp_angle = angle;
    
    for gp = 1:length(data)
        angle = embed_pie(data{gp},level+1,angle,total,colors);                     % launch draw bottom
    end
%     pause;
if length(data) > 1
    for gp = 1:length(data)                                                 % draw upper group
        cat_nbs     = embed_pie_sum(data{gp},1);                            
        new_angle   = gp_angle+(sum(cat_nbs)/total).* 360;
        g           = arcpatch(x,y,r.*(level+1),[gp_angle new_angle]);
        [x1 y1]       = pol2cart(deg2rad(gp_angle),3);
        [x2 y2]       = pol2cart(deg2rad(new_angle),3);
        hold on;plot([0 x1;0 x2]', [0 y1;0 y2]', 'color', [1 1 1].* level .* .5);hold off;
        subtotal = sum(embed_pie_sum(data,1));
        set(g,'facecolor',  colors{level+1}{gp}, 'edgealpha', 1, 'edgecolor', [1 1 1], 'LineWidth', 2+linewdith/(subtotal/total));
        gp_angle       = new_angle;        
    end
end
   
       
else                                                                        % draw bottom
    for gp = 1:length(data)
        new_angle   = angle+(data(gp)/total).* 360;
%         g           = arcpatch(x,y,r.*(level+1),[angle new_angle]);
%         set(g,'facecolor', colors{level}{gp}, 'edgealpha', 1, 'edgecolor', [1 1 1],'LineWidth', 1)%, linewdith/(sum(data)/total));
        angle       = new_angle;
    end
end

