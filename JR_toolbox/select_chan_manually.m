load('/media/DONNEES/JR KING bureau/Pro/Toolbox/JR_toolbox/my_EGI_net.mat', 'layout');
h=figure(1);clf;set(gcf,'color','k');
layout.pos = round(layout.pos .* 100);
hold on;
for chan = 1:256
    chans(chan) = scatter(layout.pos(chan,1),layout.pos(chan,2),'w');
    set(chans(chan), 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w', 'SizeData', 175);
end
axis off; axis image;
select = [];
%46 chans:
pre = [8 11 21 30 31 38 48 52 57 67 69 75 78 83 90 91 99 103 107 123 126 133 141 147 154 158 160 174 180 184 191 200 202 204 215 216 219 222 226 233 238 240 241 243 251 252];


for s = 1:length(pre)
    set(chans(pre(s)),'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    select = unique([select, pre(s)]);
    chan_sym = intersect(...
        find(layout.pos(:,1)==-layout.pos(pre(s),1)), ...
        find(layout.pos(:,2)==layout.pos(pre(s),2)));
    set(chans(chan_sym),'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
    select = unique([select, chan_sym]);
end

while 1
    % -- click
    key=waitforbuttonpress;
    index = gco(h);
    if strcmp(get(index,'Type'), 'hggroup')
        chan = find(chans==index);
        if sum(get(index,'MarkerFaceColor'))==3 % if black (unselected)
            set(index,'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
            select = unique([select chan]);
        else
            set(index,'MarkerFaceColor', 'w','MarkerEdgeColor', 'w');
            select = setdiff(select, chan);
        end
        %-- find symetrical
        chan_sym = intersect(...
            find(layout.pos(:,1)==-layout.pos(chan,1)), ...
            find(layout.pos(:,2)==layout.pos(chan,2)));
        if ~isempty(chan_sym) && layout.pos(chan_sym,1) ~= 0
            index = chans(chan_sym);
            if sum(get(index,'MarkerFaceColor'))==3 % if black (unselected)
                set(index,'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
                select = unique([select chan_sym]);
            else
                set(index,'MarkerFaceColor', 'w','MarkerEdgeColor', 'w');
                select = setdiff(select, chan_sym);
            end
        end
    end
    clc;
    disp(select);
    disp(length(select));
end