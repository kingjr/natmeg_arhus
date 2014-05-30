function myeb2(Y,mycolor,varargin)
% myeb2(Y,mycolor,varargin)



mycolor = mycolor ./ 4;

col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol=col+.8; ccol(ccol>1)=1;


if isempty(varargin);
	if length(size(Y))==2 
		m=median(Y);
		s=nanstd(Y) / size(Y,1);
		ind1=1:length(m);
		ind2=ind1(end:-1:1);
        a = [1.96 2.576];%quantile for alpha = 0.05, and alpha = 0.01
		hold on; g=fill([ind1 ind2],[m-(a(2) * s) m(ind2)+(a(2) * s(ind2))],mycolor .* 4);
        hold on; h=fill([ind1 ind2],[m-(a(1) * s) m(ind2)+(a(1) * s(ind2))],mycolor .* 3);        
		plot(ind1,m,'linewidth',2,'Color',mycolor .* 2);
        set(h,'edgecolor',mycolor .* 3, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);  
        set(g,'edgecolor',mycolor .* 4, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);         
		hold off
	elseif length(size(Y))>2 
		cla; hold on; 
		ind1=1:size(Y,2);
		ind2=ind1(end:-1:1);
		if size(Y,3)>8; col=jet(size(Y,3));ccol=col+.8; ccol(ccol>1)=1;end
		for k=1:size(Y,3)
			m=nanmean(Y(:,:,k));
			s=nanstd(Y(:,:,k));
			h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],ccol(k,:));
			set(h,'edgecolor',ccol(k,:))
		end
		for k=1:size(Y,3)
			m=nanmean(Y(:,:,k));
			s=nanstd(Y(:,:,k));
			plot(ind1,m,'linewidth',2,'color',col(k,:))
		end
		hold off 
	end

elseif length(varargin)==1;

	m=Y;
	s=varargin{1};
	if length(size(Y))>2; error;
	elseif min(size(Y))==1;
		if size(m,1)>1; m=m';s=s';end
		ind1=1:length(m);
		ind2=ind1(end:-1:1);
		hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],.6*ones(1,3));
		set(h,'edgecolor',.6*ones(1,3))
		plot(ind1,m,'linewidth',2)
		hold off
	else 
		ind1=(1:size(Y,1));
		ind2=ind1(end:-1:1);
		cla; hold on; 
		if size(Y,2)>8; col=jet(size(Y,2));ccol=col+.8; ccol(ccol>1)=1;end
		for k=1:size(Y,2)
			mm=m(:,k)';
			ss=s(:,k)';
			h=fill([ind1 ind2],[mm-ss mm(ind2)+ss(ind2)],ccol(k,:));
			set(h,'edgecolor',ccol(k,:))
		end
		for k=1:size(Y,2);
			mm=m(:,k)';
			ss=s(:,k)';
			plot(ind1,mm,'linewidth',2,'color',col(k,:))
		end
		hold off 
	end
end


