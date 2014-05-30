function [slope,intercept,STAT]=myregr(x,y,verbose)
%MYREGR: Perform a least-squares linear regression.
%This function computes a least-square linear regression suppling several
%output information.
%
% Syntax: 	myregr(x,y)
%      
%     Inputs:
%           X - Array of the independent variable 
%           Y - Dependent variable. If Y is a matrix, the i-th Y row is a
%           repeated measure of i-th X point. The mean value will be used
%           verbose - Flag to display all information (default=1)
%     Outputs:
%           - Slope with standard error an 95% C.I.
%           - Intercept with standard error an 95% C.I.
%           - Pearson's Correlation coefficient with 95% C.I. and its
%             adjusted form (depending on the elements of X and Y arrays)
%           - Spearman's Correlation coefficient
%           - Regression Standard Error
%           - Total Variability
%           - Variability due to regression
%           - Residual Variability
%           - Student's t-Test on Slope (to check if slope=0)
%           - Student's t-Test on Intercept (to check if intercept=0)
%           - Power of the regression
%           - a plot with:
%                o Data points
%                o Least squares regression line
%                o Red dotted lines: 95% Confidence interval of regression
%                o Green dotted lines: 95% Confidence interval of new y 
%                                       evaluation using this regression.
%
%   [Slope]=myregr(...) returns a structure of slope containing value, standard
%   error, lower and upper bounds 95% C.I.
%
%   [Slope,Intercept]=myregr(...) returns a structure of slope and intercept 
%   containing value, standard error, lower and upper bounds 95% C.I.
%
% Example:
%       x = [1.0 2.3 3.1 4.8 5.6 6.3];
%       y = [2.6 2.8 3.1 4.7 4.1 5.3];
%
%   Calling on Matlab the function: 
%             myregr(x,y)
%
%   Answer is:
%
%                         Slope
% -----------------------------------------------------------
%      Value          S.E.                95% C.I.  
% -----------------------------------------------------------
%    0.50107        0.09667        0.23267        0.76947
% -----------------------------------------------------------
%  
%                       Intercept
% -----------------------------------------------------------
%      Value          S.E.                95% C.I.  
% -----------------------------------------------------------
%    1.83755        0.41390        0.68838        2.98673
% -----------------------------------------------------------
%  
%             Pearson's Correlation Coefficient
% -----------------------------------------------------------
%      Value               95% C.I.                  ADJ  
% -----------------------------------------------------------
%    0.93296        0.49988        0.99281        0.91620
% -----------------------------------------------------------
% Spearman's Correlation Coefficient: 0.9429
%
%                       Other Parameters
% ------------------------------------------------------------------------
%      R.S.E.                         Variability  
%      Value        Total            by regression           Residual  
% ------------------------------------------------------------------------
%    0.44358        6.07333        5.28627 (87.0407%)     0.78706 (12.9593%)
% ------------------------------------------------------------------------
%
% Student's t-test on slope=0
% ----------------------------------------------------------------
% t = 5.1832    Critical Value = 2.7764     p = 0.0066
% Test passed: slope ~= 0
% ----------------------------------------------------------------
%  
% Student's t-test on intercept=0
% ----------------------------------------------------------------
% t = 4.4396    Critical Value = 2.7764     p = 0.0113
% Test passed: intercept ~= 0
% ----------------------------------------------------------------
%  
% Power of regression
% ----------------------------------------------------------------
% alpha = 0.05  n = 6     Zrho = 1.6807  std.dev = 0.5774
% Power of regression: 0.6046
% ----------------------------------------------------------------
%
% ...and the plot, of course.
%
% SEE also myregrinv, myregrcomp
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2007) MyRegression: a simple function on LS linear
% regression with many informative outputs. 
% http://www.mathworks.com/matlabcentral/fileexchange/15473

%Input error handling
if ~isvector(x)
    error('X must be an array.');
end

if isvector(y)
    if any(size(x)~=size(y))
        error('X and Y arrays must have the same length.');
    end
else
    if any(length(x)~=size(y,1))
        error('the length of X and the rows of Y must be equal.');
    end
end

if nargin==2
    verbose=1;
else
    verbose=logical(verbose);
end

x=x(:);
if isvector(y)
    yt=y(:); %columns vectors
else
    yt=mean(y,2);
end

if ~issorted(x) %x must be monotonically crescent 
    [x idx]=sort(x);
    yt=yt(idx);
end
xtmp=[x ones(length(x),1)]; %input matrix for regress function
ytmp=yt;

%regression coefficients
[p,pINT,R,Rint] = regress(ytmp,xtmp);

%check the presence of outliers
outl=find(ismember(sign(Rint),[-1 1],'rows')==0);
if ~isempty(outl) 
    fprintf(['These points are outliers at 95%% fiducial level: ' repmat('%i ',1,length(outl)) '\n'],outl); 
    reply = input('Do you want to delete outliers? Y/N [Y]: ', 's');
    disp(' ')
    if isempty(reply) || upper(reply)=='Y'
        ytmp(outl)=[]; xtmp(outl,:)=[];
        [p,pINT,R] = regress(ytmp,xtmp);
    end
end

xtmp(:,2)=[]; %delete column 2
%save coefficients value
m(1)=p(1); q(1)=p(2);

n=length(xtmp); 
xm=mean(xtmp); xsd=std(xtmp);

%standard error of regression coefficients
%Student's critical value
if isvector(y)
    cv=tinv(0.975,n-2); 
else
    cv=tinv(0.975,sum(size(y))-3);
end
m(2)=(pINT(3)-p(1))/cv; %slope standard error
m=[m pINT(1,:)]; %add slope 95% C.I.
q(2)=(pINT(4)-p(2))/cv; %intercept standard error
q=[q pINT(2,:)]; %add intercept 95% C.I.

%Pearson's Correlation coefficient
[rp,pr,rlo,rup]=corrcoef(xtmp,ytmp);
r(1)=rp(2); r(2)=realsqrt((1-r(1)^2)/(n-2)); r(3)=rlo(2); r(4)=rup(2); 
%Adjusted Pearson's Correlation coefficient
r(5)=sign(r(1))*(abs(r(1))-((1-abs(r(1)))/(n-2)));

%Spearman's Correlation coefficient
[rx]=tiedrank(xtmp);
[ry]=tiedrank(ytmp);
d=rx-ry;
rs=1-(6*sum(d.^2)/(n^3-n));

%Total Variability
ym=polyval(p,xm);
vtot=sum((ytmp-ym).^2);

%Regression Variability
ystar=ytmp-R;
vreg=sum((ystar-ym).^2);

%Residual Variability
vres=sum(R.^2);

%regression standard error (RSE)
if isvector(y)
    RSE=realsqrt(vres/(n-2));
else
    if ~isempty(outl) && (isempty(reply) || upper(reply)=='Y')
        y2=y; y2(outl)=[];
        RSE=realsqrt((vres+sum(sum((y2-repmat(ytmp,1,size(y,2))).^2)))/(sum(size(y2))-3));
    else
        RSE=realsqrt((vres+sum(sum((y-repmat(yt,1,size(y,2))).^2)))/(sum(size(y))-3));
    end
end

%Confidence interval at 95% of regression
sy=RSE*realsqrt(1/n+(((xtmp-xm).^2)/((n-1)*xsd^2)));
cir=[ystar+cv*sy ystar-cv*sy];

%Confidence interval at 95% of a new observation (this is the confidence
%interval that should be used when you evaluate a new y with a new observed
%x)
sy2=realsqrt(sy.^2+RSE^2);
cir2=[ystar+cv*sy2 ystar-cv*sy2];

%display results
if verbose==1
    tr=repmat('-',1,80);
    disp('                        Slope')
    disp(tr)
    disp('     Value          S.E.                95% C.I.  ')
    disp(tr)
    fprintf('%10.5f     %10.5f     %10.5f     %10.5f\n',m)
    disp(tr)
    disp(' ')
    disp('                      Intercept')
    disp(tr)
    disp('     Value          S.E.                95% C.I.  ')
    disp(tr)
    fprintf('%10.5f     %10.5f     %10.5f     %10.5f\n',q)
    disp(tr)
    disp(' ')
    disp('            Pearson''s Correlation Coefficient'); 
    disp(tr)       
    disp('     Value          S.E.                95% C.I.                 ADJ')
    disp(tr)
    fprintf('%10.5f     %10.5f     %10.5f     %10.5f       %10.5f\n',r)
    disp(tr)
    fprintf('Spearman''s Correlation Coefficient: %0.4f\n',rs)   
    disp(' ')
    disp('                      Other Parameters')
    disp(tr)
    disp('     R.S.E.                         Variability  ')
    disp('     Value        Total            By regression           Residual  ')
    disp(tr)
    fprintf('%10.5f     %10.5f       %10.5f (%0.1f%%)      %10.5f (%0.1f%%)\n',RSE,vtot,vreg,vreg/vtot*100,vres,vres/vtot*100)
    disp(tr)
    disp(' ')
    disp('Press a key to continue'); pause; disp(' ')

    %test on slope
    t=abs(m(1)/m(2)); %Student's t
    disp('Student''s t-test on slope=0')
    disp(tr)
    fprintf('t = %0.4f    Critical Value = %0.4f     p = %0.4f\n',t,cv,pr(2))
    if t>cv
        disp('Test passed: slope ~= 0')
    else
        disp('Test not passed: slope = 0')
        m(1)=0;
    end
    try
        powerStudent(t,n-1,2,0.05)
    catch ME
        disp(ME)
         disp('I am trying to download the powerStudent function by Antonio Trujillo Ortiz from FEX')
         [F,Status]=urlwrite('http://www.mathworks.com/matlabcentral/fileexchange/2907-powerstudent?controller=file_infos&download=true','powerStudent.zip');
         if Status
             unzip(F)
             powerStudent(t,n-1,2,0.05)
         end
         clear F Status
    end
    disp(tr)
    disp(' ')
    %test on intercept
    t=abs(q(1)/q(2)); %Student's t
    p=(1-tcdf(t,n-2))*2; %p-value
    disp('Student''s t-test on intercept=0')
    disp(tr)
    fprintf('t = %0.4f    Critical Value = %0.4f     p = %0.4f\n',t,cv,p)
    if t>cv
        disp('Test passed: intercept ~= 0')
    else
        disp('Test not passed: intercept = 0')
        q(1)=0;
    end
    powerStudent(t,n-1,2,0.05)
    disp(tr)
    disp(' ')
    %Power of regression
    Zrho=0.5*reallog((1+abs(r(1)))/(1-abs(r(1)))); %normalization of Pearson's correlation coefficient
    sZ=realsqrt(1/(n-3)); %std.dev of Zrho
    pwr=1-tcdf(1.96-Zrho/sZ,n-2)*2; %power of regression
    disp('Power of regression')
    disp(tr)
    fprintf('alpha = 0.05  n = %d     Zrho = %0.4f  std.dev = %0.4f\n',n,Zrho,sZ)
    fprintf('Power of regression: %0.4f\n',pwr)
    disp(tr)
    disp(' ')
    disp('Press a key to continue'); pause

    %plot regression
    if isvector(y)
        plot(x,y,'b.',x,y,'bo',xtmp,ystar,xtmp,cir,'r:',xtmp,cir2,'g:');
    else
        hold on
        plot(x',y,'LineStyle','none','Marker','o','MarkerEdgeColor','b')
        plot(xtmp,ystar,'k',xtmp,cir,'r:',xtmp,cir2,'g:');
        hold off
    end
    axis tight
    axis square

    disp('Red dotted lines: 95% Confidence interval of regression')
    disp('Green dotted lines: 95% Confidence interval of new y evaluation using this regression')
end

switch nargout
    case 1
        slope.value=m(1);
        slope.se=m(2);
        slope.lv=m(3);
        slope.uv=m(4);
    case 2
        slope.value=m(1);
        slope.se=m(2);
        slope.lv=m(3);
        slope.uv=m(4);
        intercept.value=q(1);
        intercept.se=q(2);
        intercept.lv=q(3);
        intercept.uv=q(4);
    case 3
        slope.value=m(1);
        slope.se=m(2);
        slope.lv=m(3);
        slope.uv=m(4);
        intercept.value=q(1);
        intercept.se=q(2);
        intercept.lv=q(3);
        intercept.uv=q(4);
        STAT.rse=RSE;
        STAT.cv=cv;
        STAT.n=n;
        STAT.xm=mean(x);
        STAT.ym=ym;
        STAT.sse=sum((xtmp-xm).^2);
        STAT.r=r;
end