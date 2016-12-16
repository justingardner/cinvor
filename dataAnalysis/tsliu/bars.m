function h=bars(f,s,xstr,ystr,col,LineWidth,w,we);
%
% bars(f,[s],[xstr],[ystr],[color],[linewidth],[barwidth],[errorbarwidth]);
% 
%  Draws a bar graph with heights determined by f
%  and error-bars determined by s.
%
%  if input matrices are 2d, the ROWS are plotted as
%  separate sub-bars in different colors.
%
%  xlabels set to 'xstr' (columns)
%  legend set to 'ystr' (rows)
%
%  color is a string of colors, one for each subbar
%  (default is 'brgmyc')
%
%  linewidth (default 2)
%
%  barwidth (default 0.65)
%  errorbarWidth (default 0.3)
%
%  Example:
%    f = rand(5,1)+1;
%    s = rand(5,1)/2;
%    xstr = str2mat('one','two','three','four','five');
%    ystr = '';
%    mybar(f,s,xstr,ystr,'g');
%  7/11/95  gmb Wrote it.
%  4/11/96  gmb Revised it to plot 2d data sets
%  11/19/97 gmb Converted it to Matlab version 5.0
%  12/30/97 djh Fixed bug in xlabels, for some reason it was
%               putting the 1st label under the 2nd bar, and so on.
%  1/24/97  djh Added linewidth as optional input arg
%  3/15/01  dbr Made legend size allocation a bit more intelligent.
%  05/2004  arw Modified so that if you pass in a list of colors for a
%  single data vector, it applies those colors to the individual bars.
%           EXAMPLE:
%          a=rand(4,1);b=rand(4,1)/5; h=mybar(a,b,[],[],'rgby');
% OR:      mybar(a,b,[],[],[1 0 0;0 1 0;0 0 1;1 1 0]);

if ~exist('col')
    col='brgmyckbrgmyck';
end
if ~exist('LineWidth')
    LineWidth=1;
end
if ~exist('w')
    w=0.65;
end
if ~exist('we')
    we=0.3;
end

if nargin==1
    s=zeros(size(f));
end

if(size(f,1)==1)
    f=f';
    s=conj(s');
end

if (isreal(s))
    s = s+sqrt(-1)*s;
end

nsubbars=size(f,2);
we=we/nsubbars;

heldstate=ishold;
plot(0,0)
set(gca,'XLim',[1-w,size(f,1)+w]);

yhilim = max(max(f+s))*1.2;
if ~isnan(yhilim)
    %set(gca,'YLim',[0,yhilim]);
end

subx=linspace(-w/2,w/2,nsubbars+1);

hold off;

for subbar=1:nsubbars
    x=[];
    y=[];
    for i=1:size(f,1);
        tempx= [i+subx(subbar),i+subx(subbar),i+subx(subbar+1),i+subx(subbar+1)];
        tempy= [0,f(i,subbar),f(i,subbar),0];
        
        if ((length(col)==length(f)) & (size(f,2)==1))
            
            
            % disp('.');
            if isnumeric(col) & size(col,2)==3
                thisCol=col(i,:);
            else
                thisCol=col(i);
            end
            
            h((subbar-1)*size(f,1)+i)=fill(tempx',tempy',thisCol);
            hold on;
            
        else
            if (size(col,2)==3)
                thisCol=col(subbar,:);
            else
                thisCol=col(subbar);
            end
            
            h(subbar*size(f,1)+i)=fill(tempx',tempy',thisCol);
            hold on;
            
        end
        
        x=[x';tempx]';
        y=[y';tempy]';
    end
    %   if ((length(col)==size(f,1)) & (size(f,2)==1))
    %      % do somethign!
    %   end
    %  
    %   h(subbar)={fill(x,y,col(subbar))};
    %   
    a=line(get(gca,'XLim'),[0,0]);
    set(a,'Color','k');
    
    errx=mean([subx(subbar),subx(subbar+1)]);
    for i=1:size(f,1)
        line([i+subx(subbar),i+subx(subbar),i+subx(subbar+1),i+subx(subbar+1)], ...
            [0,f(i,subbar),f(i,subbar),0],'Color','k','LineWidth',LineWidth);
        g=line([i+errx,i+errx],[f(i,subbar)-imag(s(i,subbar)),f(i,subbar)+real(s(i,subbar))],...
            'Color','k','LineWidth',LineWidth);
        g=line([i+errx-we/2,i+errx+we/2],[f(i,subbar)-imag(s(i,subbar)),f(i,subbar)-imag(s(i,subbar))],...
            'Color','k','LineWidth',LineWidth);
        g=line([i+errx-we/2,i+errx+we/2],[f(i,subbar)+real(s(i,subbar)),f(i,subbar)+real(s(i,subbar))],...
            'Color','k','LineWidth',LineWidth);
    end
    %hold on
end 					%subbars

% Xlabels
set(gca,'xLimMode','manual');
set(gca,'XTick',[1:size(f,1)]);
if (exist('xstr'))
    if ~isempty(xstr)
        set(gca,'XTickLabel',xstr);
    end
end

% Legend
if (exist('ystr'))
    if ~isempty(ystr)
        %heighten the graph
        ylim=get(gca,'Ylim');
        nLines = size(ystr, 1);
        yMult = 2.5; %1 + nLines/4; % Each legend line uses about 1/8th of the window height
        set(gca,'YLim',[ylim(1),ylim(2)*yMult]);
        bs=3;
        
        posgca=get(gca,'Position');
        posgcf=get(gcf,'Position');	
        dx=diff(get(gca,'XLim'))/(posgca(3)*posgcf(3));
        dy=diff(get(gca,'YLim'))/(posgca(4)*posgcf(4));
        
        bw=bs*dx; 				%box width
        bh=bs*dy; 				%box height	
        
        yspace=2.5*bh; 			%vertical spacing	
        
        xlim=get(gca,'Xlim');
        ylim=get(gca,'Ylim');
        
%         xc=xlim(1)+50*dx;
        xc=xlim(2)-75*dx;
        yc=ylim(2)-30*dy;
        bx=[xc-bw,xc+bw,xc+bw,xc-bw,xc-bw];
        by=[yc+bh,yc+bh,yc-bh,yc-bh,yc+bh];
        
        textx=xc+bw+10*dx;
        texty=yc;
        
        %draw box and text
        for i=1:size(ystr,1);
            fill(bx,by,col(i));
            line(bx,by,'LineWidth',LineWidth,'Color','k');
            text(textx,texty,ystr(i,:));
            by=by-yspace;
            texty=texty-yspace;
            
        end
    end
end
set(gca,'XLim',[1-w,size(f,1)+w]);

if (heldstate==0)
    hold off
end


