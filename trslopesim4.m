%拟合+获取特征斜率
clear;clc;
LoopCircles=load('SJ.txt');

%归一化,翻转
Lmax=max(LoopCircles);
Lmin=min(LoopCircles);
LoopCircles(:,1)=-2*(LoopCircles(:,1)-(Lmax(1)+Lmin(1))/2)/(Lmax(1)-Lmin(1));
LoopCircles(:,2)=-2*(LoopCircles(:,2)-(Lmax(2)+Lmin(2))/2)/(Lmax(2)-Lmin(2));

%取半周
[LMaxValue,LineMax]=max(LoopCircles);         
[LMinValue,LineMin]=min(LoopCircles); 
if LineMin(1)<LineMax(1)
HalfL=LoopCircles(LineMin(1):LineMax(1),:);
else
HalfL=LoopCircles(LineMax(1):LineMin(1),:);
HalfL=flipud(HalfL); %有的输入顺序颠倒
end
%拟合
LineNum=size(HalfL,1);
% Lp = polyfit(HalfL(:,1),HalfL(:,2),3);%次数取决于曲线
% xmesh=linspace(HalfL(1,1),HalfL(LineNum,1),1000);

%线性插值
xmesh=linspace(HalfL(1,1),HalfL(LineNum,1),37);
FeatureP=ones(37,2);
FeatureP(:,1)=xmesh;
delta=(FeatureP(37,1)-FeatureP(1,1))/36;
FeatureP(1,:)=HalfL(1,:);
FeatureP(37,:)=HalfL(LineNum,:);
for k=2:36
[Fvalue,Findex]=min(abs(HalfL(:,1)-FeatureP(1,1)-(k-1)*delta));

FeatureP1(k,:)=HalfL(Findex,:);
    if (FeatureP1(k,1)-FeatureP(1,1)-(k-1)*delta)>0
        FeatureP2(k,:)=HalfL(Findex-1,:);
    else
        FeatureP2(k,:)=HalfL(Findex+1,:);
    end
    slope=(FeatureP1(k,2)-FeatureP2(k,2))/(FeatureP1(k,1)-FeatureP2(k,1));
    FeatureP(k,2)=slope*(FeatureP(1,1)+(k-1)*delta-FeatureP2(k,1))+FeatureP2(k,2);
end
%旋转平放
angle=-pi/4;
for k=1:37
    rotate_x(k)=FeatureP(k,1)*cos(angle)-FeatureP(k,2)*sin(angle);
    rotate_y(k)=FeatureP(k,1)*sin(angle)+FeatureP(k,2)*cos(angle);
end
 FeatureP=[rotate_x.' rotate_y.'];
 for k=1:36
Fslope(k)=(FeatureP(k+1,2)-FeatureP(k,2))/(FeatureP(k+1,1)-FeatureP(k,1));
end
%圆，弓，反S，Z的斜率（常数）


Fslope1=[0.66691593	0.462717789	0.331683997	0.207405118	0.194866743	0.134870057	0.148560329	0.109604617	0.058292019	0.020511147	-0.009661363	-0.02931673	-0.038332245	-0.107176106	-0.051233374	-0.202628134	-0.145166472	-0.077774246	-0.625016178	-0.05434	-0.26849199	-0.335829899	-0.424553469	-0.432361771	-0.47027074	-0.512607911	-0.593902792	-0.667215678	-0.749719581	-0.917080533	-1.032879535	-2.153511551	-1.912207855	-3.562871254	8.89383685	2.508189237];
Fslope2=[0.683469063	0.261728142	0.01368938	-0.244893622	-0.330281346	-0.384346533	-0.387074443	-0.368339205	-0.363819629	-0.323175767	-0.268647639	-0.144151469	-0.11708353	-0.027419827	0.026205738	0.078372178	0.134898622	0.194878128	0.386818667	0.480360673	0.287908169	-0.157270374	-0.173682839	-0.190313942	-0.38609733	-0.445219596	-0.476921704	-0.483632612	-0.574608539	-0.599510787	-0.809585466	-0.930937288	-1.445406372	-2.109060055	-40.65586901	1.631624387];
Fslope3=[0.670254449	0.575021686	0.436347142	0.330641986	0.193568203	0.132761164	0.033286955	-0.021860412	-0.111201749	-0.248714819	-0.270243127	-0.253614536	-0.286780982	-0.340627014	-0.269382345	-0.245328494	-0.232965342	-0.122189587	-0.269820848	-0.314666165	-0.294739251	-0.253627287	-0.247011177	-0.242545009	-0.262016774	-0.26603924	-0.272504589	-0.252438929	-0.223734863	-0.298145658	-0.194268025	-0.160984283	-0.128608779	-0.073647861	-0.080064084	-0.090065749];
Fslope4=[0.802217952	0.365308334	0.297477006	0.167338091	0.034288555	0.048098229	-0.016609892	-0.276854186	-0.381059606	-0.520346623	-0.491912309	-0.666805806	-0.711752569	-0.749828253	-0.610276172	-0.77107621	-0.683858214	-0.496312146	-0.552886189	-0.500084759	-0.45189655	-0.347294878	-0.234028898	-0.282469713	-0.318911925	-0.246600747	-0.155200138	-0.222312838	-0.200034442	-0.208383766	0.014082417	0.034043671	0.197345463	0.220963681	0.277574091	0.303415591];    
Fsim1=simcompare(Fslope,Fslope1);
Fsim2=simcompare(Fslope,Fslope2);
Fsim3=simcompare(Fslope,Fslope3);
Fsim4=simcompare(Fslope,Fslope4);
limit=[-1,1;-1,1;-1,1;-1,1];
prefer=[0,1;0,1;0,1;0,1];
draw_radar([Fsim1 Fsim2 Fsim3 Fsim4],limit,prefer,{'梭形','弓形','反S形','Z形'});
function Fsim=simcompare(Fslope,Fslope1)%相似度计算
sumsim1=0;%分子
for k=1:36
sumsim1=sumsim1+Fslope1(k)*Fslope(k);
end
sumsim1b=sum(Fslope.*Fslope)*sum(Fslope1.*Fslope1);
Fsim=sumsim1/sqrtm(sumsim1b);
end

function draw_radar(data,lim,prefer_range,labels) %雷达图
    n=length(data);
    adj_data=zeros(n,1);
    point=zeros(n,2);
    adj_preferl=zeros(n,1);
    preferl_point=zeros(n,2);
    adj_preferu=zeros(n,1);
    preferu_point=zeros(n,2);
    
    set(gca,'units','normal','pos',[0 0 1 1]);
    axis off
    axis equal
    hold on
    theta_last=pi/2;
    for i=1:n
        theta=2*pi/n*i+pi/2;
        plot([0,500*cos(theta)],[0,500*sin(theta)],'k-','linewidth',2);
        for j=1:4
           plot([j*125*cos(theta_last),j*125*cos(theta)],[j*125*sin(theta_last),j*125*sin(theta)],'--','linewidth',2,'color',[0.5,0.5,0.5]);
        end
        
        theta_last=theta;
        if data(i)<lim(i,1)
            adj_data(i)=0;
        elseif data(i)>lim(i,2)
            adj_data(i)=500;
        else
            adj_data(i)=(data(i)-lim(i,1))/(lim(i,2)-lim(i,1))*500;
        end
        point(i,1:2)=[adj_data(i)*cos(theta);adj_data(i)*sin(theta)];
        adj_preferl(i)=(prefer_range(i,1)-lim(i,1))/(lim(i,2)-lim(i,1))*500;
        preferl_point(i,1:2)=[adj_preferl(i)*cos(theta);adj_preferl(i)*sin(theta)];
        adj_preferu(i)=(prefer_range(i,2)-lim(i,1))/(lim(i,2)-lim(i,1))*500;
        preferu_point(i,1:2)=[adj_preferu(i)*cos(theta);adj_preferu(i)*sin(theta)];
        text_around(510*cos(theta),510*sin(theta),labels{i},theta,18);
    end
    
    plot([preferl_point(:,1);preferl_point(1,1)],[preferl_point(:,2);preferl_point(1,2)],'r-','linewidth',3);
    plot([preferu_point(:,1);preferu_point(1,1)],[preferu_point(:,2);preferu_point(1,2)],'b-','linewidth',3);
    for i=1:n
        theta=2*pi/n*i+pi/2;
        for j=1:4
            text_around(j*125*cos(theta),j*125*sin(theta),num2str(lim(i,1)+(lim(i,2)-lim(i,1))/4*j),theta+pi/2,12);
        end
    end
    plot([point(:,1);point(1,1)],[point(:,2);point(1,2)],'k-','linewidth',2);
    fill(point(:,1),point(:,2),[0.9 0.9 0.7])
    alpha(0.5);
    texts=findobj(gca,'Type','Text');
    minx=-300;
    maxx=300;
    miny=-300;
    maxy=300;
    for i=1:length(texts)
        rect=get(texts(i),'Extent');
        x=rect(1);
        y=rect(2);
        dx=rect(3);
        dy=rect(4);
        if x<minx
            minx=x;
        elseif x+dx>maxx
            maxx=x+dx;
        end
        if y<miny
            miny=y;
        elseif y+dy>maxy
            maxy=y+dy;
        end
    end
    axis([minx-50,maxx+50,miny-20,maxy+20]);
end

function text_around(x,y,txt,theta,fontsize)
    if nargin==4
        fontsize=10;
    end
    section=mod(theta+pi/12,2*pi);
    if section>pi+pi/6
        %上对齐
        if section>1.5*pi+pi/6
            %左对齐
            text(x,y,txt,'VerticalAlignment','cap','HorizontalAlignment','left','Fontsize',fontsize);
        elseif section>1.5*pi
            %中对齐
            text(x,y,txt,'VerticalAlignment','cap','HorizontalAlignment','center','Fontsize',fontsize);
        else
            %右对齐
            text(x,y,txt,'VerticalAlignment','cap','HorizontalAlignment','right','Fontsize',fontsize);
        end
    elseif section>pi
        %中、右对齐
        text(x,y,txt,'VerticalAlignment','middle','HorizontalAlignment','right','Fontsize',fontsize);
    elseif section>pi/6
        %下对齐
        if section>0.5*pi+pi/6
            %右对齐
            text(x,y,txt,'VerticalAlignment','bottom','HorizontalAlignment','right','Fontsize',fontsize);
        elseif section>0.5*pi
            %中对齐
            text(x,y,txt,'VerticalAlignment','bottom','HorizontalAlignment','center','Fontsize',fontsize);
        else
            %左对齐
            text(x,y,txt,'VerticalAlignment','bottom','HorizontalAlignment','left','Fontsize',fontsize);
        end
    else
        %中、左对齐
        text(x,y,txt,'VerticalAlignment','middle','HorizontalAlignment','left','Fontsize',fontsize);
    end
end


