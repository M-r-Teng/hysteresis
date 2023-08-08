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
for k=1:36
Fslope(k)=(FeatureP(k+1,2)-FeatureP(k,2))/(FeatureP(k+1,1)-FeatureP(k,1));
end
Fslope1=[5.004490093	2.722438523	1.992596305	1.523357198	1.484060845	1.311791443	1.348962666	1.246193138	1.123800627	1.041881328	0.980862171	0.943036523	0.926165744	0.806397364	0.902527116	0.663024457	0.746470971	0.855676185	0.230756978	0.896921297	0.576675309	0.497196613	0.403948706	0.396295294	0.360293683	0.322219714	0.254781666	0.199604842	0.143040303	0.043252991	-0.016173873	-0.365786373	-0.313235834	-0.561679502	-1.253362216	-2.326093537];
Fslope2=[5.318497705	1.70902917	1.02775876	0.606562974	0.503441363	0.444724967	0.441883679	0.461625884	0.466469581	0.511514985	0.576481868	0.748020306	0.79037641	0.946623909	1.053821919	1.170073377	1.311867778	1.484095974	2.261677898	2.8488234	1.808626517	0.728204614	0.704037866	0.680228997	0.442900117	0.383872738	0.354167925	0.348042625	0.270156963	0.250382314	0.105225499	0.035766419	-0.182140023	-0.356718762	-0.951987558	-4.166438852];
Fslope3=[5.065282737	3.706122488	2.548283258	1.987937634	1.480060938	1.306169785	1.068866258	0.957214484	0.799853179	0.601646725	0.57450173	0.595386734	0.554266055	0.49183925	0.575569416	0.606001958	0.622105611	0.782230046	0.575025329	0.521298755	0.544712573	0.595370507	0.603834863	0.60959964	0.584764989	0.579729866	0.571703566	0.59688425	0.634340951	0.540659161	0.674665952	0.722676207	0.77209325	0.862808164	0.851741975	0.834751712];
Fslope4=[9.112141218	2.151136383	1.846881904	1.401935261	1.071012009	1.101057128	0.967322977	0.566349566	0.448163418	0.315489488	0.340561364	0.199899828	0.168393164	0.142969315	0.242022974	0.129256883	0.18774846	0.336619505	0.287924391	0.333257996	0.377508611	0.484456025	0.620707589	0.559491019	0.516401484	0.604362908	0.731301732	0.636242325	0.666618832	0.655103334	1.028567128	1.070486978	1.491732005	1.567274402	1.76844999	1.871152405];
Fsim1=simcompare(Fslope,Fslope1)
Fsim2=simcompare(Fslope,Fslope2)
Fsim3=simcompare(Fslope,Fslope3)
Fsim4=simcompare(Fslope,Fslope4)

function Fsim=simcompare(Fslope,Fslope1)%相似度计算
sumsim1=0;%分子
for k=1:36
sumsim1=sumsim1+Fslope1(k)*Fslope(k);
end
sumsim1b=sum(Fslope.*Fslope)*sum(Fslope1.*Fslope1);
Fsim=sumsim1/sqrtm(sumsim1b);
end



