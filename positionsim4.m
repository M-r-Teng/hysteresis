%拟合+获取特征斜率
clear;clc;
LoopCircles=load('SJ4.txt');

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
Fpos=FeatureP(:,2);
%圆，弓，反S，Z的位置（常数）
Fpos1=[-0.664380947 -0.38635372 -0.235107135 -0.12440734 -0.039776385 0.04267144 0.115548742 0.190491112 0.259724065 0.322157433 0.380039729 0.434532072 0.486922989 0.538376642 0.583176495 0.633316891 0.670151583 0.711622192 0.759159758 0.77197959 0.821808551 0.853846068 0.881468102 0.903909697 0.925926102 0.945942418 0.963843513 0.97799805 0.989087208 0.997033892 0.999436836 0.998538287 0.978216822 0.960814831 0.929610415 0.85997918 0.730751762];
Fpos2=[-0.711716307 -0.416244212 -0.321298147 -0.264200438 -0.230502495 -0.202533531 -0.177826588 -0.153277495 -0.127631612 -0.101716636 -0.073299137 -0.041272366 0.000284318 0.044194118 0.096784335 0.155329997 0.220334074 0.293215617 0.375665394 0.501314166 0.659582132 0.760061383 0.800517195 0.83963041 0.87742091 0.902026472 0.923352735 0.943028731 0.962364432 0.977373152 0.991283281 0.997129142 0.999116165 0.988997275 0.969179566 0.916291368 0.684822543];
Fpos3=[-1 -0.718595404 -0.51269971 -0.371128418 -0.260687438 -0.17846183 -0.105896842 -0.046515383 0.006663199 0.051099487 0.084524305 0.116441068 0.149518108 0.180310667 0.20763507 0.239611148 0.273277924 0.307839347 0.351296572 0.383242423 0.412203465 0.442465275 0.475541414 0.509087795 0.542954442 0.575441386 0.607648601 0.63940991 0.672570146 0.70781131 0.73784793 0.775329372 0.81547805 0.858372119 0.906305906 0.953624905 1];
Fpos4=[-1 -0.493769932 -0.374262356 -0.271657805 -0.193772513 -0.134271846 -0.073102005 -0.01936184 0.012102025 0.036999992 0.054527186 0.073447262 0.084552808 0.093907984 0.101850724 0.115296444 0.122477382 0.132907852 0.151608936 0.167604735 0.186119068 0.207091769 0.234005993 0.268489748 0.299572582 0.328261553 0.36183727 0.402465144 0.43781194 0.47484632 0.51124095 0.568383568 0.627855067 0.710729067 0.797799867 0.896047089 1];

Fsim1=simcompare(Fpos,Fpos1)
Fsim2=simcompare(Fpos,Fpos2)
Fsim3=simcompare(Fpos,Fpos3)
Fsim4=simcompare(Fpos,Fpos4)

function Fsim=simcompare(Fpos,Fpos1)%相似度计算
sumsim1=0;%分子
for k=1:37
sumsim1=sumsim1+Fpos1(k)*Fpos(k);
end
sumsim1b=sum(Fpos.*Fpos)*sum(Fpos1.*Fpos1);
Fsim=sumsim1/sqrtm(sumsim1b);
end




