clear;clc;
fid=fopen('SkeletonCurve.txt','wt');
coordinate=[0.0 0.0];
dlmwrite('SkeletonCurve.txt',coordinate,'delimiter','\t','newline','pc','precision','%.2f')
sta=fclose(fid);
outdata=load('SJ.txt');
LineNum=size(outdata,1);
EndCircle=1;
k=1;
%--------------------找各级循环分割点位置----------------------按位移，负->正or负->零
for i=3:LineNum
if (outdata(i,1)*outdata(i-1,1)<=0)&&(outdata(i-1,1)<0) %--默认初始加载从正半周开始
       EndCirclePoints(EndCircle)=i;  
       EndCircle=EndCircle+1;                              
   end
end
%按力正负变化找，现用于找残余变形
for i=3:LineNum
if (outdata(i,2)*outdata(i-1,2)<=0)&&(outdata(i-1,2)>0) %--默认初始加载从正半周开始
    canyu(k)=outdata(i,1);
    k=k+1;
   end
end

%--------------------检查最后是否存在不完整滞回环----------------------
if EndCirclePoints(EndCircle-1)<LineNum
   LoopNum=length(EndCirclePoints)+1;                    
else
   LoopNum=length(EndCirclePoints);
end
%-------------------------分割各个滞回环---------------------------
for k=1:LoopNum
   if  k==1
       LoopCircles(k)={outdata(1:EndCirclePoints(k),:)};
   elseif k<LoopNum
       LoopCircles(k)={outdata(EndCirclePoints(k-1):EndCirclePoints(k),:)};
   else
       LoopCircles(k)={outdata(EndCirclePoints(k-1):LineNum,:)};%最后一段应区分是少于半周、少于整周
   end
end
%--------提取位移控制的加载制度
for k=1:LoopNum-1 
   A=LoopCircles{k};            
        [ColMaxValue,LineMax]=max(A);         
        [ColMinValue,LineMin]=min(A);              
       WeiyiPositive(k)=ColMaxValue(1,1);
       WeiyiNegative(k)=ColMinValue(1,1);
       Weiyi(2*k-1)=ColMaxValue(1,1);
       Weiyi(2*k)=ColMinValue(1,1);
 end 
x=1:1:2*LoopNum-2;
plot(x,Weiyi,'r');

%--------------------提取各个滞回环荷载骨架曲线点---------------
for k=1:LoopNum %若舍去最后几圈则减几
   A=LoopCircles{k};    
   if k<LoopNum  
        [ColMaxValue,LineMax]=max(A);         
        [ColMinValue,LineMin]=min(A);              
       SkeletonPointMax=A(LineMax(2),:);             
       SkeletonPointMin=A(LineMin(2),:);
       SkeletonPointsPositive(k,:)=SkeletonPointMax;      %骨架曲线点储存在矩阵SkeletonPoints中
       SkeletonPointsNegative(k,:)=SkeletonPointMin;
       Gang(k)= (abs(SkeletonPointMax(2))+abs(SkeletonPointMin(2)))/(abs(SkeletonPointMax(1))+abs(SkeletonPointMin(1))); %第k圈刚度退化，不计最后2圈
%-------------------区分最后一段是完整周，少于半周，还是多于半周---------------
   else
       if EndCirclePoints(EndCircle-1)==LineNum   %完整周--------------------------------------------------
            [ColMaxValue,LineMax]=max(A);          
            [ColMinValue,LineMin]=min(A);         
           SkeletonPointMax=A(LineMax(2),:);         
           SkeletonPointMin=A(LineMin(2),:);
           SkeletonPointsPositive(k,:)=SkeletonPointMax; %骨架曲线点储存在矩阵SkeletonPoints中
           SkeletonPointsNegative(k,:)=SkeletonPointMin;        
       elseif outdata(LineNum,1)<0             %仅正半周，仅一个最大值
            [ColMinValue,LineMin]=min(A);          %获取最小值行号LineMin
           SkeletonPointMin=A(LineMin(2),:);
           SkeletonPointsPositive(k,:)=SkeletonPointMax;
           SkeletonPointsNegative(k,:)=SkeletonPointMin;
       end
   end
end


%从骨架线求延性系数
A=SkeletonPointsPositive;
x=A(:,1); %位移
y=A(:,2); %力
a0=max(y);
[m0,n0]=find(A==a0); %找图中峰值点G,n0用不上
B0=A(1:m0,:); %提出前半段
x1=B0(:,1);
y1=B0(:,2);
lineb=polyfit(x1,y1,5); %5阶多项式拟合
m0=m0(1,1);
c0=A(m0,1);
t=0:0.001:c0; 
d0=polyval(lineb,t); %计算在b拟合曲线上，t处的值。
s=trapz(t,d0); %算OG曲线包围梯形积分面积
u1=2*(B0(m0,1)*B0(m0,2)-s)/a0;  %反算此Dy。图中A1=A2，则同时加上曲边梯形OBGDmax也相等，即梯形=曲线包围面积。s=B0(m0,1)*B0(m0,2)-u1*a/2。
%F1(k)=polyval(lineb,u1); %得出对应屈服力
e0=length(A);
C=A(m0:e0,:); %提出后半段，其中点m是BC共用
x2=C(:,1);
y2=C(:,2);
linef=polyfit(x2,y2,5); %再单独拟合
k0=c0; %c0=B0(m0,1)
q=polyval(linef,k0); 
while abs(0.85*a0-q)/a0>0.01 %循环找到容差1%的那个“交点”
k0=k0+0.1;
q=polyval(linef,k0);
end
u2=k0; %此Du
YanxingSkeleton=u2/u1

for k=LoopNum:-1:2 %步长-1，倒过来一样,筛选骨架线
   A=LoopCircles{k};
   B=LoopCircles{k-1};
    [ColMaxValueA,LineMaxA]=max(A);         
    [ColMinValueA,LineMinA]=min(A);          
    [ColMaxValueB,LineMaxB]=max(B);        
    [ColMinValueB,LineMinB]=min(B);  
   if abs(ColMaxValueA(1)-ColMaxValueB(1))<4 && abs(ColMaxValueA(2)-ColMaxValueB(2))<3000         %判断为同一荷载级，且非控制位移，容差取为4mm
           SkeletonPointsPositive(k,:)=[]; %同级取第一个
      
   end
   if abs(ColMinValueA(1)-ColMinValueB(1))<4 && abs(ColMaxValueA(2)-ColMaxValueB(2))<3000 
           SkeletonPointsNegative(k,:)=[];
      
   end
end
SkeletonNump=size(SkeletonPointsPositive,1);
SkeletonNumn=size(SkeletonPointsNegative,1);

SkeletonPointsPositive=sortrows(SkeletonPointsPositive,1);
SkeletonPointsNegative=sortrows(SkeletonPointsNegative,1); %按位移从小到大排序
%剔除曲线中忽然的“沉降”点，若沉降点非“个别”，则无效，易超出索引。
for k=SkeletonNumn-1:-1:2
    n1=SkeletonPointsNegative(k,2);
    n2=SkeletonPointsNegative(k-1,2);
    n3=SkeletonPointsNegative(k+1,2);
    if n1 > n2 && n1>n3
    SkeletonPointsNegative(k,:)=[];
    end    
end
C=[SkeletonPointsPositive;SkeletonPointsNegative];
G=sortrows(C,1);
dlmwrite('SkeletonCurve.txt',G);                           %骨架曲线数值存储在SkeletonCurve.txt文档中
d=G(:,1);
e=G(:,2);
a=outdata(:,1);
c=outdata(:,2);
plot(a,c,'r');
hold on;
plot(d,e,'s-k');                                       %画出滞回曲线与骨架曲线
currentFolder = pwd;%保存到当前文件夹（未必是.m文件所在文件夹）
fileName = 'skeleton.jpg';
saveas(gcf,fullfile(currentFolder, fileName));
%求耗能
TotalE=0;
for k=1:LoopNum-1
 A=LoopCircles{k};
 [ColMaxValueA,LineMaxA]=max(A);
 [ColMinValueA,LineMinA]=min(A); 
 f=A(:,1);
 g=A(:,2);
 sum=0;
 LoopNumA=size(A,1);
 for j=1:LoopNumA-1
 sum=sum+f(j)*g(j+1)-f(j+1)*g(j); %求有向面积
 end
 LoopArea(k)=-(sum+f(LoopNumA)*g(1)-f(1)*g(LoopNumA))/2;  %加上最初的形成闭环，其实之前划分Loop时也重复用了新起点，不加也问题不大。
 Nianzhi(k)=(LoopArea(k)/3.1416)/(ColMaxValueA(1,1)*ColMaxValueA(1,2)+ColMinValueA(1,1)*ColMinValueA(1,2));
 TotalE=TotalE+LoopArea(k);
 LeijiE(k)=TotalE;
end
x=1:1:LoopNum-1;
%title('Energy');
%xlabel('LoopNum');
close all; % 关闭所有的图形视窗
currentFolder = pwd;
plot(x,LoopArea,'b:+');
fileName = 'Energy.jpg';
saveas(gcf,fullfile(currentFolder, fileName));
close all; 
plot(x,LeijiE,'b:+');
fileName = 'TotalEnergy.jpg';
saveas(gcf,fullfile(currentFolder, fileName));
close all; 
plot(x,Nianzhi,'r:+');
fileName = 'Nianzhi.jpg';
saveas(gcf,fullfile(currentFolder, fileName));
close all; 
plot(x,canyu,'b:+');
fileName = 'canyu.jpg';
saveas(gcf,fullfile(currentFolder, fileName));
close all; 
x=1:1:LoopNum-2;
Gang(1)=0;
Gang(2)=0;%初始刚度一般过大，且难以体现刚度退化。此处设为0或剔除。也可以加一个判断，如果超出后边数据太多则删去，数值合理则保留。
plot(x,Gang,'b:+');
fileName = 'Gang.jpg';
saveas(gcf,fullfile(currentFolder, fileName));
close all;

