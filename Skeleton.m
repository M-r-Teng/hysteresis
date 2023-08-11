%按峰值荷载取点，再筛去“不合理”的突变点
clear;clc;
fid=fopen('SkeletonCurve.txt','wt');
coordinate=[0.0 0.0];
dlmwrite('SkeletonCurve.txt',coordinate,'delimiter','\t','newline','pc','precision','%.2f')
sta=fclose(fid);
outdata=load('SJ.txt');
%outdata = unique(outdata, 'rows');
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

%close all;

%--------------------提取各个滞回环荷载骨架曲线点---------------
choice = input('请输入1、2、3或4,选择要运行的程序部分:峰值荷载法（优化版），峰值位移法，峰值极半径法；凹包法');

switch choice
case 1
disp('运行峰值荷载法（优化版）');
for k=1:LoopNum
   A=LoopCircles{k};    
   if k<LoopNum  
        [ColMaxValue,LineMax]=max(A);         
        [ColMinValue,LineMin]=min(A);              
       SkeletonPointMax=A(LineMax(2),:);             
       SkeletonPointMin=A(LineMin(2),:);
       SkeletonPointsPositive(k,:)=SkeletonPointMax;      %骨架曲线点储存在矩阵SkeletonPoints中
       SkeletonPointsNegative(k,:)=SkeletonPointMin;
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

case 2 
       disp('运行峰值位移法');
for k=1:LoopNum
   A=LoopCircles{k};    
   if k<LoopNum  
        [ColMaxValue,LineMax]=max(A);         
        [ColMinValue,LineMin]=min(A);              
       SkeletonPointMax=A(LineMax(1),:);             
       SkeletonPointMin=A(LineMin(1),:);
       SkeletonPointsPositive(k,:)=SkeletonPointMax;      %骨架曲线点储存在矩阵SkeletonPoints中
       SkeletonPointsNegative(k,:)=SkeletonPointMin;
%-------------------区分最后一段是完整周，少于半周，还是多于半周---------------
   else
       if EndCirclePoints(EndCircle-1)==LineNum   %完整周--------------------------------------------------
            [ColMaxValue,LineMax]=max(A);          
            [ColMinValue,LineMin]=min(A);         
           SkeletonPointMax=A(LineMax(1),:);         
           SkeletonPointMin=A(LineMin(1),:);
           SkeletonPointsPositive(k,:)=SkeletonPointMax; %骨架曲线点储存在矩阵SkeletonPoints中
           SkeletonPointsNegative(k,:)=SkeletonPointMin;        
       elseif outdata(LineNum,1)<0             %仅正半周，仅一个最大值
            [ColMinValue,LineMin]=min(A);          %获取最小值行号LineMin
           SkeletonPointMin=A(LineMin(1),:);
           SkeletonPointsPositive(k,:)=SkeletonPointMax;
           SkeletonPointsNegative(k,:)=SkeletonPointMin;
       end
   end
end

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

    
case 3
    disp('运行峰值极半径法');
    
    for k=1:LoopNum
   A=LoopCircles{k};    
   max_x = max(abs(A(:,1)));
max_y = max(abs(A(:,2)));

% 归一化
A_norm = A ./ [max_x max_y]; 
   % 分别提取第一、第三象限的数据
quad1 = A_norm(A(:,1)>=0 & A(:,2)>=0, :);
quad3 = A_norm(A(:,1)<=0 & A(:,2)<=0, :);
% 计算距离并找出最大值
dist1 = sqrt(quad1(:,1).^2 + quad1(:,2).^2);
[maxDist1, idx1] = max(dist1);
dist3 = sqrt(quad3(:,1).^2 + quad3(:,2).^2); 
[maxDist3, idx3] = max(dist3);
quad1 = A(A(:,1)>=0 & A(:,2)>=0, :);
quad3 = A(A(:,1)<=0 & A(:,2)<=0, :);%对应回原来的未归一化的值
    if idx1 >0
       SkeletonPointsPositive(k,:)=quad1(idx1,:); 
    end
    if idx3 >0%idx可能不存在
       SkeletonPointsNegative(k,:)=quad3(idx3,:); 
    end 
    end

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
    
case 4
disp('运行凹包法');
angleThreshold = pi/2;  % 设置夹角阈值为90度
A=unique(outdata, 'rows');
 Posi = A(A(:,1)>=0 & A(:,2)>=0, :);%或设置A(:,1)>=峰值承载力对应位移
 Nega = A(A(:,1)<=0 & A(:,2)<=0, :);
SkeletonPointsPositive = concaveHullFilter(Posi, angleThreshold);
SkeletonPointsNegative = concaveHullFilter(Nega, angleThreshold);
   
otherwise
    disp('输入有误,请输入1、2、3或4');
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
fileName = 'Skeleton.jpg';
currentFolder = pwd;
saveas(gcf,fullfile(currentFolder, fileName));
%close all;

function [concaveHull] = concaveHullFilter(A, angleThreshold)
sortedA = unique(A, 'rows');   

    % 遍历排序后的点，计算并检查向量夹角
     add=1;
 while add ==1
     add =0;
        for i = size(sortedA,1):-1:3 
             v1 = sortedA(i, :)-sortedA(i-1, :);
             v2 = sortedA(i-2, :)-sortedA(i-1, :);

        cosTheta = dot(v1, v2) / (norm(v1) * norm(v2));
        theta = acos(cosTheta);
%         
%         % 如果夹角大于阈值，暂时保留，小于阈值，则忽略内部点
        if theta <  angleThreshold
        
[minDist, minId] = min([abs(sortedA(i, 2)),abs(sortedA(i-1, 2)), abs(sortedA(i-2, 2))]);
if minId == 2
    sortedA(i-1, :)=[];
add=1;  
end    
        end
        end
 end
 [theta, rho] = cart2pol(sortedA(:,1), sortedA(:,2));
    [~, idx] = sort(theta);
    sortedA = sortedA(idx,:);%避免出现连续两个内部点，使得角度满足限值，但是看起来不像包络
 add=1;
 while add ==1
     add =0;
        for i = size(sortedA,1):-1:3 
             v1 = sortedA(i, :)-sortedA(i-1, :);
             v2 = sortedA(i-2, :)-sortedA(i-1, :);

        cosTheta = dot(v1, v2) / (norm(v1) * norm(v2));
        theta = acos(cosTheta);
%         
%         % 如果夹角大于阈值，暂时保留，小于阈值，则忽略内部点
        if theta <  angleThreshold
        
[minDist, minId] = min([abs(sortedA(i, 2)),abs(sortedA(i-1, 2)), abs(sortedA(i-2, 2))]);
if minId == 2
    sortedA(i-1, :)=[];
add=1;  
end    
        end
        end
 end
 concaveHull=sortedA
end