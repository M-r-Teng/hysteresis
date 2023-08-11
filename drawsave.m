clear;clc;
fid=fopen('SkeletonCurve.txt','wt');
coordinate=[0.0 0.0];
dlmwrite('SkeletonCurve.txt',coordinate,'delimiter','\t','newline','pc','precision','%.2f')
sta=fclose(fid);
outdata=load('SJ.txt');
LineNum=size(outdata,1);
EndCircle=1;
k=1;
%--------------------�Ҹ���ѭ���ָ��λ��----------------------��λ�ƣ���->��or��->��
for i=3:LineNum
if (outdata(i,1)*outdata(i-1,1)<=0)&&(outdata(i-1,1)<0) %--Ĭ�ϳ�ʼ���ش������ܿ�ʼ
       EndCirclePoints(EndCircle)=i;  
       EndCircle=EndCircle+1;                              
   end
end
%���������仯�ң��������Ҳ������
for i=3:LineNum
if (outdata(i,2)*outdata(i-1,2)<=0)&&(outdata(i-1,2)>0) %--Ĭ�ϳ�ʼ���ش������ܿ�ʼ
    canyu(k)=outdata(i,1);
    k=k+1;
   end
end

%--------------------�������Ƿ���ڲ������ͻػ�----------------------
if EndCirclePoints(EndCircle-1)<LineNum
   LoopNum=length(EndCirclePoints)+1;                    
else
   LoopNum=length(EndCirclePoints);
end
%-------------------------�ָ�����ͻػ�---------------------------
for k=1:LoopNum
   if  k==1
       LoopCircles(k)={outdata(1:EndCirclePoints(k),:)};
   elseif k<LoopNum
       LoopCircles(k)={outdata(EndCirclePoints(k-1):EndCirclePoints(k),:)};
   else
       LoopCircles(k)={outdata(EndCirclePoints(k-1):LineNum,:)};%���һ��Ӧ���������ڰ��ܡ���������
   end
end
%--------��ȡλ�ƿ��Ƶļ����ƶ�
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

%--------------------��ȡ�����ͻػ����عǼ����ߵ�---------------
for k=1:LoopNum %����ȥ���Ȧ�����
   A=LoopCircles{k};    
   if k<LoopNum  
        [ColMaxValue,LineMax]=max(A);         
        [ColMinValue,LineMin]=min(A);              
       SkeletonPointMax=A(LineMax(2),:);             
       SkeletonPointMin=A(LineMin(2),:);
       SkeletonPointsPositive(k,:)=SkeletonPointMax;      %�Ǽ����ߵ㴢���ھ���SkeletonPoints��
       SkeletonPointsNegative(k,:)=SkeletonPointMin;
       Gang(k)= (abs(SkeletonPointMax(2))+abs(SkeletonPointMin(2)))/(abs(SkeletonPointMax(1))+abs(SkeletonPointMin(1))); %��kȦ�ն��˻����������2Ȧ
%-------------------�������һ���������ܣ����ڰ��ܣ����Ƕ��ڰ���---------------
   else
       if EndCirclePoints(EndCircle-1)==LineNum   %������--------------------------------------------------
            [ColMaxValue,LineMax]=max(A);          
            [ColMinValue,LineMin]=min(A);         
           SkeletonPointMax=A(LineMax(2),:);         
           SkeletonPointMin=A(LineMin(2),:);
           SkeletonPointsPositive(k,:)=SkeletonPointMax; %�Ǽ����ߵ㴢���ھ���SkeletonPoints��
           SkeletonPointsNegative(k,:)=SkeletonPointMin;        
       elseif outdata(LineNum,1)<0             %�������ܣ���һ�����ֵ
            [ColMinValue,LineMin]=min(A);          %��ȡ��Сֵ�к�LineMin
           SkeletonPointMin=A(LineMin(2),:);
           SkeletonPointsPositive(k,:)=SkeletonPointMax;
           SkeletonPointsNegative(k,:)=SkeletonPointMin;
       end
   end
end


%�ӹǼ���������ϵ��
A=SkeletonPointsPositive;
x=A(:,1); %λ��
y=A(:,2); %��
a0=max(y);
[m0,n0]=find(A==a0); %��ͼ�з�ֵ��G,n0�ò���
B0=A(1:m0,:); %���ǰ���
x1=B0(:,1);
y1=B0(:,2);
lineb=polyfit(x1,y1,5); %5�׶���ʽ���
m0=m0(1,1);
c0=A(m0,1);
t=0:0.001:c0; 
d0=polyval(lineb,t); %������b��������ϣ�t����ֵ��
s=trapz(t,d0); %��OG���߰�Χ���λ������
u1=2*(B0(m0,1)*B0(m0,2)-s)/a0;  %�����Dy��ͼ��A1=A2����ͬʱ������������OBGDmaxҲ��ȣ�������=���߰�Χ�����s=B0(m0,1)*B0(m0,2)-u1*a/2��
%F1(k)=polyval(lineb,u1); %�ó���Ӧ������
e0=length(A);
C=A(m0:e0,:); %������Σ����е�m��BC����
x2=C(:,1);
y2=C(:,2);
linef=polyfit(x2,y2,5); %�ٵ������
k0=c0; %c0=B0(m0,1)
q=polyval(linef,k0); 
while abs(0.85*a0-q)/a0>0.01 %ѭ���ҵ��ݲ�1%���Ǹ������㡱
k0=k0+0.1;
q=polyval(linef,k0);
end
u2=k0; %��Du
YanxingSkeleton=u2/u1

for k=LoopNum:-1:2 %����-1��������һ��,ɸѡ�Ǽ���
   A=LoopCircles{k};
   B=LoopCircles{k-1};
    [ColMaxValueA,LineMaxA]=max(A);         
    [ColMinValueA,LineMinA]=min(A);          
    [ColMaxValueB,LineMaxB]=max(B);        
    [ColMinValueB,LineMinB]=min(B);  
   if abs(ColMaxValueA(1)-ColMaxValueB(1))<4 && abs(ColMaxValueA(2)-ColMaxValueB(2))<3000         %�ж�Ϊͬһ���ؼ����ҷǿ���λ�ƣ��ݲ�ȡΪ4mm
           SkeletonPointsPositive(k,:)=[]; %ͬ��ȡ��һ��
      
   end
   if abs(ColMinValueA(1)-ColMinValueB(1))<4 && abs(ColMaxValueA(2)-ColMaxValueB(2))<3000 
           SkeletonPointsNegative(k,:)=[];
      
   end
end
SkeletonNump=size(SkeletonPointsPositive,1);
SkeletonNumn=size(SkeletonPointsNegative,1);

SkeletonPointsPositive=sortrows(SkeletonPointsPositive,1);
SkeletonPointsNegative=sortrows(SkeletonPointsNegative,1); %��λ�ƴ�С��������
%�޳������к�Ȼ�ġ��������㣬��������ǡ����𡱣�����Ч���׳���������
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
dlmwrite('SkeletonCurve.txt',G);                           %�Ǽ�������ֵ�洢��SkeletonCurve.txt�ĵ���
d=G(:,1);
e=G(:,2);
a=outdata(:,1);
c=outdata(:,2);
plot(a,c,'r');
hold on;
plot(d,e,'s-k');                                       %�����ͻ�������Ǽ�����
currentFolder = pwd;%���浽��ǰ�ļ��У�δ����.m�ļ������ļ��У�
fileName = 'skeleton.jpg';
saveas(gcf,fullfile(currentFolder, fileName));
%�����
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
 sum=sum+f(j)*g(j+1)-f(j+1)*g(j); %���������
 end
 LoopArea(k)=-(sum+f(LoopNumA)*g(1)-f(1)*g(LoopNumA))/2;  %����������γɱջ�����ʵ֮ǰ����LoopʱҲ�ظ���������㣬����Ҳ���ⲻ��
 Nianzhi(k)=(LoopArea(k)/3.1416)/(ColMaxValueA(1,1)*ColMaxValueA(1,2)+ColMinValueA(1,1)*ColMinValueA(1,2));
 TotalE=TotalE+LoopArea(k);
 LeijiE(k)=TotalE;
end
x=1:1:LoopNum-1;
%title('Energy');
%xlabel('LoopNum');
close all; % �ر����е�ͼ���Ӵ�
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
Gang(2)=0;%��ʼ�ն�һ��������������ָն��˻����˴���Ϊ0���޳���Ҳ���Լ�һ���жϣ���������������̫����ɾȥ����ֵ����������
plot(x,Gang,'b:+');
fileName = 'Gang.jpg';
saveas(gcf,fullfile(currentFolder, fileName));
close all;

