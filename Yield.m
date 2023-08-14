
clear;clc;
Skeleton=load('SkeletonCurve.txt');
y_mid = (max(Skeleton(:,2)) + min(Skeleton(:,2)))/2;
A(:,2) = Skeleton(:,2) - y_mid;
A(:,1) = Skeleton(:,1);%也可以再居中
A = A(A(:,1)>0 & A(:,2)>0, :);
x=A(:,1); %位移
y=A(:,2); %力
u0=max(x);
a0=max(y);
[m0,n0]=find(A==a0); %找峰值点
Pr0=A(1:m0,:); %提出前半段
lineb=polyfit(Pr0(:,1),Pr0(:,2),5); %5阶多项式拟合前半段，前半段骨架曲线不可太稀疏
m0=m0(1,1);%是必要的，因为“峰值”可能不止一个

c0=A(m0,1);%峰值点对应位移
t=0:0.001*c0:c0; 
d0=[t; polyval(lineb,t)].'; %计算在b拟合曲线上，t处的值。

choice = input('请输入1、2、3,选择连线端点:峰值荷载点，峰值位移点，（峰值位移，峰值荷载）点');
switch choice
    case 1
        Max=[c0,a0];
    case 2
   [max_x, idx] = max(A(:,1)); 
   Max = A(idx,:);
    case 3
        n0=n0(1,1);u0=u0(1,1);
        Max=[u0,a0];
end

%计算最远点
k = Max(2)/Max(1);
dist = abs(d0(:,2) - k*d0(:,1)) ./ sqrt(1+k^2);
[maxDist, idx] = max(dist);
u1 = d0(idx,1)

e0=length(A);
P=A(m0:e0,:); %提出后半段，其中点m是BC共用
x2=P(:,1);
y2=P(:,2);
linef=polyfit(x2,y2,5); %再单独拟合
k0=c0; %c0=B0(m0,1)
q=polyval(linef,k0); 
while abs(0.85*a0-q)/a0>0.01 %循环找到容差1%的那个“交点”
k0=k0+0.01*c0;
q=polyval(linef,k0);
end
u2=k0 %此Du
YanxingSkeleton=u2/u1
