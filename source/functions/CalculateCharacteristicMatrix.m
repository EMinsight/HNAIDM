function [Mxx,Myx,Mxy,Myy,Mxz,Myz] = CalculateCharacteristicMatrix(f_calc,g_calc,fgSquare,NA,indexImage)
% 计算偏振传输矩阵M
% NA/indexImage对sin(theta_obj)
%《optical imaging in projection microlithography》式（6.13）
alpha = (NA/indexImage)*f_calc;  % x方向余弦
beta = (NA/indexImage)*g_calc;   % y方向余弦
gamma = sqrt(1-(NA/indexImage)^2*fgSquare);%z方向余弦
vectorRho2 = 1-gamma.^2;

Pxsx = (beta.^2)./vectorRho2; %偏振光向s方向映射后在y方向的投影分量
Pysx = (-alpha.*beta)./vectorRho2;
Pxsy = Pysx;
Pysy = (alpha.^2)./vectorRho2;
%         Pxsz = 0;
%         Pysz = 0;

Pxpx = gamma.*Pysy;  %偏振光在向p方向映射后在y方向的投影分量
Pypx = -gamma.*Pysx;
Pxpy = Pypx;
Pypy = gamma.*Pxsx;
Pxpz = -alpha;
Pypz = -beta;

centerpoint = find(fgSquare<eps);
if length(centerpoint)>eps
    Pxsx(centerpoint) = 0;
    Pysx(centerpoint) = 0;
    Pxsy(centerpoint) = 0;
    Pysy(centerpoint) = 0;
    
    Pxpx(centerpoint) = 0;
    Pypx(centerpoint) = 0;
    Pxpy(centerpoint) = 0;
    Pypy(centerpoint) = 0;
    Pxpz(centerpoint) = 0;
    Pypz(centerpoint) = 0;
end

Mxx = Pxsx + Pxpx;
Myx = Pysx + Pypx;
Mxy = Pxsy + Pxpy;
Myy = Pysy + Pypy;
Mxz = Pxpz;
Myz = Pypz;


end
