function [Mxx,Myx,Mxy,Myy,Mxz,Myz] = CalculateCharacteristicMatrix(f_calc,g_calc,fgSquare,NA,indexImage)
% ����ƫ�������M
% NA/indexImage��sin(theta_obj)
%��optical imaging in projection microlithography��ʽ��6.13��
alpha = (NA/indexImage)*f_calc;  % x��������
beta = (NA/indexImage)*g_calc;   % y��������
gamma = sqrt(1-(NA/indexImage)^2*fgSquare);%z��������
vectorRho2 = 1-gamma.^2;

Pxsx = (beta.^2)./vectorRho2; %ƫ�����s����ӳ�����y�����ͶӰ����
Pysx = (-alpha.*beta)./vectorRho2;
Pxsy = Pysx;
Pysy = (alpha.^2)./vectorRho2;
%         Pxsz = 0;
%         Pysz = 0;

Pxpx = gamma.*Pysy;  %ƫ�������p����ӳ�����y�����ͶӰ����
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
