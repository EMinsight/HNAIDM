classdef Source
    properties
        PntNum = 101;
        Wavelength = 193.368;
        Shape = 'annular';  
        %Sourc shape =>  'annular'  'multipole' 'dipolecirc'  'quadrupole' 'MMA' 'pixel' 'quasar'
        SigmaIn = 0.5;
        SigmaOut = 0.7;
        SigmaCenter = 0.5;
        SigmaRadius = 0.1;
        OpenAngle = pi/6;
        RotateAngle = 0;
        PoleNumber = 2;
        
        AdvancedParametersEnable = 0;
        AdvancedParameters = [];
        Source_Mask = [];
        %点光源参数
        SPointX
        SPointY
        SPointValue 
        % PSF parameter
        PSFEnable
        PSFSigma
        %% 光源偏振
        PolarizationType 
        %Polorization mode  'x_pol'   'y_pol'  'unpolarized'   'r_pol'  't_pol' 'fun' 
        %TODO 'line' 
        PolarizationParameters
        
        %% MMA model
        MMA_Number = 128;
        MMA_Coordinate = 'Polar'; %'Cartesian'
        MMA_Rho = sqrt(rand(1,128));
        MMA_Theta =  pi/2*rand(1,128);
        MMA_X =  rand(1,128)*0.5;
        MMA_Y =  rand(1,128)*0.5;
        
        MMA_PSFMAX = 0.01;
        MMA_PSFSigma = 0.03;
        MMA_PSFType ='gaussian'; % gaussian / sinc / none
        MMA_Normalization = 0;
        MMA_SourceConstructType = 'quarter'; % full / quarter
    end
    
    %% functions
    methods
        function s = Source()
            s.SigmaOut = 0.8;
            s.SigmaIn =0.6;
            % 初始化高级参数
            s.AdvancedParameters.A = []; %相对光强最大值，为null则不参与计算
            s.AdvancedParameters.BG = 0; %背景光强
            s.AdvancedParameters.sli = 50;%内环斜率参数
            s.AdvancedParameters.slo = 50;%外环斜率参数
            s.AdvancedParameters.sla = 50;%角向光强斜率
            s.AdvancedParameters.cA(1) = 0;
            s.AdvancedParameters.sA(1) = 0;
            s.AdvancedParameters.cA(2) = 0;
            s.AdvancedParameters.sA(2) = 0;
            
            s.AdvancedParameters.x0 = 0; %环中心x偏移
            s.AdvancedParameters.y0 = 0; %环中心y位置
            s.AdvancedParameters.dx0 = 0;
            s.AdvancedParameters.dy0 = 0;
            s.AdvancedParameters.cri2 =0;
            s.AdvancedParameters.sri2 = 0;
            s.AdvancedParameters.cro2 = 0;
            s.AdvancedParameters.sro2 = 0;
            
            % 光源Blur参数初始化
            s.PSFEnable = false;
            s.PSFSigma = 0.02;
            
            % 偏振参数初始化
            s.PolarizationType  = 't_pol';
            s.PolarizationParameters.Degree = 1;%尚不支持部分偏振光
            s.PolarizationParameters.Angle = 0; %偏振方向 x轴正半轴起逆时针
            s.PolarizationParameters.PhaseDifference = 0; %
            % 仅在PolarizationType参数设置fun时有用
            s.PolarizationParameters.PolFun_X = []; % xpol = fun(theta,rho)
            s.PolarizationParameters.PolFun_Y = [];% ypol = fun(theta,rho)
            
            %             % MMA initialization
            %             s.MMA_Number = 128;
            %             s.MMA_Rho = rand(1,128);
            %             s.MMA_Theta =  pi/2*rand(1,128);
            %             s.MMA_PSFMAX = 0.01;
            %             s.MMA_PSFSigma = 0.03;
            %             s.MMA_PSFType ='gaussian'; % gaussian / sinc / none
            %             s.MMA_Normalization = 0;
            %             s.MMA_SourceConstructType = 'quarter'; % full / quarter
        end
        
        function [Source_Data,Source_Matrix,X,Y] = Calc_SourceSimple(s) %去除值较小后的序列化光源点
            [ Source_Data, Source_Matrix,X,Y] = s.Calc_SourceAll();
            Low_Weight = Source_Data.Value<1e-5;%remove invalid point
            % r2 = Source_Data(:,1).^2 + Source_Data(:,2).^2;
            % delRaw = (r2>1) | Low_Weight;
            Source_Data.X(Low_Weight)=[];        
            Source_Data.Y(Low_Weight)=[];        
            Source_Data.Value(Low_Weight)=[]; 
        end
        
        function [Source_Data,Source_Matrix,X,Y] = Calc_SourceValid(s) % seralize source point
            [Source_Data,Source_Matrix,X,Y] = s.Calc_SourceAll();
                       
            Radius = Source_Data.X.^2 + Source_Data.Y.^2;
            Radius_Valid = find(Radius<=1);
            
            Source_Data.X=Source_Data.X(Radius_Valid);%
            Source_Data.Y=Source_Data.Y(Radius_Valid);%
            Source_Data.Value=Source_Data.Value(Radius_Valid);%
        end
        
        function [Source_Data,Source_Matrix,X,Y] = Calc_SourceAll(s) 

            Source_Coordinate_X = linspace(-1,1,s.PntNum);
            
            %% illumination mode
            switch lower(s.Shape) 
                case 'pixel'%像素化点光源照明
                    % 点光源需要对光源网格进行检查以保证网格尺寸符合定义
                    % FIXME:光瞳坐标和数值需要进行定义并且对数据格式进行检查

                    Source_Data.X = reshape(s.SPointX,[],1);%第一列存储x坐标
                    Source_Data.Y = reshape(s.SPointY,[],1);%第二列存储y坐标
                    Source_Data.Value = reshape(s.SPointValue,[],1);%第三列存储光强值
                    if ~isequal(length(Source_Data.X), s.PntNum.^2)
                        error('Source Matrix Size X and difined PntNum are not matched ');
                    end
                    if ~isequal(length(Source_Data.X), s.PntNum.^2)
                        error('Source Matrix Size Y and difined PntNum are not matched ');
                    end
                    if ~isequal(length(Source_Data.Value), s.PntNum.^2)
                        error('Source Matrix Size Value and difined PntNum are not matched ');
                    end
                    X = s.SPointX;
                    Y = s.SPointY;
                    Source_Matrix = s.SPointValue;

                case 'mma' % Micro-Mirror-Array model
                    if  strcmpi(s.MMA_Coordinate,'polar' )
                        % % %    [ Source_Matrix, X, Y ] = s.CalculateSourceMatrixWithMMA(s.MMA_Number,s.MMA_Rho, s.MMA_Theta,  Source_Coordinate_X, Source_Coordinate_X,s.MMA_PSFType, s.MMA_PSFSigma,s.MMA_PSFMAX, s.MMA_SourceConstructType);
                        [ Source_Matrix, X, Y ] = CalculateSourceMatrixWithMMA_V3(s.MMA_Number,s.MMA_Rho, s.MMA_Theta,  Source_Coordinate_X, Source_Coordinate_X,s.MMA_PSFType, s.MMA_PSFSigma,s.MMA_PSFMAX, s.MMA_SourceConstructType);
                    elseif  strcmpi(s.MMA_Coordinate,'Cartesian' )
                        [ Source_Matrix, X, Y ] = CalculateSourceMatrixWithMMA_Cartesian (s.MMA_Number,s.MMA_X, s.MMA_Y,  Source_Coordinate_X, Source_Coordinate_X,s.MMA_PSFType, s.MMA_PSFSigma,s.MMA_PSFMAX, s.MMA_SourceConstructType);
                    else
                        error('Unsupported coordinate system')
                    end
                    if s.MMA_Normalization
                        Source_Matrix = Source_Matrix/max(Source_Matrix(:));
                    end
                case 'annular'  %环形照明
                    if s.AdvancedParametersEnable %是否开启高级参数设置
                        [Source_Matrix,X,Y]=CalculateAnnularSourceMatrixWithAP(s.SigmaOut,s.SigmaIn,Source_Coordinate_X,Source_Coordinate_X,s.AdvancedParameters);
                    else
                        [Source_Matrix,X,Y]=CalculateAnnularSourceMatrix(s.SigmaOut,s.SigmaIn,Source_Coordinate_X,Source_Coordinate_X);
                    end
                case 'quasar'  %四极扇形照明
                    
                    openAngle = s.OpenAngle;
                    rotateAngle = s.RotateAngle;
                    if (rotateAngle>pi/2) ||  (rotateAngle<-pi/2)
                        error('error: roate angle must be in the range of [-pi/2,pi/2] ');
                    end
                    if (openAngle<=0 ) || openAngle>=pi/2
                        error('error: open angle must be in the range of [0,pi/2] ');
                    end
                    if s.AdvancedParametersEnable %是否开启高级参数设置
                        [Source_Matrix,X,Y]=CalculateQuasarSourceMatrixWithAP(s.SigmaOut,s.SigmaIn,openAngle,rotateAngle,Source_Coordinate_X,Source_Coordinate_X,s.AdvancedParameters);
                    else
                        [Source_Matrix,X,Y]=CalculateQuasarSourceMatrix(s.SigmaOut,s.SigmaIn,openAngle,Source_Coordinate_X,Source_Coordinate_X);
                    end
                    
                case 'dipolecirc' %二极扇形照明
                    openAngle = s.OpenAngle;
                    rotateAngle = s.RotateAngle;
                    if (rotateAngle>pi/2) ||  (rotateAngle<-pi/2)
                        error('error: roate angle must be in the range of [-pi/2,pi/2] ');
                    end
                    if (openAngle<0 ) || openAngle>pi
                        error('error: open angle must be in the range of [0,pi] ');
                    end
                    
                    if s.AdvancedParametersEnable %是否开启高级参数设置
                        [Source_Matrix,X,Y] = CalculateDipoleSourceMatrixWithAP(s.SigmaOut,s.SigmaIn,openAngle,rotateAngle,Source_Coordinate_X,Source_Coordinate_X,s.AdvancedParameters);
                    else
                        [Source_Matrix,X,Y] = CalculateDipoleSourceMatrix(s.SigmaOut,s.SigmaIn,openAngle,rotateAngle,Source_Coordinate_X,Source_Coordinate_X);
                    end                   
                    
                case 'multipole' %多极照明
                    rotateAngle = s.RotateAngle;
                    if (rotateAngle>pi/2) ||  (rotateAngle<-pi/2)
                        error('error: roate angle must be in the range of [-pi/2,pi/2] ');
                    end
                    [Source_Matrix,X,Y] = CalculateMultiCircSourceMatrix( s.SigmaCenter,s.SigmaRadius,s.PoleNumber,rotateAngle,Source_Coordinate_X,Source_Coordinate_X );
                    if s.AdvancedParametersEnable %是否开启高级参数设置
                        warning('Advanced Parameters are not supported for multipole source !');
                    end
                otherwise
                    error( 'unsupported illumination');
            end

            
            [sizeX, sizeY] = size(Source_Matrix);

            if ~isempty( s.Source_Mask )
                Source_Matrix = Source_Matrix .* s.Source_Mask;
            end
            % 数据序列化和添加光源blur
            if sizeX==sizeY
                % 添加光源Blur
                if s.PSFEnable
                    %                     s.PntNum = 501; % 测试参数
                    %                     s.PSFSigma = 0.02;
                    kernelSize = round(s.PntNum/10)*2+1;
                    kernelEdge = 1/(s.PntNum-1)*(kernelSize-1);
                    [kernelX,KernelY] = meshgrid(linspace(-kernelEdge,kernelEdge,kernelSize));
                    kernel = 1/sqrt(2*pi)/s.PSFSigma *exp(-(kernelX.^2+KernelY.^2)/s.PSFSigma.^2);
                    kernel(:,all(kernel<1e-6,1)) = [];
                    kernel(all(kernel<1e-6,2),:) = [];
                    
                    kernel = kernel/sum(kernel(:));
                    Source_Matrix = conv2(Source_Matrix,kernel,'same');
                end
                %中心点必须置0
                Source_Matrix((sizeX+1)/2,(sizeY+1)/2 ) = 0;
                Source_Data = ConvertSourceMatrix2SourceData(Source_Matrix,X,Y);
            end

            %             Radius = Source_Data.X.^2+Source_Data.Y.^2;
            %             Radius_Out=(Radius-1)>1e-10;
            %             Source_Data.X(Radius_Out) = [];
            %             Source_Data.Y(Radius_Out) = [];
            %             Source_Data.Value(Radius_Out) = [];
            %             Low_Weight=Source_Data(:,3)<1e-5;
            %             Source_Data(Low_Weight,3)=0;%
        end

        %% 光源偏振
        function [ PolarizedX, PolarizedY ] = Calc_PolarizationMap(s,theta,rho)
            %% 目前矢量模型无法计算完全非偏振光照明
                switch s.PolarizationType
                    case 'x_pol'
                        PolarizedX = ones(size(theta));
                        PolarizedY = zeros(size(theta));
                    case 'y_pol'
                        PolarizedX = zeros(size(theta));
                        PolarizedY = ones(size(theta));
                    case 'r_pol'
                        PolarizedX = cos(theta);
                        PolarizedY = sin(theta);
                    case 't_pol'
                        PolarizedX = sin(theta);
                        PolarizedY = -cos(theta);
                    case 'line_pol'
                        PolarizedX = sin(s.PolarizationParameters.Angle) .*ones(size(theta));
                        PolarizedY = cos(s.PolarizationParameters.Angle).*ones(size(theta));
                    case 'fun'
                        PolarizedX = s.PolarizationParameters.PolFun_X(theta,rho);
                        PolarizedY = s.PolarizationParameters.PolFun_Y(Theta,rho);
                    otherwise
                end
            
            biz = rho<eps;
            if length(biz)>eps
                PolarizedX(biz) = 0;
                PolarizedY(biz) = 0;
            end
        end
    end

end

% functions for generate source shape
% MMA卷积模型光源 连续MMA中心点 速度较快
function [ Source_Matrix, X, Y ] = CalculateSourceMatrixWithMMA_V3(MMA_Number,MMA_Rho, MMA_Theta,  Source_Coordinate_X, Source_Coordinate_Y,PSFType, PSFSigma,PSFMAX, SourceConstructType)

[ X, Y ] = meshgrid( single(Source_Coordinate_X), single(Source_Coordinate_Y ));
Source_Matrix = single(zeros(size(X)));

[ core_X, core_Y ] = pol2cart( MMA_Theta, MMA_Rho );
pntNum = length(Source_Coordinate_X);
core_Yr = round(core_X*(pntNum+1)/2+pntNum/2);
core_Xr = round(core_Y*(pntNum+1)/2+pntNum/2);
switch lower(PSFType)
    case 'gaussian'
        PSFRange = ceil(PSFSigma*length(Source_Coordinate_X)+1);
        PSFSigma2 = single(PSFSigma.^2);
        for iMMA =1:MMA_Number
            
            Xcol = max(core_Xr(iMMA)-PSFRange,1):min(core_Xr(iMMA)+PSFRange,pntNum);
            Yrow = max(core_Yr(iMMA)-PSFRange,1):min(core_Yr(iMMA)+PSFRange,pntNum);
            
            % 边界光源点单独处理
            %                 if  core_Xr(iMMA) == (pntNum+1)/2 || core_Yr(iMMA) == (pntNum+1)/2
            %                     Source_Matrix(Xcol, Yrow) = Source_Matrix(Xcol,Yrow)+ exp(-( (X(Xcol,Yrow)-core_X(iMMA)).^2+(Y(Xcol,Yrow)-core_Y(iMMA)).^2)/PSFSigma2)/2;
            %                 else
            %                     Source_Matrix(Xcol, Yrow) = Source_Matrix(Xcol,Yrow)+ exp(-( (X(Xcol,Yrow)-core_X(iMMA)).^2+(Y(Xcol,Yrow)-core_Y(iMMA)).^2)/PSFSigma2);
            %                 end
            Source_Matrix(Xcol, Yrow) = Source_Matrix(Xcol,Yrow)+ exp(-( (X(Xcol,Yrow)-core_X(iMMA)).^2+(Y(Xcol,Yrow)-core_Y(iMMA)).^2)/PSFSigma2);
            %                 %调试代码
            %                 xx = linspace(-1,1,101);
            %                 [ xx, yy ] = meshgrid( xx );
            %                 pcolor(xx,yy,Source_Matrix);hold
            % %                 scatter(core_X(iMMA),core_Y(iMMA));
            % %                 Source_Matrix(core_Xr(iMMA),core_Yr(iMMA)) = 1;
        end
    case 'gaussianrect'
        PSFRange = ceil(PSFSigma*length(Source_Coordinate_X)+1);
        PSFSigma2 = single(PSFSigma.^2);
        for iMMA =1:MMA_Number
            
            Xcol = max(core_Xr(iMMA)-PSFRange,1):min(core_Xr(iMMA)+PSFRange,pntNum);
            Yrow = max(core_Yr(iMMA)-PSFRange,1):min(core_Yr(iMMA)+PSFRange,pntNum);
            Source_Matrix(Xcol, Yrow) = Source_Matrix(Xcol,Yrow)+ exp(-( (X(Xcol,Yrow)-core_X(iMMA)).^2)/PSFSigma2).*exp(-((Y(Xcol,Yrow)-core_Y(iMMA)).^2)/PSFSigma2);
        end
    case 'sinc2'
        PSFRange = ceil(PSFSigma*length(Source_Coordinate_X)+1);
        for iMMA =1:MMA_Number
            Xcol = max(core_Xr(iMMA)-PSFRange,1):min(core_Xr(iMMA)+PSFRange,pntNum);
            Yrow = max(core_Yr(iMMA)-PSFRange,1):min(core_Yr(iMMA)+PSFRange,pntNum);
            Source_Matrix(Xcol, Yrow) = Source_Matrix(Xcol,Yrow)+ (sinc( hypot(X(Xcol,Yrow)-core_X(iMMA),Y(Xcol,Yrow)-core_Y(iMMA))/PSFSigma)).^2;
        end
    case 'sinc2rect'
        PSFRange = ceil(PSFSigma*length(Source_Coordinate_X)+1);
        for iMMA =1:MMA_Number
            Xcol = max(core_Xr(iMMA)-PSFRange,1):min(core_Xr(iMMA)+PSFRange,pntNum);
            Yrow = max(core_Yr(iMMA)-PSFRange,1):min(core_Yr(iMMA)+PSFRange,pntNum);
            Source_Matrix(Xcol, Yrow) = Source_Matrix(Xcol,Yrow)+ (sinc((X(Xcol,Yrow)-core_X(iMMA))/PSFSigma).*sinc((Y(Xcol,Yrow)-core_Y(iMMA))/PSFSigma)).^2;
        end
        
    case 'sincrect'
        PSFRange = ceil(PSFSigma*length(Source_Coordinate_X)*2+1);
        for iMMA =1:MMA_Number
            Xcol = max(core_Xr(iMMA)-PSFRange,1):min(core_Xr(iMMA)+PSFRange,pntNum);
            Yrow = max(core_Yr(iMMA)-PSFRange,1):min(core_Yr(iMMA)+PSFRange,pntNum);
            Source_Matrix(Xcol, Yrow) = Source_Matrix(Xcol,Yrow)+ abs(sinc((X(Xcol,Yrow)-core_X(iMMA))/PSFSigma).*sinc((Y(Xcol,Yrow)-core_Y(iMMA))/PSFSigma));
        end
    case 'bessel'
        
    otherwise
        error('LithoModel:Source:ParameterError','Unsupported PSF Shape');
end
Source_Matrix = PSFMAX * Source_Matrix;

%该部分可以加速
% 四分之一正向编码 获得光源序列编码map

switch lower(SourceConstructType)
    case 'full' %完全光源
        
    case 'quarter' % 四分之一解码恢复
        Source_Matrix = Source_Matrix + fliplr(Source_Matrix);
        Source_Matrix = Source_Matrix + flipud(Source_Matrix);
end



end

%% MMA卷积模型光源 连续MMA中心点 直角坐标系版本
function [ Source_Matrix, X, Y ] = CalculateSourceMatrixWithMMA_Cartesian(MMA_Number,MMA_X, MMA_Y,  Source_Coordinate_X, Source_Coordinate_Y,PSFType, PSFSigma,PSFMAX, SourceConstructType)

[ X, Y ] = meshgrid( Source_Coordinate_X, Source_Coordinate_Y );
Source_Matrix = zeros(size(X));

pntNum = length(Source_Coordinate_X);
core_Yr = round(MMA_X*(pntNum+1)/2+pntNum/2);
core_Xr = round(MMA_Y*(pntNum+1)/2+pntNum/2);

PSFRange = ceil(pntNum/20);
PSFSigma2 = PSFSigma.^2;
for iMMA =1:MMA_Number
    
    Xcol = max(core_Xr(iMMA)-PSFRange,1):min(core_Xr(iMMA)+PSFRange,pntNum);
    Yrow = max(core_Yr(iMMA)-PSFRange,1):min(core_Yr(iMMA)+PSFRange,pntNum);
    
    % 边界光源点单独处理
    if  core_Xr(iMMA) == (pntNum+1)/2 || core_Yr(iMMA) == (pntNum+1)/2
        Source_Matrix(Xcol, Yrow) = Source_Matrix(Xcol,Yrow)+ exp(-( (X(Xcol,Yrow)-MMA_X(iMMA)).^2+(Y(Xcol,Yrow)-MMA_Y(iMMA)).^2)/PSFSigma2) /2;
    else
        Source_Matrix(Xcol, Yrow) = Source_Matrix(Xcol,Yrow)+ exp(-( (X(Xcol,Yrow)-MMA_X(iMMA)).^2+(Y(Xcol,Yrow)-MMA_Y(iMMA)).^2)/PSFSigma2);
    end
    %                 %调试代码
    %                 xx = linspace(-1,1,101);
    %                 [ xx, yy ] = meshgrid( xx );
    %                 pcolor(xx,yy,Source_Matrix);hold
    % %                 scatter(core_X(iMMA),core_Y(iMMA));
    % %                 Source_Matrix(core_Xr(iMMA),core_Yr(iMMA)) = 1;
end
Source_Matrix = PSFMAX * Source_Matrix;

%该部分可以加速
% 四分之一正向编码 获得光源序列编码map

switch lower(SourceConstructType)
    case 'full' %完全光源
        
    case 'quarter' % 四分之一解码恢复
        Source_Matrix = Source_Matrix + fliplr(Source_Matrix);
        Source_Matrix = Source_Matrix + flipud(Source_Matrix);
end
end

%% 带高级参数的环形照明光源
% UNDERSTANDING ILLUMINATION EFFECTS FOR CONTROL OF OPTICAL PROXIMITY EFFECTS (OPE)
% Proc. of SPIE Vol. 6924, 69241U, (2008)  0277-786X
function [ Source_Matrix, X, Y ] = CalculateAnnularSourceMatrixWithAP(SigmaOut,SigmaIn,Source_Coordinate_X,Source_Coordinate_Y,AdvancedParameters)
ap = AdvancedParameters;
[X,Y]=meshgrid(Source_Coordinate_X,Source_Coordinate_Y);
R2 = X.^2 + Y.^2;
% Source_Matrix = AdvancedParameters.BG * ones(size(Source_Coordinate_X));
[theta, rho] = cart2pol(X-ap.x0,Y-ap.y0);
ri = SigmaIn + ap.dx0/2*cos(theta)+ap.dy0/2*sin(theta)+ap.cri2*cos(2*theta)+ap.sri2*sin(2*theta);
ro = SigmaOut - ap.dx0/2*cos(theta)-ap.dy0/2*sin(theta)+ap.cro2*cos(2*theta)+ap.sro2*sin(2*theta);
Iintensityazimuth = 1-erf(ap.slo/SigmaOut*(rho-ro)).*erf(ap.sli/SigmaIn*(rho-ri));
IntensityA = ones(size(X));
for iazm = 1:2
    IntensityA = IntensityA +(ap.cA(iazm)*cos(iazm*theta)+ap.sA(iazm)*sin(iazm*theta)).*(rho/SigmaOut).^iazm;
end
Source_Matrix = Iintensityazimuth .* IntensityA;
Source_Matrix = Source_Matrix/max(Source_Matrix(:));
Source_Matrix(Source_Matrix<ap.BG) = ap.BG;
Source_Matrix(R2>1) = 0;
end
%% 环形照明光源
function [Source_Matrix,X,Y] = CalculateAnnularSourceMatrix(SigmaOut, SigmaIn, Source_Coordinate_X, Source_Coordinate_Y)

[X,Y]=meshgrid(Source_Coordinate_X,Source_Coordinate_Y);
Radius=sqrt(X.^2+Y.^2);
Radius_Out=(Radius-SigmaOut)<=1e-10;
Radius_In=(Radius-SigmaIn)>=-1e-10;
Source_Matrix=(Radius_Out+Radius_In)-1;
end

%% 带高级参数的四级照明光源
% Predictive modeling of advanced illumination pupils used as imaging enhancement for low k1 applications
% Proceedings of SPIE Vol. 5377 (SPIE, Bellingham, WA, 2004) 0277-786X
% 极平衡需要进一步考虑
function [Source_Matrix,X,Y] = CalculateQuasarSourceMatrixWithAP(SigmaOut, SigmaIn, openAngle, rotateAngle, Source_Coordinate_X, Source_Coordinate_Y, AdvancedParameters)
ap = AdvancedParameters;
[X,Y]=meshgrid(Source_Coordinate_X,Source_Coordinate_Y);
R2 = X.^2+Y.^2;

% Source_Matrix = AdvancedParameters.BG * ones(size(Source_Coordinate_X));
[theta, rho] = cart2pol(X-ap.x0,Y-ap.y0);
newtheta = theta-rotateAngle;
ri = SigmaIn + ap.dx0/2*cos(newtheta)+ap.dy0/2*sin(newtheta)+ap.cri2*cos(2*newtheta)+ap.sri2*sin(2*newtheta);
ro = SigmaOut - ap.dx0/2*cos(newtheta)-ap.dy0/2*sin(newtheta)+ap.cro2*cos(2*newtheta)+ap.sro2*sin(2*newtheta);
Iintensityaxial = (1-erf(ap.slo/SigmaOut*(rho-ro)).*erf(ap.sli/SigmaIn*(rho-ri)))/2;

% 计算角向参数
newThetaFlag = newtheta<-pi;
newtheta(newThetaFlag) = newtheta(newThetaFlag)+2*pi;
newThetaFlag = newtheta>pi;
newtheta(newThetaFlag) = newtheta(newThetaFlag)-2*pi;
newtheta = pi/4 - abs(abs(abs(newtheta)-pi/2)-pi/4);
Iintensityazimuth = (1-erf(ap.sla*(newtheta+openAngle/2)).*erf(ap.sla*(newtheta-openAngle/2)))/2;
IntensityA = ones(size(X));
for iazm = 1:2
    IntensityA = IntensityA +(ap.cA(iazm)*cos(iazm*(theta-rotateAngle))+ap.sA(iazm)*sin(iazm*(theta-rotateAngle))).*(rho/SigmaOut).^iazm;
end
Source_Matrix = Iintensityaxial .* IntensityA.*Iintensityazimuth;
if ap.A>eps
    Source_Matrix = ap.A * Source_Matrix/max(Source_Matrix(:));
end
Source_Matrix(Source_Matrix<ap.BG) = ap.BG;
Source_Matrix(R2>1) = 0;

end
%% 四极光源
function [Source_Matrix,X,Y] = CalculateQuasarSourceMatrix(SigmaOut,SigmaIn,openAngle,Source_Coordinate_X,Source_Coordinate_Y)
[X,Y]=meshgrid(Source_Coordinate_X,Source_Coordinate_Y);
Source_Matrix=zeros(size(X));
Radius=sqrt(X.^2+Y.^2);
theta=atan2(Y,X);
Indextheta1=(abs(theta)<=openAngle/2 )   |  ( (abs(theta)>=pi-openAngle/2 ) & (abs(theta)<=pi ) )  ;
Indextheta2=(abs(theta-pi/2)<=openAngle/2 )   |    (abs(theta+pi/2)<=openAngle/2 );
Indextheta=Indextheta1 |  Indextheta2;
Index=(Radius<=SigmaOut) & (Radius>=SigmaIn) & Indextheta ;
Source_Matrix(Index)=1;
end

%% 带高级参数的二级照明光源
% Predictive modeling of advanced illumination pupils used as imaging enhancement for low k1 applications
% Proceedings of SPIE Vol. 5377 (SPIE, Bellingham, WA, 2004) 0277-786X
function [Source_Matrix,X,Y]=CalculateDipoleSourceMatrixWithAP(SigmaOut, SigmaIn, openAngle, rotateAngle, Source_Coordinate_X, Source_Coordinate_Y, AdvancedParameters)
ap = AdvancedParameters;
[X,Y]=meshgrid(Source_Coordinate_X,Source_Coordinate_Y);
R2 = X.^2+Y.^2;
% Source_Matrix = AdvancedParameters.BG * ones(size(Source_Coordinate_X));
[theta, rho] = cart2pol(X-ap.x0,Y-ap.y0);
newtheta = theta-rotateAngle;
ri = SigmaIn + ap.dx0/2*cos(newtheta)+ap.dy0/2*sin(newtheta)+ap.cri2*cos(2*newtheta)+ap.sri2*sin(2*newtheta);
ro = SigmaOut - ap.dx0/2*cos(newtheta)-ap.dy0/2*sin(newtheta)+ap.cro2*cos(2*newtheta)+ap.sro2*sin(2*newtheta);
Iintensityaxial = (1-erf(ap.slo/SigmaOut*(rho-ro)).*erf(ap.sli/SigmaIn*(rho-ri)))/2;
newtheta = theta-rotateAngle;
newThetaFlag = newtheta<-pi;
newtheta(newThetaFlag) = newtheta(newThetaFlag)+2*pi;
newThetaFlag = newtheta>pi;
newtheta(newThetaFlag) = newtheta(newThetaFlag)-2*pi;
newtheta = pi/2 - abs(abs(newtheta)-pi/2);
Iintensityazimuth = (1-erf(ap.sla*(newtheta+openAngle/2)).*erf(ap.sla*(newtheta-openAngle/2)))/2;
IntensityA = ones(size(X));
for iazm = 1:2
    IntensityA = IntensityA +(ap.cA(iazm)*cos(iazm*(theta-rotateAngle))+ap.sA(iazm)*sin(iazm*(theta-rotateAngle))).*(rho/SigmaOut).^iazm;
end
Source_Matrix = Iintensityaxial .* IntensityA.*Iintensityazimuth;
if ap.A~=0
    Source_Matrix = ap.A * Source_Matrix/max(Source_Matrix(:));
end
Source_Matrix(Source_Matrix<ap.BG) = ap.BG;
Source_Matrix(R2>1) = 0;

end

%% 二极光源
function  [Source_Matrix,X,Y]=CalculateDipoleSourceMatrix(SigmaOut,SigmaIn,openAngle, rotateAngle, Source_Coordinate_X,Source_Coordinate_Y)
[X,Y]=meshgrid(Source_Coordinate_X,Source_Coordinate_Y);
Source_Matrix=zeros(size(X));
Radius=sqrt(X.^2+Y.^2);
theta=atan2(Y,X);
%Indextheta=(abs(theta)<=angle/2 )   |  ( (abs(theta)>=pi-angle/2 ) & (abs(theta)<=pi ) )  ;
Indextheta=( abs(cos(theta-rotateAngle ))>=cos(openAngle/2) );
Index=(Radius<=SigmaOut) & (Radius>=SigmaIn) &Indextheta;
Source_Matrix(Index)=1;
end



%% 离散光源
function  [Source_Matrix,X,Y] = CalculateMultiCircSourceMatrix( SigmaCenter,SigmaRadius,PoleNumber,RotateAngle,Source_Coordinate_X,Source_Coordinate_Y )

[X,Y] = meshgrid(Source_Coordinate_X,Source_Coordinate_Y);
Source_Matrix = zeros(size(X));
rotateStep = 2*pi / PoleNumber;
for ii=1:PoleNumber
    
    [ xCenter, yCenter ] = pol2cart(RotateAngle+(ii-1)*rotateStep, SigmaCenter);
    Radius2 = (X-xCenter).^2+(Y-yCenter).^2;
    Index = Radius2 <= SigmaRadius^2;
    Source_Matrix(Index)=1;
end

end

%% 光源点数据的转换无论如何都要将作出的matrix形式的光源转化为source_data的形式。
function [Source_Data] = ConvertSourceMatrix2SourceData(Source_Matrix,X,Y)
rSquare = X.^2 +Y.^2;
sizeSourceX = size(Source_Matrix,1);
Source_Matrix(rSquare> 1) = 0;

Source_Data.X = reshape(X,sizeSourceX*sizeSourceX,1); %将X化为一列存储到Source_Data的第一列
Source_Data.Y = reshape(Y,sizeSourceX*sizeSourceX,1); %将Y化为一列存储到Source_Data的第二列
Source_Data.Value = reshape(Source_Matrix,sizeSourceX*sizeSourceX,1);%将Source_Matrix化为一列存储到Source_Data的第三列
end