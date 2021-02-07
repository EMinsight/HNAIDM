classdef Source
    properties
        PntNum = 101;
        Wavelength = 193.368;
        Shape = 'annular';  
        %光源形状  'annular'  'multipole' 'dipolecirc'  'quadrupole' 'MMA' 'pixel' 'quasar'
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
        ModulationType % Null / ASML / Nikon 该参数设置照明光瞳调制方案
        %% Nikon 照明光瞳调制方案参数
        ZIM
        ZDM
        %% ASML照明光瞳调制方案参数
        ModulationPara
        %% PSF参数
        PSFEnable
        PSFSigma
        %% 光源偏振
        PolarizationType 
        %光源偏振模式  'x_pol'   'y_pol'  'unpolarized'   'r_pol'  't_pol' 'fun' 
        %TODO 'line' 
        PolarizationParameters
        
        %% MMA描述模型参数
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
    
    %% 构造函数 光源生成函数
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
            s.ModulationType = 'null';
            
            %Nikon 照明光瞳调制方案初始化
            s.ZIM = zeros(1,37);
            s.ZDM = zeros(1,55);
            % ASML 照明光瞳调制方案初始化
            s.ModulationPara.A = 1;
            s.ModulationPara.cA = zeros(1,4);
            s.ModulationPara.sA = zeros(1,4);
            s.ModulationPara.apo_2 = 0;
            s.ModulationPara.apo_10 = 0;
            s.ModulationPara.smear_x = 1;
            s.ModulationPara.smear_y = 1;
            s.ModulationPara.BGshort  = 0; %TODO: 需要进一步分析实现算法
            s.ModulationPara.BGlong = 0;   %TODO: 需要进一步分析实现算法
            s.ModulationPara.BGdark = 0;   %TODO: 需要进一步分析实现算法 目前仅有该BG参数可用
            s.ModulationPara.scal_r = 1;
            s.ModulationPara.offset_r = 0;
            s.ModulationPara.x0 = 0;
            s.ModulationPara.y0 = 0;
            s.ModulationPara.cr = zeros(1,3);
            s.ModulationPara.sr = zeros(1,3);
            s.ModulationPara.phi0 = 0;
            s.ModulationPara.cphi = zeros(1,3);
            s.ModulationPara.sphi = zeros(1,3);
            
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
            
            %             % MMA 模块初始化
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
            if  strcmpi(s.ModulationType,'null')
                [ Source_Data, Source_Matrix,X,Y] = s.Calc_SourceAll();
            elseif strcmpi(s.ModulationType,'nikon')
                [ Source_Data, Source_Matrix,X,Y] = s.ModulatedSourcWithNikonModel();
            elseif strcmpi(s.ModulationType,'asml')
                [ Source_Data, Source_Matrix,X,Y] = s.ModulatedSourcWithASMLModel();
            else
                error('Unsupported Modulation Type!');
            end
            Low_Weight = Source_Data.Value<1e-5;%去除光强小于某值的点
            % r2 = Source_Data(:,1).^2 + Source_Data(:,2).^2;
            % delRaw = (r2>1) | Low_Weight;
            Source_Data.X(Low_Weight)=[];        
            Source_Data.Y(Low_Weight)=[];        
            Source_Data.Value(Low_Weight)=[]; 
        end
        
        function [Source_Data,Source_Matrix,X,Y] = Calc_SourceValid(s) %照明光瞳范围内的序列化光源点
            if  strcmpi(s.ModulationType,'null')
                [Source_Data,Source_Matrix,X,Y] = s.Calc_SourceAll();
            elseif strcmpi(s.ModulationType,'nikon')
                [Source_Data,Source_Matrix,X,Y] = s.ModulatedSourcWithNikonModel();
            elseif strcmpi(s.ModulationType,'asml')
                [Source_Data,Source_Matrix,X,Y] = s.ModulatedSourcWithASMLModel();
            else
                error('Unsupported Modulation Type!');
            end
            
            Radius = Source_Data.X.^2 + Source_Data.Y.^2;
            Radius_Valid = find(Radius<=1);
            
            Source_Data.X=Source_Data.X(Radius_Valid);%
            Source_Data.Y=Source_Data.Y(Radius_Valid);%
            Source_Data.Value=Source_Data.Value(Radius_Valid);%
        end
        
        function [Source_Data,Source_Matrix,X,Y] = Calc_SourceAll(s) %完整的序列化光源点

            Source_Coordinate_X = linspace(-1,1,s.PntNum);
            
            %% 其他照明形式
            switch lower(s.Shape) %变成小写
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

                case 'point' %兼容设置
                    % warning('Please use "pixel" to replace the "point" source shape');
                    [sizeX,sizeY] = size(s.SPointValue);
                    if ~isequal(sizeX, s.PntNum)
                        error('Source Matrix Size X and difined PntNum are not matched ');
                    end
                    if ~isequal(sizeY, s.PntNum)
                        error('Source Matrix Size Y and difined PntNum are not matched ');
                    end
                    Source_Data.X = reshape(s.SPointX,[],1);%第一列存储x坐标
                    Source_Data.Y = reshape(s.SPointY,[],1);%第二列存储y坐标
                    Source_Data.Value = reshape(s.SPointValue,[],1);%第三列存储光强值
                    X = s.SPointX;
                    Y = s.SPointY;
                    Source_Matrix = s.SPointValue; 
                case 'mma' % MMA模型自由照明
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
            %             Source_Data(Low_Weight,3)=0;%【】后面的[]表示将该项除去
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

    %% 光源调制函数
    methods
        %% 将光源按照指定方式和参数调制
        
        %% ASML 调制方案
        % Generation of arbitrary freeform source shapes using advanced illumination systems in high-NA immersion scanners
        % Proc. of SPIE Vol. 7640, 764005
        function [ Source_Data, Source_Matrix, X, Y ] = ModulatedSourcWithASMLModel(s)
            mp = s.ModulationPara;
            [ ~, Source_Matrix, X, Y ] = Calc_SourceAll(s);
            
            r_High = 1;
            r_Low = 0;
            rXY = hypot( X, Y );
            
            % 二分法计算 r_shift
            for i=1:16
                r_shift = (r_Low+r_High)/2;
                Source_In = Source_Matrix(rXY<r_shift);
                if  sum(Source_In(:))<sum(Source_Matrix(:))/2
                    r_Low = r_shift;
                else
                    r_High = r_shift;
                end
            end
            
            % 分布调制
            [theta, rho] = cart2pol(X-mp.x0, Y-mp.y0);
            
            RFTerm = zeros(size(rho)); % 径向傅里叶分解项
            for k = 2:3
                RFTerm = RFTerm + mp.cphi(k).*cos(k*theta) + mp.sphi(k).*sin(k*theta);
            end
            rhoNew = (mp.offset_r+r_shift+(rho-r_shift)*mp.scal_r).*(1+RFTerm);
            
            TFTerm = zeros(size(theta));
            for k = 1:3
                TFTerm = TFTerm + mp.cr(k).*cos(k*theta) + mp.sr(k).*sin(k*theta);
            end
            thetaNew = theta + mp.phi0 + TFTerm;
            
            % 光强调制
            [thetaRaw, rhoRaw] = cart2pol(X, Y);
            
            AFTerm = zeros(size(thetaRaw));% 角向傅里叶分解项
            for k = 1:4
                AFTerm = AFTerm + mp.cA(k)*rhoRaw^k.*cos(k*thetaRaw) + mp.sA(k)*rhoRaw^k.*sin(k*thetaRaw);
            end
            IntensityA = mp.A * ( 1 + AFTerm + mp.apo_2 * rhoRaw^2 + mp.apo_10 * rhoRaw^10 );
            
            [X,Y] = pol2cart(thetaNew,rhoNew);
            
            Source_Matrix = Source_Matrix .*  IntensityA;
            
            %增加背景光
            Source_Matrix( Source_Matrix<mp.BGdark ) = mp.BGdark;
            
            %序列化
            Source_Data = ConvertSourceMatrix2SourceData(Source_Matrix,X,Y);
        end
        
        %% Nikon 调制方案
        % Tolerancing analysis of customized illumination for practical applications of source and mask optimization
        % Proc. of SPIE Vol. 7640, 764007
        function [ Source_Data, Source_Matrix, X, Y ] = ModulatedSourcWithNikonModel(s)
            [ ~, Source_Matrix, X, Y ] = Calc_SourceAll(s);
            [theta, rho] = cart2pol(X, Y);
            ZernikePoly = CalculateZernikePolynomial(rho, theta, length(s.ZIM));
            
            % 强度调制
            Source_Intesnity_Factor = zeros(size(X));
            for iOrder = 1:37
                Source_Intesnity_Factor = Source_Intesnity_Factor + s.ZIM(iOrder).* ZernikePoly{iOrder};
            end
            Source_Matrix = Source_Matrix.* exp(Source_Intesnity_Factor);
            maxValue = max(Source_Matrix(:));
            if maxValue>1
                s.ZIM(1) = log(1/maxValue);
                Source_Matrix = Source_Matrix/maxValue;
            end
            
            % 形变调制
            [X,Y] = CalculateModulatedCoordinate(s.ZDM,X,Y, ZernikePoly);
            
            %序列化
            Source_Data=ConvertSourceMatrix2SourceData(Source_Matrix,X,Y);
        end
        
    end
    %% 光源相关功能
    methods % MAA
        %% 基于模拟退火算法的MAA
        function [sr, convC] = MAA_Anneal(sr, TargetSourceMatrix,options_MMA,options_SA)
            %输出为优化后光源和适应度值随评价函数调用次数的关系
            if ~strcmpi(sr.Shape,'MMA')
                error('LithoModel:Source:PropertiesError','Source Shape must be MMA');
            end
            if nargin <2
                    error('LithoModel:Source:InputErroe','No Target Source Map!');
            elseif nargin<3
                options_MMA = struct('RemodInvalidMMA',0);
            end
            if nargin<4
                options_SA = struct(...
                    'CoolSched',@(T) (.98*T),...
                    'Generator',[],...
                    'InitTemp',0.0001,...
                    'MaxConsRej',1000,...
                    'MaxSuccess',20,...
                    'MaxTries',100,...
                    'StopTemp',1e-9,...
                    'StopVal',0,...
                    'Verbosity',0);
            end
            
%             xx = linspace(-1,1,sr.PntNum);
%             [ xx, yy ] = meshgrid( xx );
            
            srTemp = sr;
            srTemp.MMA_Number = 1;
            srTemp.MMA_Coordinate = 'polar';
            srTemp.MMA_Rho = 0;
            srTemp.MMA_Theta = 0;
            srTemp.MMA_SourceConstructType = 'full';
            [~,kernel,~,~] = srTemp.Calc_SourceAll;
            kernelsum = sum(kernel(:))/srTemp.MMA_PSFMAX;
            
            rho = sqrt(rand(1,sr.MMA_Number));
            b = rand(1,sr.MMA_Number);
            
            switch  lower(sr.MMA_SourceConstructType)
                case 'full' %完全光源
                    PSFMAX = sum( TargetSourceMatrix(:) )/ ( sr.MMA_Number.*kernelsum );
                    theta_Range = [0,pi*2];
                    theta = pi*b*2;
                case 'quarter' % 四分之一解码恢复
                    PSFMAX = sum( TargetSourceMatrix(:) )/ ( 4*sr.MMA_Number.*kernelsum );
                    theta_Range = [0,pi/2];
                    theta = pi*b/2;
            end
            
            rho_Range = [0,1];
            
            sr.MMA_PSFMAX = PSFMAX;

            xInit = [rho,theta];
            
            options_SA.Generator = @(x) Source.GenerateNewXSA(x,sr.MMA_Number,rho_Range,theta_Range);
            
            [xBest, ~,convC] = Optimizer.Anneal(@(xxx)Source.fitnessSAMMA(xxx,sr,TargetSourceMatrix),xInit,options_SA);
            if ~isfield(options_MMA,'RemodInvalidMMA')
                options_MMA.RemodInvalidMMA = 0;
            end
            
            if options_MMA.RemodInvalidMMA
                
            end
            
            sr.MMA_Rho = xBest(1:sr.MMA_Number);
            sr.MMA_Theta = xBest(sr.MMA_Number+1:end);

        end
        
        %% 基于随机翻转的MAA
        function sr = MAA_RandomFlip(sr, TargetSourceMatrix, options_MMA, options_RF)
            %输出为优化后光源
            if ~strcmpi(sr.Shape,'MMA')
                error('LithoModel:Source:PropertiesError','Source Shape must be MMA');
            end
            if nargin <2
                error('LithoModel:Source:InputErroe','No Target Source Map!');
            elseif nargin<3
                options_MMA = struct('RemodInvalidMMA',1);
            end
            if nargin<4
                options_RF = struct('MaxTries',100);
            end
            
            if ~isfield(options_RF,'MaxTries')
                options_RF.MaxTries = 100;
            end
            
            
            %%初始化参数
            pntNum = sr.PntNum;
            
            
%             xx = linspace(-1,1,sr.PntNum);
%             [ xx, yy ] = meshgrid( xx );
            %             kernel = exp(-( xx.^2 + yy.^2) / sr.MMA_PSFSigma.^2);
            %             kernelsum = sum(kernel(:));
            srTemp = sr;
            srTemp.MMA_Number = 1;
            srTemp.MMA_Coordinate = 'polar';
            srTemp.MMA_Rho = 0;
            srTemp.MMA_Theta = 0;
            srTemp.MMA_SourceConstructType = 'full';
            [~,kernel,~,~] = srTemp.Calc_SourceAll;
            kernelsum = sum(kernel(:))/srTemp.MMA_PSFMAX;
            switch  lower(sr.MMA_SourceConstructType)
                case 'full' %完全光源
                    PSFMAX = sum( TargetSourceMatrix(:) )/ ( sr.MMA_Number.*kernelsum );
                case 'quarter' % 四分之一解码恢复
                    PSFMAX = sum( TargetSourceMatrix(:) )/ ( 4*sr.MMA_Number.*kernelsum );
            end
            
            sr.MMA_PSFMAX = PSFMAX;
            
            %   调试代码
            %   srView = sr;
            
            FreeMMA_Number = sr.MMA_Number;
            MMA_State = false(1,sr.MMA_Number);
            %% 迭代计算
            rho = zeros(1,sr.MMA_Number);
            theta = zeros(1,sr.MMA_Number);
            
            for iter = 1:options_RF.MaxTries
                if FreeMMA_Number<1
                    break;
                end
                startMMA = min(sr.MMA_Number - FreeMMA_Number+1,sr.MMA_Number);
                a = rand(1,FreeMMA_Number);
                b = rand(1,FreeMMA_Number);
                rho(startMMA:end) = sqrt(a);
                theta(startMMA:end) = pi*b/2;
                
                
                sr.MMA_Rho = rho;
                sr.MMA_Theta = theta;

                [~,TempSourceMatrix,~,~] = sr.Calc_SourceAll();
                [ core_X, core_Y ] = pol2cart( theta, rho );
                core_Yr = round(core_X*(pntNum+1)/2+pntNum/2);
                core_Xr = round(core_Y*(pntNum+1)/2+pntNum/2);
                
                for iMMA = startMMA:sr.MMA_Number
                    if TempSourceMatrix(core_Xr(iMMA),core_Yr(iMMA)) ...
                            < TargetSourceMatrix(core_Xr(iMMA),core_Yr(iMMA))-max(kernel(:))/kernelsum
                        MMA_State(iMMA) = true;
                    end
                end
                
                %           调试代码
                %                 theta2 = theta(MMA_State);
                %                 rho2 = rho(MMA_State);
                %                 MMA_Number2 = sum(MMA_State);
                %                 srView.MMA_Number = MMA_Number2;
                %                 srView.MMA_Rho = rho2;
                %                 srView.MMA_Theta = theta2;
                %                 figure(200);
                %                 LithoView.PlotSource(srView);
                %                 snapnow;
                
                Lockedindex = MMA_State == true;
                LockedMMA_Number = sum(Lockedindex(:));
                FreeMMA_Number = sr.MMA_Number - LockedMMA_Number;
                rho(1:LockedMMA_Number) = rho (Lockedindex);
                theta(1:LockedMMA_Number) = theta(Lockedindex);
                MMA_State = sort(MMA_State,'descend');
            end
            %             sr.MMA_Number = LockedMMA_Number;
            %             sr.MMA_Rho = rho (1:LockedMMA_Number);
            %             sr.MMA_Theta = theta(1:LockedMMA_Number);
            sr.MMA_Rho = rho;
            sr.MMA_Theta = theta;
        end
    end
    
    
    methods (Static)
        %% MAA解卷积算法 维纳滤波 像素化光源
        function SourceMatrix_deconved = deconv_Pixel( TargetSourceMatrix,PSFSigma )
            [ pntNum, ~ ] =  size(TargetSourceMatrix);

            %% 解卷积
            xx = linspace(-1,1,pntNum);
            [ xx, yy ] = meshgrid( xx );
            kernel = 1/sqrt(2*pi)/PSFSigma *exp(-( xx.^2 + yy.^2) / PSFSigma.^2);
            kernel = kernel/sum(kernel(:));
            SourceMatrix_deconved = deconvwnr (TargetSourceMatrix, kernel);
            SourceMatrix_deconved(SourceMatrix_deconved>1) = 1;
            SourceMatrix_deconved(SourceMatrix_deconved<0) = 0;
        end
        
        %% MAA解卷积算法 维纳滤波 MMA光源
        function [ MMA_Theta, MMA_Rho,MMA_Number, PSFMAX] = deconv_MAA( TargetSourceMatrix, MMA_Number,PSFSigma, SourceConstructType )
            %  SourceConstructType = 'quarter'; %%%%%'quarter' or 'full'
            [ pntNum, ~ ] =  size(TargetSourceMatrix);

            %% 解卷积
            xx = linspace(-1,1,pntNum);
            [ xx, yy ] = meshgrid( xx );
            
            
            kernel = exp(-( xx.^2 + yy.^2) / PSFSigma.^2);
            kernelsum = sum(kernel(:));
            switch  lower(SourceConstructType)
                case 'full' %完全光源
                    PSFMAX = sum( TargetSourceMatrix(:) )/ ( MMA_Number.*kernelsum );
                case 'quarter' % 四分之一解码恢复
                    PSFMAX = sum( TargetSourceMatrix(:) )/ ( 4*MMA_Number.*kernelsum );
            end
            
            kernel = PSFMAX * kernel;
            
            deconvMatrix = deconvwnr(TargetSourceMatrix, kernel, 0);
            roundSourceMatrixCenter = round(abs(deconvMatrix));
            
            switch lower(SourceConstructType)
                case 'full' %完全光源
                    
                    SourceMatrixCenter = roundSourceMatrixCenter;
                case 'quarter' % 四分之一解码恢复
                    xFlag = xx>=0;
                    yFlag = yy>=0;
                    xyFlag = xFlag & yFlag;
                    SourceMatrixCenterRaw = abs(deconvMatrix);
                    SourceMatrixCenter = SourceMatrixCenterRaw .* xyFlag;
                    if mod(pntNum,2) == 1
                        SourceMatrixCenter((pntNum+1)/2,:) =  SourceMatrixCenter((pntNum+1)/2,:)/2;
                        SourceMatrixCenter(:,(pntNum+1)/2) =  SourceMatrixCenter(:,(pntNum+1)/2)/2;
                    end
                    SourceMatrixCenter = round(SourceMatrixCenter);
                otherwise
            end
            
            [ cX, cY, Value ] = find( SourceMatrixCenter);
            
            TotalNumber = sum(SourceMatrixCenter(:));
            centerX = zeros(1, TotalNumber);
            centerY = centerX;
            
            ic = 1;
            for iPoint = 1 : length( cX )
                for iSubPoint = 1 : Value( iPoint )
                    centerX(ic) = xx( cX(iPoint), cY(iPoint) );
                    centerY(ic) = yy( cX(iPoint), cY(iPoint) );
                    ic = ic + 1;
                end
            end
            MMA_Number = ic-1;
            [ MMA_Theta, MMA_Rho ] = cart2pol( centerX, centerY );
            PSFMAX = sum( TargetSourceMatrix(:) ) /( kernelsum * sum(roundSourceMatrixCenter(:)) );
        end
        
        %% SA MAA 相关函数
        function f = fitnessSAMMA(xx,srv,refSM)
            D = length(xx);
            srv.MMA_Rho = xx(1:D/2);
            srv.MMA_Theta = xx(D/2+1:end);
            [~,TempSourceMatrix,~,~] = srv.Calc_SourceAll();
            f = rms(TempSourceMatrix(:)-refSM(:));
            if max(TempSourceMatrix(:)>1)
                f = f*2;
            end
        end
        
        function newX = GenerateNewXSA(oldX,MMA_Number,rho_Range,theta_Range)
            newX = oldX;
            index = find(randperm(MMA_Number*2) == MMA_Number);
            if index <= MMA_Number
                deltaRho = (rho_Range(2)-rho_Range(1));
                newValue = oldX(index) + rand * deltaRho;
                if newValue>rho_Range(2)
                    newValue = newValue - deltaRho;
                elseif newValue<rho_Range(1)
                    newValue = newValue + deltaRho;
                end
                newX(index) = newValue;
            else
                deltaTheta = (theta_Range(2)-theta_Range(1));
                newValue = oldX(index) + rand *deltaTheta;
                if newValue>theta_Range(2)
                    newValue = newValue - deltaTheta;
                elseif newValue<theta_Range(1)
                    newValue = newValue + deltaTheta;
                end
                newX(index) = newValue;
            end
        end
    end
end

% 光源函数
%% MMA卷积模型光源 离散MMA中心点
function [ Source_Matrix, X, Y ] = CalculateSourceMatrixWithMMA(MMA_Number,MMA_Rho, MMA_Theta,  Source_Coordinate_X, Source_Coordinate_Y,PSFType, PSFSigma,PSFMAX, SourceConstructType)

[X,Y]=meshgrid(Source_Coordinate_X,Source_Coordinate_Y);
Source_Partial = zeros(size(X));
[ core_X, core_Y ] = pol2cart(MMA_Theta, MMA_Rho);
pntNum = length(Source_Coordinate_X);
core_Yr = round((core_X+1)*pntNum/2+0.5);
core_Xr = round((core_Y+1)*pntNum/2+0.5);

% 该部分可以加速
% core_XY = [core_Xr;core_Yr];
% Source_Matrix(core_XY) = ones(2,MMA_Number);


for iMMA =1:MMA_Number
    Source_Partial(core_Xr(iMMA),core_Yr(iMMA)) = Source_Partial(core_Xr(iMMA),core_Yr(iMMA)) + 1;
    
    % 调试代码
    %                 xx = linspace(-1,1,101);
    %                 [ xx, yy ] = meshgrid( xx );
    %                 pcolor(xx,yy,Source_Partial);hold
    %                 scatter(core_X(iMMA),core_Y(iMMA));
    %     Source_Matrix(core_Xr(i),core_Yr(i)) = 1;
end

%该部分可以加速
% 四分之一正向编码 获得光源序列编码编码map

switch lower(SourceConstructType)
    case 'full' %完全光源
        
    case 'quarter' % 四分之一解码恢复
        xF = X>=0;
        yF = Y>=0;
        IndexMap = xF&yF;
        Weight = IndexMap;
        Source_Full  = Source_Partial;
        Source_Full = Source_Full + fliplr(Source_Full);
        Source_Full = Source_Full + flipud(Source_Full);
        
        Weight = Weight + fliplr(Weight);
        Weight = Weight + flipud(Weight);
        Weight = abs(Weight-0.5)+0.5;
        Source_Full =Source_Full./Weight;
end

switch lower(PSFType)
    case 'gaussian'
        %                     s.PntNum = 501; % 测试参数
        %                     s.PSFSigma = 0.02;
        kernelSize = round(pntNum/10)*2+1;
        kernelEdge = 1/(pntNum-1)*(kernelSize-1);
        [kernelX,KernelY] = meshgrid(linspace(-kernelEdge,kernelEdge,kernelSize));
        kernel = exp(-(kernelX.^2+KernelY.^2)/PSFSigma.^2);%归一化
        %                     kernel = exp(-(kernelX.^2+KernelY.^2)/PSFSigma.^2);% 非归一化
        kernel(:,all(kernel<1e-6,1)) = [];
        kernel(all(kernel<1e-6,2),:) = [];
        
        Source_Matrix = PSFMAX*conv2(Source_Full,kernel,'same');
    case 'sinc'
    case 'none'
        Source_Matrix = Source_Full;
    otherwise
end


end

%% MMA卷积模型光源 连续MMA中心点 速度慢
function [ Source_Matrix, X, Y ] = CalculateSourceMatrixWithMMA_V2(MMA_Number,MMA_Rho, MMA_Theta, Source_Coordinate_X, Source_Coordinate_Y,PSFType, PSFSigma,PSFMAX, SourceConstructType)


[X,Y]=meshgrid(Source_Coordinate_X,Source_Coordinate_Y);

Source_Matrix = zeros(length(Source_Coordinate_X));
[ core_X, core_Y ] = pol2cart(MMA_Theta, MMA_Rho);

pntNum = length(Source_Coordinate_X);
xx = linspace(-1,1,pntNum);
[ xx, yy ] = meshgrid( xx );
subx = (pntNum+1)/2-5:pntNum;
suby = (pntNum+1)/2-5:pntNum;
switch lower(PSFType)
    case 'gaussian'
        for iMMA =1 : MMA_Number
            %     temp = 1/sqrt(2*pi)/PSFSigma *exp(-( (xx-core_X(iMMA)).^2+(yy-core_Y(iMMA)).^2)/PSFSigma.^2);
            temp = exp(-( (xx(subx,suby)-core_X(iMMA)).^2+(yy(subx,suby)-core_Y(iMMA)).^2)/PSFSigma.^2);
            Source_Matrix(subx,suby) = Source_Matrix(subx,suby) + temp;
        end
    otherwise
end

switch lower(SourceConstructType)
    case 'full'
    case 'quarter'
        Source_Matrix = Source_Matrix + fliplr(Source_Matrix);
        Source_Matrix = PSFMAX*(Source_Matrix + flipud(Source_Matrix));
end


end

%% MMA卷积模型光源 连续MMA中心点 速度较快
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

%% Nikon照明光瞳形变泽尼克像生成函数
function [X, Y] = CalculateModulatedCoordinate(ZDM,X,Y,ZernikePoly)

Factor  = [1,0,2,2,3,4,0,5,6,5,6,7,7,8,10,11,10,11,9,0,12,13,12,13,17,18,17,18,14,14,15,19,20,19,20,26,27,16,0,21,22,21,22,28,29,23,23,24,30,31,25,0,32,33,34; ...
    0,1,3,-3,2,0,4,6,-5,-6,5,8,-8,7,11,-10,-11,10,0,9,13,-12,-13,12,18,-17,-18,17,15,-15,14,20,-19,-20,19,27,-26,0,16,22,-21,-22,21,29,-28,24,-24,23,31,-30,0,25,33,-32,35];
%             [theta, rho] = cart2pol(X, Y);
dx = zeros(size(X));
dy = zeros(size(Y));
for i = 1:size(Factor,2)
    if Factor(1,i)~=0
        dx = dx + ZDM(i)*ZernikePoly{Factor(1,i)};
    end
    if Factor(2,i)~=0
        dy = dy + ZDM(i)*sign(Factor(2,i))*ZernikePoly{abs(Factor(2,i))};
    end
end

X = X + dx;
Y = Y + dy;

end