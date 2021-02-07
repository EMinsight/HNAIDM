function normalIntensity = CalculateNormalImage(source, mask, projector, receipe, numerics)
% 优化版本参考面为

source.PntNum = numerics.SampleNumber_Source;
[ sourceData, ~ ] =  source.Calc_SourceSimple();
weight = sum( sourceData.Value );
wavelength = source.Wavelength; % 193
NA = projector.NA;
M = projector.Reduction;
if strcmpi(projector.LensType,'Immersion')
    indexImage = projector.Index_ImmersionLiquid;
elseif strcmpi(projector.LensType,'Dry')
     indexImage = 1;   % 1.44
     if NA>=1
         error('Wrong NA!');
     end
else
    error('Unsupported Lens Type')
end

% 掩模旋转
if mask.Orientation == 0
    % 无变化
elseif abs( mask.Orientation - pi/2)<eps
    tenpY = sourceData.Y;
    sourceData.Y = -sourceData.X;
    sourceData.X = tenpY;
else
    error('Not supported orientation angle')
end

% 掩模频谱
normalized_Frequency = projector.NA/source.Wavelength; % 进行傅里叶变换之前的系数表示
switch lower(mask.MaskType)
    case '1d'
        SpectrumCalc = normalized_Frequency*mask.Period_X;
        normalized_Period_X = mask.Period_X * normalized_Frequency; % 归一化采样间距
        dfmdg = 1/normalized_Period_X;
    case '1dpixel'
        error('To be implement');
        %  Hk = PixelMask1D2Spectrum(mk,NA,Wavelength);
    case '2d'
        SpectrumCalc = normalized_Frequency^2 * mask.Period_X*mask.Period_Y;
        normalized_Period_X = mask.Period_X * normalized_Frequency; % 归一化采样间距
        normalized_Period_Y = mask.Period_Y * normalized_Frequency; % 归一化采样间距
        dfmdg = 1/normalized_Period_X/normalized_Period_Y;
    case '2dpixel'
        SpectrumCalc = normalized_Frequency^2 * mask.Period_X*mask.Period_Y;
        normalized_Period_X = mask.Period_X * normalized_Frequency; % 归一化采样间距
        normalized_Period_Y = mask.Period_Y * normalized_Frequency; % 归一化采样间距
        dfmdg = 1/normalized_Period_X/normalized_Period_Y;
    case '3d'
        % TODO: 三维掩模算法待完成
end

%% 开始计算

ExyzCalculateNumber = 1;

f0_s = sourceData.X;
g0_s = sourceData.Y;

fgSquare = f0_s.^2 + g0_s.^2;

[theta_calc, rho_calc] = cart2pol(f0_s,g0_s);

% 倾斜因子
obliquityFactor = sqrt(sqrt((1-(M^2*NA^2).*fgSquare)./(1-(NA/indexImage)^2.*fgSquare))); % 对强度进行修正的倾斜因子
% 计算波像差
if mask.Orientation == 0% 无变化
    aberration = projector.CalculateAberrationFast(rho_calc,theta_calc, 0);
elseif abs( mask.Orientation - pi/2)<eps
    aberration = projector.CalculateAberrationFast(rho_calc,theta_calc, pi/2);%已经修改
else
    error('Not supported orientation angle')
end


TempH0Aber = SpectrumCalc.* obliquityFactor.*exp(1j*2*pi*aberration); 
%掩模频谱、修正后的倾斜因子、像差，还缺少光瞳函数。（掩模、光源离轴照明、投影物镜像差）

% 计算偏振传输矩阵M
obliqueRaysMatrix = ones(length(rho_calc),ExyzCalculateNumber); % 倾斜光线矩阵
if strcmp( numerics.ImageCalculationMode,'vector')
    ExyzCalculateNumber = 3;
    
    [ nts, nrs ] = cart2pol(f0_s,g0_s);
    
    [ PolarizedX, PolarizedY ] = source.Calc_PolarizationMap( nts, nrs );
    
    % NA/indexImage对应文献中的sin(theta_obj)，《optical imaging in projection microlithography》式（6.13）
    alpha = (NA/indexImage)*f0_s;  % x方向余弦
    beta = (NA/indexImage)*g0_s;   % y方向余弦
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
    
    %中心点单独处理
    centerpoint = find(fgSquare<eps);
    if length(centerpoint)>eps
        Pxsx(centerpoint) = 1;
        Pysx(centerpoint) = 0;
        Pxsy(centerpoint) = 0;
        Pysy(centerpoint) = 0;
        
        Pxpx(centerpoint) = 0;
        Pypx(centerpoint) = 0;
        Pxpy(centerpoint) = 0;
        Pypy(centerpoint) = 1;
        Pxpz(centerpoint) = 0;
        Pypz(centerpoint) = 0;
    end
    
    M0xx = Pxsx + Pxpx;
    M0yx = Pysx + Pypx;
    M0xy = Pxsy + Pxpy;
    M0yy = Pysy + Pypy;
    M0xz = Pxpz;
    M0yz = Pypz;
    
    obliqueRaysMatrix = ones(length(rho_calc),ExyzCalculateNumber);  % 根据矢量模型来建立的倾斜光线矩阵
    obliqueRaysMatrix(:,1) = PolarizedX.*M0xx + PolarizedY.*M0yx;
    obliqueRaysMatrix(:,2) = PolarizedX.*M0xy + PolarizedY.*M0yy;
    obliqueRaysMatrix(:,3) = PolarizedX.*M0xz + PolarizedY.*M0yz;
    
end

if ~strcmp(projector.PupilFilter.Type,'none') % 此时光瞳函数存在低通滤波特性，还可能是pupilFilterGauss
    filter = projector.PupilFilter.Type;
    parameter = projector.PupilFilter.Parameter;
    f_pupil = f0_s(validPupil);
    g_pupil = g0_s(validPupil);
    % 光瞳函数。有了这个变量，就为后面计算完整的频谱数据提供基础
    pupilFilterData = feval(filter,parameter,f_pupil,g_pupil); 
    TempH0Aber = pupilFilterData.*TempH0Aber; % 完整的频谱计算
end

TempFocus = -1j*2*pi/wavelength.*sqrt(indexImage^2 - NA*NA.*fgSquare); 

tempF = exp(TempFocus.*receipe.Focus);
intensityBlank = 0;
for iEM = 1:ExyzCalculateNumber
    ExyzFrequency = obliqueRaysMatrix(:,iEM).*TempH0Aber .*tempF;    
    Exyz = fft(ExyzFrequency,[],2);
    IntensityCon = real(Exyz).^2 + imag(Exyz).^2;
    IntensityTemp = indexImage * dfmdg.^2*fftshift(sourceData.Value' * IntensityCon);
    intensityBlank = intensityBlank + IntensityTemp;
end

normalIntensity = intensityBlank / weight;
end