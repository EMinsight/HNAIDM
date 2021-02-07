function aerailLithoImage = Calculate1DAerialImage(sr, mk, po, aerailLithoImage, rp, numerics)
% 优化版本
maskNf = numerics.SampleNumber_Mask_X;  % 81
waferNf = numerics.SampleNumber_Wafer_X; % 81

sr.PntNum = numerics.SampleNumber_Source;
[ sourceData, ~ ] =  sr.Calc_SourceSimple();

weight = sum( sourceData.Value );
wavelength = sr.Wavelength; % 193
NA = po.NA;
M = po.Reduction;
if strcmpi(po.LensType,'Immersion')
    indexImage = po.Index_ImmersionLiquid; % Default = 1.44
elseif strcmpi(po.LensType,'Dry')
     indexImage = 1;  
     if NA>=1
         error('Wrong NA!');
     end
else
    error('Unsupported Lens Type')
end

% 掩模旋转
if mk.Orientation == 0
    % 无变化
elseif abs( mk.Orientation - pi/2)<eps
    tenpY = sourceData.Y;
    sourceData.Y = -sourceData.X;
    sourceData.X = tenpY;
else
    error('Not supported orientation angle')
end

% 掩模频谱
mk.Nf = maskNf;
[ spectrum, mask_fs, ~ ] = mk.CalculateMaskSpectrum(po, sr);

%
SimulationRange = rp.FocusRange - rp.Focus ;

%% 开始计算
Intensity = zeros(length(SimulationRange),waferNf-1);%%%%%%%%

ExyzCalculateNumber = 1;

fm_s = repmat(sourceData.X,1,maskNf-1);
gm_s = repmat(sourceData.Y,1,maskNf-1);

fm_sm = fm_s + mask_fs(1,1:end-1); % NOTE:   %%需要MATLAB 2017以上版本
gm_sm = gm_s;  % 

rho2 = fm_sm.^2 + gm_sm.^2;

validPupil = (rho2<=1);     % 通过归一化的孔径来限定有效孔径
f_calc = fm_sm(validPupil); 
g_calc = gm_sm(validPupil);

[theta_calc, rho_calc] = cart2pol(f_calc,g_calc);
fgSquare = rho_calc.^2 ;

% 倾斜因子
obliquityFactor = sqrt(sqrt((1-(M^2*NA^2).*fgSquare)./(1-(NA/indexImage)^2.*fgSquare))); 
% 计算波像差
if mk.Orientation == 0
    % 无变化
    aberration = po.CalculateAberrationFast(rho_calc,theta_calc, 0);
elseif abs( mk.Orientation - pi/2)<eps
    aberration = po.CalculateAberrationFast(rho_calc,theta_calc, pi/2);%已经修改
else
    error('Not supported orientation angle')
end

SpectrumCalc = ones(length(sourceData.Value),1)*spectrum(1:end-1); % TODO: 需要优化
SpectrumCalc = SpectrumCalc(validPupil); % 有效孔径范围内的掩模频谱

TempH0Aber = SpectrumCalc.* obliquityFactor.*exp(1j*2*pi*aberration); 
% 掩模频谱、修正后的倾斜因子、像差，还缺少光瞳函数。（掩模、光源离轴照明、投影物镜像差）


obliqueRaysMatrix = ones(length(rho_calc),ExyzCalculateNumber); % 倾斜光线矩阵
if strcmp( numerics.ImageCalculationMode,'vector')
    ExyzCalculateNumber = 3;
    
    f_calc_s = fm_s(validPupil);
    g_calc_s = gm_s(validPupil);
    [ nts, nrs ] = cart2pol(f_calc_s,g_calc_s);
    
    %偏振分解
    [ PolarizedX, PolarizedY ] = sr.Calc_PolarizationMap( nts, nrs );
    
    % 计算偏振传输矩阵M
    [M0xx,M0yx,M0xy,M0yy,M0xz,M0yz] = CalculateCharacteristicMatrix(f_calc,g_calc,fgSquare,NA,indexImage);
        
    obliqueRaysMatrix = ones(length(rho_calc),ExyzCalculateNumber);  % 根据矢量模型来建立的倾斜光线矩阵
    obliqueRaysMatrix(:,1) = PolarizedX.*M0xx + PolarizedY.*M0yx;
    obliqueRaysMatrix(:,2) = PolarizedX.*M0xy + PolarizedY.*M0yy;
    obliqueRaysMatrix(:,3) = PolarizedX.*M0xz + PolarizedY.*M0yz;
    
end


if ~strcmp(po.PupilFilter.Type,'none') % 此时光瞳函数存在低通滤波特性，还可能是pupilFilterGauss
    filter = po.PupilFilter.Type;
    parameter = po.PupilFilter.Parameter;
    f_pupil = fm_sm(validPupil);
    g_pupil = gm_sm(validPupil);
    pupilFilterData = feval(filter,parameter,f_pupil,g_pupil); % 光瞳函数。有了这个变量，就为后面计算完整的频谱数据提供基础
    TempH0Aber = pupilFilterData.*TempH0Aber; % 完整的频谱计算
    
end

TempFocus = -1j*2*pi/wavelength.*sqrt(indexImage^2 - NA*NA.*fgSquare); % 这个变量是表示离焦量吗？

for iFocus = 1:length(SimulationRange)
    intensity = zeros(1,waferNf-1);

    
    tempF = exp(TempFocus.*SimulationRange(iFocus));
    for iEM = 1:ExyzCalculateNumber
        rho2(:) = 0;
        rho2(validPupil) = obliqueRaysMatrix(:,iEM).*TempH0Aber .*tempF;
        if waferNf == maskNf
            ExyzFrequency = rho2;
        elseif waferNf >maskNf
            ExyzFrequency = [zeros( length(sourceData.Value), (waferNf-maskNf)/2),rho2,zeros( length(sourceData.Value), (waferNf-maskNf)/2)];
        else
            ExyzFrequency = rho2;
            ExyzFrequency(:, 1:(maskNf-waferNf)/2) = [ ];
            ExyzFrequency(:, (end-(maskNf-waferNf)/2+1):end) = [ ];
        end

        Exyz = fft(ExyzFrequency,[],2); 
        IntensityCon = real(Exyz).^2 + imag(Exyz).^2;
        
        IntensityTemp = indexImage * (mask_fs(2) - mask_fs(1)).^2 * fftshift(sourceData.Value' * IntensityCon);
        intensity = intensity + IntensityTemp;

    end

    Intensity(iFocus,:) = intensity./weight;

end
ImageX = linspace(-mk.Period_X/2,mk.Period_X/2,waferNf);
ImageY = 0;
ImageZ = rp.FocusRange;

aerailLithoImage.ImageType = '1d';
aerailLithoImage.SimulationType = 'aerial';
aerailLithoImage.Intensity = fliplr([Intensity, Intensity(:,1)]);
aerailLithoImage.ImageX = ImageX;
aerailLithoImage.ImageY = ImageY;
aerailLithoImage.ImageZ = ImageZ;

end