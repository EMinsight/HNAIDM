function resistLithoImage = Calculate2DResistImageV3(source, mask, projector, filmStack, resistLithoImage, receipe, numerics)

%% 数据初始化
maskNf = numerics.SampleNumber_Mask_X;  % 默认81
maskNg = numerics.SampleNumber_Mask_Y;  % 默认81
waferNf = numerics.SampleNumber_Wafer_X; % 默认81
waferNg = numerics.SampleNumber_Wafer_Y; % 默认81
waferNz = numerics.SampleNumber_Wafer_Z; %默认101



source.PntNum = numerics.SampleNumber_Source;
[ sourceData, ~ ] =  source.Calc_SourceSimple();

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

% 掩模频谱
mask.Nf = maskNf;
mask.Ng = maskNg;
[ spectrumMask, mask_fs, mask_gs ] = mask.CalculateMaskSpectrum(projector ,source);


%计算实际需要的网格数
index_calc_fs = abs(mask_fs)<2;
index_calc_gs = abs(mask_gs)<2;
mask_calc_fs = mask_fs(index_calc_fs);
mask_calc_gs = mask_gs(index_calc_gs);
spectrumCalc = spectrumMask(index_calc_gs,index_calc_fs);
calcNf = length(mask_calc_fs);
calcNg = length(mask_calc_gs);


%% 开始计算

Orientation = mask.Orientation;

if ~strcmp( numerics.ImageCalculationMode,'vector')
    error('Projector model must be vector');
end

if length(receipe.Focus)>1
    error('Not support multi focus');
end
% if length(receipe.FocusRange)>1
%     error('Not support multi focus');
% end

[ mask_fm, mask_gm ] = meshgrid( mask_calc_fs(1:end-1), mask_calc_gs(1:end-1) );

sourceX = sourceData.X;
sourceY = sourceData.Y;
sourceV = sourceData.Value;
weight = sum( sourceV );

[ source_theta, source_rho] = cart2pol(sourceX,  sourceY);
[ PolarizedX, PolarizedY ] = source.Calc_PolarizationMap( source_theta, source_rho );
new_spectrum = spectrumCalc(1:end-1,1:end-1);
mask_fg2m = mask_fm.^2 + mask_gm.^2;
sourceXY2 =  sourceX.^2+ sourceY.^2;

%初始数据
thcknessResist = filmStack.GetResistThickness;
indexResist = filmStack.GetResistIndex();
indexSubstrate = filmStack.GetSubstrateIndex();
TARCs = filmStack.GetTARCLayers();
BARCs = filmStack.GetBARCLayers();
eta0 = sqrt(Numerics.Mu0 / Numerics.Epsilon0);


if isempty(numerics.SimulationRange_Resist)
    ImageZ_Resist = linspace(0,thcknessResist,waferNz);
else
    ImageZ_Resist = numerics.SimulationRange_Resist;
end
Exyzp = single(zeros(calcNg-1,calcNf-1,length(ImageZ_Resist)));
intensity3D = zeros(calcNg-1,calcNf-1,length(ImageZ_Resist));
for iSource = 1:length(sourceV)  %%parfor
    
    rho2 = mask_fg2m +2*(sourceX(iSource)*mask_fm + sourceY(iSource)*mask_gm) + sourceXY2(iSource);
    
    validPupil = find(rho2<=1);
    f_calc = mask_fm(validPupil) + sourceX(iSource);
    g_calc = mask_gm(validPupil) + sourceY(iSource);
    
    [theta_calc, rho_calc] = cart2pol(f_calc,g_calc);
    fgSquare = rho_calc.^2 ;
    %计算方向余弦
    alpha = (NA/indexImage)*f_calc;
    beta = (NA/indexImage)*g_calc;
    gamma = sqrt(1-(NA/indexImage)^2*fgSquare);
    
    obliquityFactor = sqrt(sqrt((1-(M^2*NA^2).*fgSquare)./ (1-((NA/indexImage)^2).*fgSquare)));%倾斜因子
    aberration = projector.CalculateAberrationFast( rho_calc,theta_calc, Orientation);%波像差
    pupilFilter = projector.CalculatePupilFilter( rho_calc, theta_calc );%光瞳函数
    tempFocus = exp(-1j*2*pi/wavelength .* sqrt( indexImage^2 - NA*NA.*fgSquare).*-receipe.Focus);%离焦
    SpectrumCalc = new_spectrum(validPupil);
    
    TempHAber = SpectrumCalc.*obliquityFactor.*exp(1j*2*pi.*aberration).*pupilFilter.* tempFocus;%传递函数
    
    %% 计算偏振矩阵
    vectorRho2 = alpha.^2+beta.^2;
    Pxsx = (beta.^2)./vectorRho2;  % 此处修改运算符号 20171108
    Pysx = (-alpha.*beta)./vectorRho2;
    Pxsy = Pysx;
    Pysy = (alpha.^2)./vectorRho2;
    
    Pxpx = gamma.*Pysy;
    Pypx = -gamma.*Pysx;
    Pxpy = Pypx;
    Pypy = gamma.*Pxsx;
    Pxpz = -alpha;
    Pypz = -beta;
    
    %中心点单独处理
    centerpoint = find(vectorRho2<eps);
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
    
    %% 计算透射系数和反射系数
    cosThetaInc = gamma;%z方向余弦;
    sinThetaInc = sqrt(1-gamma.^2);%
    % 满足各项同性介质假设
    sinThetaResist = sinThetaInc * indexImage/indexResist;
    cosThetaResist = sqrt(1-sinThetaResist.^2);
    sinThetaSubstrate = sinThetaInc * indexImage/indexSubstrate;
    cosThetaSubstrate = sqrt(1-sinThetaSubstrate.^2);
    
    etaInc = eta0/indexImage;
    etaResist = eta0/indexResist;
    etaSubstrate = eta0/indexSubstrate;
    
    % TE偏振透射率系数与反射系数
    % TopARCs
    faiInc = cosThetaInc/etaInc;
    faiResist = cosThetaResist/etaResist;
    [M11,M12,M21,M22] = CalculateTransferMatrix(TARCs,sinThetaInc,indexImage,wavelength,'TE');
    rhoStacks = ( faiInc.*(M11-faiResist.*M12) + (M21-faiResist.*M22) ) ./ ( faiInc.*(M11-faiResist.*M12) - (M21-faiResist.*M22) );%反射率
    tauStacks = 2*faiInc ./ ( faiInc.*(M11-faiResist.*M12) - (M21-faiResist.*M22) );%透射率
    
    % BottomArcs
    faiSubstrate = cosThetaSubstrate/etaSubstrate;
    [M11,M12,M21,M22] = CalculateTransferMatrix(BARCs,sinThetaInc,indexImage,wavelength,'TE');
    rhoSubStacks = ( faiResist.*(M11-faiSubstrate.*M12) + (M21-faiSubstrate.*M22) ) ./ ( faiResist.*(M11-faiSubstrate.*M12) - (M21-faiSubstrate.*M22) );
    % tauSubStacks = 2*faiResist ./ ( faiResist.*(M11-faiSubstrate.*M12) - (M21-faiSubstrate.*M22) );
    
    % TM偏振透射率系数与反射系数
    % TopARCs
    faiIncTM = -cosThetaInc*etaInc;
    faiResistTM = -cosThetaResist*etaResist;
    [M11,M12,M21,M22] = CalculateTransferMatrix(TARCs,sinThetaInc,indexImage,wavelength,'TM');
    rhoStackp = ( faiIncTM.*(M11-faiResistTM.*M12) + (M21-faiResistTM.*M22) ) ./ ( faiIncTM.*(M11-faiResistTM.*M12) - (M21-faiResistTM.*M22) );
    tauStackp = -2*etaResist.*cosThetaInc ./ ( faiIncTM.*(M11-faiResistTM.*M12) - (M21-faiResistTM.*M22) );
    
    % BottomARCs
    faiSubstrateTM = -cosThetaSubstrate*etaSubstrate;
    [M11,M12,M21,M22] = CalculateTransferMatrix(BARCs,sinThetaInc,indexImage,wavelength,'TM');
    
    rhoSubStackp = ( faiResistTM.*(M11-faiSubstrateTM.*M12) + (M21-faiSubstrateTM.*M22) ) ./ ( faiResistTM.*(M11-faiSubstrateTM.*M12) - (M21-faiSubstrateTM.*M22) );
    % tauSubStackp = -2*etaSubstrate*cosThetaResist ./ ( faiResistTM.*(M11-faiSubstrateTM.*M12) - (M21-faiSubstrateTM.*M22) );
    
    %% 计算胶内像
    expWaveD2 = exp(2*1j*2*pi/wavelength*indexResist.*cosThetaResist.*thcknessResist);
    A2FsH = tauStacks./(1 + rhoStacks.*rhoSubStacks.*expWaveD2) .*TempHAber;
    A2FpxyH = tauStackp./(1 + rhoStackp.*rhoSubStackp.*expWaveD2) .*TempHAber;
    %     A2FpzH = tauStackp./(1 + rhoStackp.*rhoSubStackp.*expWaveD2) .*TempHAber;
    A2FpzH = A2FpxyH;
    
    expPosWaveZ = exp(1j*2*pi/wavelength*indexResist.*cosThetaResist.*(thcknessResist-ImageZ_Resist));
    expNegWaveD2Z = expWaveD2./expPosWaveZ;
    
    Fs = A2FsH.*( expPosWaveZ + rhoSubStacks.*expNegWaveD2Z);
    Fpxy = A2FpxyH.*( expPosWaveZ - rhoSubStackp.*expNegWaveD2Z);
    Fpz = A2FpzH.*( expPosWaveZ + rhoSubStackp.*expNegWaveD2Z);
    
    MSxx = Fs.*Pxsx + Fpxy.*Pxpx;
    MSyx = Fs.*Pysx + Fpxy.*Pypx;
    MSxy = Fs.*Pxsy + Fpxy.*Pxpy;
    MSyy = Fs.*Pysy + Fpxy.*Pypy;
    MSxz = Fpz.*Pxpz;
    MSyz = Fpz.*Pypz;
    
    ExValid = PolarizedX(iSource).*MSxx + PolarizedY(iSource).*MSyx;
    EyValid = PolarizedX(iSource).*MSxy + PolarizedY(iSource).*MSyy;
    EzValid = PolarizedX(iSource).*MSxz + PolarizedY(iSource).*MSyz;
    
    zIndex = 0:((calcNf-1)*(calcNg-1)):(length(ImageZ_Resist)-1)*((calcNf-1)*(calcNg-1));%构建索引
    validPupilIndexFull = reshape( validPupil + zIndex,[],1);
    
    Ex = Exyzp;
    Ey = Exyzp;
    Ez = Exyzp;
    
    Ex(validPupilIndexFull) = reshape( ExValid,[],1);
    Ey(validPupilIndexFull) = reshape( EyValid,[],1);
    Ez(validPupilIndexFull) = reshape( EzValid,[],1);
    
    Ex = fft2(Ex);
    Ey = fft2(Ey);
    Ez = fft2(Ez);
    intensityTemp = real(Ex).^2 + imag(Ex).^2 + real(Ey).^2 + imag(Ey).^2+real(Ez).^2 + imag(Ez).^2;
    intensity3D = intensity3D + intensityTemp*sourceV(iSource);
end

dfmdg = (mask_fs(2)-mask_fs(1))*(mask_gs(2)-mask_gs(1));
intensity3D = dfmdg^2*fftshift(fftshift(intensity3D,1),2);

intensity3D = real(indexResist)*intensity3D/weight;

if waferNf == calcNf && waferNg == calcNg
    intensityOutPut = intensity3D;

elseif waferNf < calcNf || waferNg < calcNg
    error('Wafer grid must be larger than calculate grid !');
else
    intensityOutPut = zeros(waferNg-1,waferNf-1,length(ImageZ_Resist));
    
    intensity3DFrequency = fftshift(fftshift(fft2(intensity3D),1),2);
    intensity3DFrequency(:,calcNf,:) = intensity3DFrequency(:,1,:);
    intensity3DFrequency(calcNg,:,:) = intensity3DFrequency(1,:,:);
    
    rangeNg = (waferNg+1)/2 -(calcNg-1)/2 :(waferNg+1)/2 + (calcNg-1)/2;
    rangeNf = (waferNf+1)/2 -(calcNf-1)/2 :(waferNf+1)/2 + (calcNf-1)/2;
    
    intensityOutPut(rangeNg,rangeNf,:) = intensity3DFrequency;
    intensityOutPut = fftshift(fftshift(intensityOutPut,1),2);
    intensityOutPut = abs(ifft2(intensityOutPut)) * (waferNf-1)/(calcNf-1) * (waferNg-1)/(calcNg-1);
    
end
    intensityOutPut(:,waferNf,:) =  intensityOutPut(:,1,:);
    intensityOutPut(waferNg,:,:) =  intensityOutPut(1,:,:);
% 不能用FFT插值
% intensity3DWafer = interpft(intensity3D,waferNf,2);
% intensity3DWafer = interpft(intensity3DWafer,waferNg,1);

resistLithoImage.ImageType = '2d';
resistLithoImage.SimulationType = 'resist';
resistLithoImage.Intensity = intensityOutPut;
resistLithoImage.ImageX = linspace(-mask.Period_X/2,mask.Period_X/2,waferNf);
resistLithoImage.ImageY = linspace(-mask.Period_Y/2,mask.Period_Y/2,waferNg);
resistLithoImage.ImageZ = ImageZ_Resist;

end


