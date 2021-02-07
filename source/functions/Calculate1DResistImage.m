function lithoImage = Calculate1DResistImage(source, mask, projector, filmStack, lithoImage, receipe, numerics)

maskNf = numerics.SampleNumber_Mask_X;  % Default 81
waferNf = numerics.SampleNumber_Wafer_X; % 81

source.PntNum = numerics.SampleNumber_Source;
[ sourceData, ~ ] =  source.Calc_SourceSimple();
weightSource = sum( sourceData.Value );
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

% ��ģ��ת
if mask.Orientation == 0
    % �ޱ仯
elseif abs( mask.Orientation - pi/2)<eps
    tenpY = sourceData.Y;
    sourceData.Y = -sourceData.X;
    sourceData.X = tenpY;
else
    error('Not supported orientation angle')
end

% ��ģƵ��
mask.Nf = maskNf;
[ spectrum, mask_fs, ~] = mask.CalculateMaskSpectrum(projector ,source);

%% ��ʼ����


if ~strcmp( numerics.ImageCalculationMode,'vector')
    error('Projector model must be vector');
end

if length(receipe.Focus)>1
    error('Not support multi focus');
end
fm_s = repmat(sourceData.X,1,maskNf-1);
gm_s = repmat(sourceData.Y,1,maskNf-1);

fm_sm = fm_s + mask_fs(1,1:end-1); % NOTE:   %%��ҪMATLAB 2017���ϰ汾
gm_sm = gm_s;  % ��%Slow

rho2 = fm_sm.^2 + gm_sm.^2;

validPupil = find(rho2<=1); % ͨ����һ���Ŀ׾����޶���Ч�׾�
f_calc = fm_sm(validPupil);
g_calc = gm_sm(validPupil);

[theta_calc, rho_calc] = cart2pol(f_calc,g_calc);
fgSquare = rho_calc.^2 ;

% ��б����
obliquityFactor = sqrt(sqrt((1-(M^2*NA^2).*fgSquare)./(1-(NA/indexImage)^2.*fgSquare))); 
% ��ǿ�Ƚ�����������б����
% ���㲨���
if mask.Orientation == 0
    % �ޱ仯
    aberration = projector.CalculateAberrationFast(rho_calc,theta_calc, 0);
elseif abs( mask.Orientation - pi/2)<eps
    aberration = projector.CalculateAberrationFast(rho_calc,theta_calc, pi/2);%�Ѿ��޸�
else
    error('Not supported orientation angle')
end

SpectrumCalc = ones(length(sourceData.Value),1)*spectrum(1:end-1); % TODO: ��Ҫ�Ż�
SpectrumCalc = SpectrumCalc(validPupil); % ��Ч�׾���Χ�ڵ���ģƵ��

TempH0Aber = SpectrumCalc.* obliquityFactor.*exp(1j*2*pi*aberration); %��ģƵ�ס����������б���ӡ�����ȱ�ٹ�ͫ����������ģ����Դ����������ͶӰ�ﾵ��

f_calc_s = fm_s(validPupil);
g_calc_s = gm_s(validPupil);
[ nts, nrs ] = cart2pol(f_calc_s,g_calc_s);

[ PolarizedX, PolarizedY ] = source.Calc_PolarizationMap( nts, nrs );

% NA/indexImage��Ӧ�����е�sin(theta_obj)����optical imaging in projection microlithography��ʽ��6.13��
alpha = (NA/indexImage)*f_calc;  % x��������
beta = (NA/indexImage)*g_calc;   % y��������
gamma = sqrt(1-(NA/indexImage)^2*fgSquare);%z��������
vectorRho2 = 1-gamma.^2; % ��ĸ��������ڲ�ʿ���ģ�ʽ��2.30���ͣ�2.31�������˾���ı������

%% ����P����
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

%���ĵ㵥������
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


%% �����ͫ
if ~strcmp(projector.PupilFilter.Type,'none') % ��ʱ��ͫ�������ڵ�ͨ�˲����ԣ���������pupilFilterGauss
    filter = projector.PupilFilter.Type;
    parameter = projector.PupilFilter.Parameter;
    f_pupil = fm_sm(validPupil);
    g_pupil = gm_sm(validPupil);
    pupilFilterData = feval(filter,parameter,f_pupil,g_pupil); % ��ͫ���������������������Ϊ�������������Ƶ�������ṩ����
    TempH0Aber = pupilFilterData.*TempH0Aber; % ������Ƶ�׼���
end

%% ����͸��ϵ���÷���ϵ��
%%%%%%%% TODO: �����ɺ��������� ���Ϊ���㴫����󼰹��ڸĲ��ڵ����������ֵ
thcknessResist = filmStack.GetResistThickness;
indexResist = filmStack.GetResistIndex();
indexSubstrate = filmStack.GetSubstrateIndex();
TARCs = filmStack.GetTARCLayers();
BARCs = filmStack.GetBARCLayers();
eta0 = sqrt(Numerics.Mu0 / Numerics.Epsilon0);

%%%%%     mu0 = 4*pi*10-7; %��մŵ���
%%%%%     epsilon0 = 8.854187817*10-12; %��ս�糣��
% NOTE: ��Դŵ��ʽ���Ϊ1

cosThetaInc = gamma;%z��������;
sinThetaInc = sqrt(1-gamma.^2);%
% �������ͬ�Խ��ʼ���
sinThetaResist = sinThetaInc * indexImage/indexResist;
cosThetaResist = sqrt(1-sinThetaResist.^2);
sinThetaSubstrate = sinThetaInc * indexImage/indexSubstrate;
cosThetaSubstrate = sqrt(1-sinThetaSubstrate.^2);

etaInc = eta0/indexImage;
etaResist = eta0/indexResist;
etaSubstrate = eta0/indexSubstrate;

% TEƫ��͸����ϵ���뷴��ϵ��
% TopARCs
faiInc = cosThetaInc/etaInc;
faiResist = cosThetaResist/etaResist;
[M11,M12,M21,M22] = CalculateTransferMatrix(TARCs,sinThetaInc,indexImage,wavelength,'TE');
rhoStacks = ( faiInc.*(M11-faiResist.*M12) + (M21-faiResist.*M22) ) ./ ( faiInc.*(M11-faiResist.*M12) - (M21-faiResist.*M22) );%������
tauStacks = 2*faiInc ./ ( faiInc.*(M11-faiResist.*M12) - (M21-faiResist.*M22) );%͸����

% BottomArcs
faiSubstrate = cosThetaSubstrate/etaSubstrate;
[M11,M12,M21,M22] = CalculateTransferMatrix(BARCs,sinThetaInc,indexImage,wavelength,'TE');
rhoSubStacks = ( faiResist.*(M11-faiSubstrate.*M12) + (M21-faiSubstrate.*M22) ) ./ ( faiResist.*(M11-faiSubstrate.*M12) - (M21-faiSubstrate.*M22) );
% tauSubStacks = 2*faiResist ./ ( faiResist.*(M11-faiSubstrate.*M12) - (M21-faiSubstrate.*M22) );

% TMƫ��͸����ϵ���뷴��ϵ��
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

%% ���㽺����

tempF = exp(-1j*2*pi/wavelength .* sqrt( indexImage^2 - NA*NA.*fgSquare).* (-receipe.Focus) );
TempHAber = TempH0Aber .* tempF;
if isempty(numerics.SimulationRange_Resist)
     ImageZ_Resist = linspace(0,thcknessResist,numerics.SampleNumber_Wafer_Z);
    Sample_RIZ = numerics.SampleNumber_Wafer_Z;
elseif min(numerics.SimulationRange_Resist)<0 || max(numerics.SimulationRange_Resist)>thcknessResist 
    error('Wrong focus range');
else
    Sample_RIZ = length(numerics.SimulationRange_Resist);
    ImageZ_Resist = numerics.SimulationRange_Resist;
end
%     ImageZ_Resist = linspace(0,thcknessResist,numerics.SampleNumber_Wafer_Z);
IntensityResist = zeros(Sample_RIZ,maskNf-1);%%%%%%%%
expWaveD = exp(2*1j*2*pi/wavelength*indexResist.*cosThetaResist.*thcknessResist);
rho2(:) = 0;
for idepth = 1:Sample_RIZ

    waveZ = 1j*2*pi/wavelength*indexResist.*cosThetaResist.*(thcknessResist-ImageZ_Resist(idepth));

    expPosWaveZ = exp(waveZ);
    expNegWaveDZ = expWaveD./expPosWaveZ;
    %     %%%%%%%% Optical Imaging in ProjectionMicrolithography Chapter 6 & 8
    %     ȱ��ע��
    %     Fs = tauStacks./tauSubStacks.*( exp(waveZ) + rhoSubStacks.*exp(-waveZ));
    %     Fpxy = tauStackp./tauSubStackp.*( exp(waveZ) - rhoSubStackp.*exp(-waveZ));
    %     Fpz = tauStackp./tauSubStackp.*( exp(waveZ) + rhoSubStackp.*exp(-waveZ));
    
    %%%% Inside Prolith
    %%%% Chris A. Mack, Analytical expression for the standing wave
    %%%% intensity in photoresist, APPLIED OPTICS
    %%%% / Vol. 25, No. 12 / 15 June 1986 
    Fs = tauStacks./(1 + rhoStacks.*rhoSubStacks.*expWaveD).*( expPosWaveZ + rhoSubStacks.*expNegWaveDZ);
    Fpxy = tauStackp./(1 + rhoStackp.*rhoSubStackp.*expWaveD).*( expPosWaveZ - rhoSubStackp.*expNegWaveDZ);
    Fpz = tauStackp./(1 + rhoStackp.*rhoSubStackp.*expWaveD).*( expPosWaveZ + rhoSubStackp.*expNegWaveDZ);
    
    %% ������������
    MSxx = Fs.*Pxsx + Fpxy.*Pxpx;
    MSyx = Fs.*Pysx + Fpxy.*Pypx;
    MSxy = Fs.*Pxsy + Fpxy.*Pxpy;
    MSyy = Fs.*Pysy + Fpxy.*Pypy;
    MSxz = Fpz.*Pxpz;
    MSyz = Fpz.*Pypz;
    
    obliqueRaysMatrixX = PolarizedX.*MSxx + PolarizedY.*MSyx;
    obliqueRaysMatrixY = PolarizedX.*MSxy + PolarizedY.*MSyy;
    obliqueRaysMatrixZ = PolarizedX.*MSxz + PolarizedY.*MSyz;
    
    rho2(validPupil) = obliqueRaysMatrixX.*TempHAber;
    Ex = fft(rho2,[],2);
    rho2(validPupil) = obliqueRaysMatrixY.*TempHAber;
    Ey = fft(rho2,[],2);
    rho2(validPupil) = obliqueRaysMatrixZ.*TempHAber;
    Ez = fft(rho2,[],2);

    IntensityCon = real(Ex).^2 + imag(Ex).^2 + real(Ey).^2 + imag(Ey).^2+real(Ez).^2 + imag(Ez).^2;
    
    IntensityTemp = real(indexResist) * (mask_fs(2) - mask_fs(1)).^2 /weightSource * fftshift(sourceData.Value' * IntensityCon,2);

    IntensityResist(idepth,:) = IntensityTemp;
end
ImageX = linspace(-mask.Period_X/2,mask.Period_X/2,waferNf);
ImageY = 0;
ImageZ = ImageZ_Resist;

calcNf = maskNf;

if waferNf == calcNf
    intensityOutPut =  [IntensityResist, IntensityResist(:,1)];
elseif waferNf > calcNf
    IntensityFrequency = fftshift(fft(IntensityResist,[],2),2);
    IntensityFrequency = [IntensityFrequency,IntensityFrequency(:,1)];
    IntensityFrequencyWafer = [zeros( size(IntensityResist,1), (waferNf-calcNf)/2),IntensityFrequency,zeros( size(IntensityResist,1), (waferNf-calcNf)/2-1)];
    IntensityFrequencyWafer = fftshift(IntensityFrequencyWafer,2);
    IntensityWafer = abs(ifft(IntensityFrequencyWafer,[],2)) * (waferNf-1)/(calcNf-1);
    intensityOutPut =  [IntensityWafer, IntensityWafer(:,1)];
else 
    error('WaferNf must be more than calcNf !');
end

lithoImage.ImageType = '1d';
lithoImage.SimulationType = 'resist';
lithoImage.Intensity = intensityOutPut;
lithoImage.ImageX = ImageX;
lithoImage.ImageY = ImageY;
lithoImage.ImageZ = ImageZ;

end

