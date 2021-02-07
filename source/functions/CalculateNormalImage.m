function normalIntensity = CalculateNormalImage(source, mask, projector, receipe, numerics)
% �Ż��汾�ο���Ϊ

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
normalized_Frequency = projector.NA/source.Wavelength; % ���и���Ҷ�任֮ǰ��ϵ����ʾ
switch lower(mask.MaskType)
    case '1d'
        SpectrumCalc = normalized_Frequency*mask.Period_X;
        normalized_Period_X = mask.Period_X * normalized_Frequency; % ��һ���������
        dfmdg = 1/normalized_Period_X;
    case '1dpixel'
        error('To be implement');
        %  Hk = PixelMask1D2Spectrum(mk,NA,Wavelength);
    case '2d'
        SpectrumCalc = normalized_Frequency^2 * mask.Period_X*mask.Period_Y;
        normalized_Period_X = mask.Period_X * normalized_Frequency; % ��һ���������
        normalized_Period_Y = mask.Period_Y * normalized_Frequency; % ��һ���������
        dfmdg = 1/normalized_Period_X/normalized_Period_Y;
    case '2dpixel'
        SpectrumCalc = normalized_Frequency^2 * mask.Period_X*mask.Period_Y;
        normalized_Period_X = mask.Period_X * normalized_Frequency; % ��һ���������
        normalized_Period_Y = mask.Period_Y * normalized_Frequency; % ��һ���������
        dfmdg = 1/normalized_Period_X/normalized_Period_Y;
    case '3d'
        % TODO: ��ά��ģ�㷨�����
end

%% ��ʼ����

ExyzCalculateNumber = 1;

f0_s = sourceData.X;
g0_s = sourceData.Y;

fgSquare = f0_s.^2 + g0_s.^2;

[theta_calc, rho_calc] = cart2pol(f0_s,g0_s);

% ��б����
obliquityFactor = sqrt(sqrt((1-(M^2*NA^2).*fgSquare)./(1-(NA/indexImage)^2.*fgSquare))); % ��ǿ�Ƚ�����������б����
% ���㲨���
if mask.Orientation == 0% �ޱ仯
    aberration = projector.CalculateAberrationFast(rho_calc,theta_calc, 0);
elseif abs( mask.Orientation - pi/2)<eps
    aberration = projector.CalculateAberrationFast(rho_calc,theta_calc, pi/2);%�Ѿ��޸�
else
    error('Not supported orientation angle')
end


TempH0Aber = SpectrumCalc.* obliquityFactor.*exp(1j*2*pi*aberration); 
%��ģƵ�ס����������б���ӡ�����ȱ�ٹ�ͫ����������ģ����Դ����������ͶӰ�ﾵ��

% ����ƫ�������M
obliqueRaysMatrix = ones(length(rho_calc),ExyzCalculateNumber); % ��б���߾���
if strcmp( numerics.ImageCalculationMode,'vector')
    ExyzCalculateNumber = 3;
    
    [ nts, nrs ] = cart2pol(f0_s,g0_s);
    
    [ PolarizedX, PolarizedY ] = source.Calc_PolarizationMap( nts, nrs );
    
    % NA/indexImage��Ӧ�����е�sin(theta_obj)����optical imaging in projection microlithography��ʽ��6.13��
    alpha = (NA/indexImage)*f0_s;  % x��������
    beta = (NA/indexImage)*g0_s;   % y��������
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
    
    %���ĵ㵥������
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
    
    obliqueRaysMatrix = ones(length(rho_calc),ExyzCalculateNumber);  % ����ʸ��ģ������������б���߾���
    obliqueRaysMatrix(:,1) = PolarizedX.*M0xx + PolarizedY.*M0yx;
    obliqueRaysMatrix(:,2) = PolarizedX.*M0xy + PolarizedY.*M0yy;
    obliqueRaysMatrix(:,3) = PolarizedX.*M0xz + PolarizedY.*M0yz;
    
end

if ~strcmp(projector.PupilFilter.Type,'none') % ��ʱ��ͫ�������ڵ�ͨ�˲����ԣ���������pupilFilterGauss
    filter = projector.PupilFilter.Type;
    parameter = projector.PupilFilter.Parameter;
    f_pupil = f0_s(validPupil);
    g_pupil = g0_s(validPupil);
    % ��ͫ���������������������Ϊ�������������Ƶ�������ṩ����
    pupilFilterData = feval(filter,parameter,f_pupil,g_pupil); 
    TempH0Aber = pupilFilterData.*TempH0Aber; % ������Ƶ�׼���
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