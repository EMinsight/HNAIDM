function [TCCMatrix_Stacked,FG_ValidSize] = Calculate2DTCCMatrix(source, mask, projector, receipe,numerics)

%% 数据初始化
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

wavelength = source.Wavelength;

xPitch = mask.Period_X;
yPitch = mask.Period_Y;

% 计算偏移的光瞳函数
if strcmpi( numerics.ImageCalculationMode,'scalar')
    [ SP, F_Valid, G_Valid, sourceData ] = CalculateShiftedPupilS(wavelength,projector,source,xPitch,yPitch,indexImage, receipe.Focus);
    TCCMatrix_Stacked = GetTCCMatrix(sourceData,SP);
elseif strcmpi( numerics.ImageCalculationMode,'vector')
    [ SPXX, SPXY, SPYX,SPYY, SPXZ, SPYZ, F_Valid, G_Valid, sourceData ] = CalculateShiftedPupilV(wavelength,projector,source,xPitch,yPitch,indexImage, receipe.Focus);
    TCCMatrixX = GetTCCMatrix(sourceData,SPXX+SPYX);
    TCCMatrixY = GetTCCMatrix(sourceData,SPXY+SPYY);
    TCCMatrixZ = GetTCCMatrix(sourceData,SPXZ+SPYZ);
    TCCMatrix_Stacked = TCCMatrixX+TCCMatrixY+TCCMatrixZ;
end
TCCMatrix_Stacked = projector.Index_ImmersionLiquid*TCCMatrix_Stacked;
FG_ValidSize = [length(G_Valid), length(F_Valid)];
end


function [shiftedPupil, f,g,sourceData] = CalculateShiftedPupilS(wavelength, projector,source,xPitch,yPitch,indexImage,focus)
[ sourceData, ~ ] = source.Calc_SourceSimple();
M = projector.Reduction;
NA = projector.NA;
normalized_xPitch = xPitch/(wavelength/NA);
normalized_yPitch = yPitch/(wavelength/NA);
Nf = ceil(2*normalized_xPitch); 
Ng = ceil(2*normalized_yPitch);
f = (1/normalized_xPitch) * (-Nf:1:Nf);
g = (1/normalized_yPitch) * (-Ng:1:Ng);
[ff,gg] = meshgrid(f,g);
new_f = reshape(ff,[],1) + reshape(sourceData.X,1,[]);
new_g = reshape(gg,[],1) + reshape(sourceData.Y,1,[]);
[ theta, rho ] = cart2pol(new_f,new_g);

validPupil = (rho<=1);
validRho = rho(validPupil);
validTheta = theta(validPupil);
validRhoSquare = validRho.^2;

obliquityFactor = sqrt(sqrt((1-(M^2*projector.NA^2).*validRhoSquare)./(1-((projector.NA/indexImage)^2).*validRhoSquare )));
Orientation = 0;
aberration = projector.CalculateAberrationFast(validRho,validTheta,Orientation);
% shiftedPupil = validPupil.*obliquityFactor.*exp(1j*2*pi*aberration);
shiftedPupil = zeros(size(validPupil));
TempFocus = 1j*2*pi/wavelength.*sqrt(indexImage^2-NA*NA.*validRhoSquare);
shiftedPupil(validPupil) = obliquityFactor .* exp(1j*2*pi*aberration) .* exp(TempFocus.*focus);

end

function [SPXX, SPXY, SPYX,SPYY, SPXZ, SPYZ, f,g,sourceData] = CalculateShiftedPupilV(wavelength, projector,source,xPitch,yPitch,indexImage,focus)
[ sourceData, ~ ] = source.Calc_SourceSimple();
M = projector.Reduction;
NA = projector.NA;
normalized_xPitch = xPitch/(wavelength/NA);
normalized_yPitch = yPitch/(wavelength/NA);
Nf = ceil(2*normalized_xPitch); 
Ng = ceil(2*normalized_yPitch);
f = (1/normalized_xPitch) * (-Nf:1:Nf);
g = (1/normalized_yPitch) * (-Ng:1:Ng);
[ff,gg] = meshgrid(f,g);

new_f = reshape(ff,[],1) + reshape(sourceData.X,1,[]);
new_g = reshape(gg,[],1) + reshape(sourceData.Y,1,[]);
[ theta, rho ] = cart2pol(new_f,new_g);
rhoSquare = rho.^2;
validPupil = (rho<=1);
validRho = rho(validPupil);
validTheta = theta(validPupil);
validRhoSquare = rhoSquare(validPupil);

obliquityFactor = sqrt(sqrt((1-(M^2*projector.NA^2).*validRhoSquare)./(1-((projector.NA/indexImage)^2).*validRhoSquare )));
Orientation = 0;
aberration = projector.CalculateAberrationFast(validRho,validTheta,Orientation);

shiftedPupil = zeros(size(validPupil));
TempFocus = 1j*2*pi/wavelength.*sqrt(indexImage^2-NA*NA.*validRhoSquare);
shiftedPupil(validPupil) = obliquityFactor .* exp(1j*2*pi*aberration) .* exp(TempFocus.*focus);

[M0xx,M0yx,M0xy,M0yy,M0xz,M0yz] = CalculateCharacteristicMatrix(new_f,new_g,rhoSquare,NA,indexImage);
% 提取光源偏振
[ theta_s, rho_s ] = cart2pol(sourceData.X,  sourceData.Y);
[ PolarizedX, PolarizedY ] = source.Calc_PolarizationMap( theta_s, rho_s );
PolarizedX = reshape(PolarizedX,1,[]);
PolarizedY = reshape(PolarizedY,1,[]);

SPXX = PolarizedX.*M0xx .* shiftedPupil;
SPXY = PolarizedX.*M0xy .* shiftedPupil;
SPXZ = PolarizedX.*M0xz .* shiftedPupil;

SPYX = PolarizedY.*M0yx .* shiftedPupil;
SPYY = PolarizedY.*M0yy .* shiftedPupil;
SPYZ = PolarizedY.*M0yz .* shiftedPupil;
end

function TCCMatrix = GetTCCMatrix(sourceData,shiftedPupil)
% J 指有效光源
% shiftedPupil 偏移的光瞳函数
i = 1:length(sourceData.Value);
j = i;
S = sparse(i,j,sourceData.Value);
TCCMatrix = shiftedPupil * S *shiftedPupil' ; % 利用矩阵的共轭转置得到HSH*
TCCMatrix = TCCMatrix./sum(sourceData.Value);% 体现积分作用，光源强度归一化到1
end
