function TCCMatrix_Kernel = DecomposeTCC_SOCS(TCCMatrix_Stacked,FG_ValidSize, numerics)

maskNf = numerics.SampleNumber_Mask_X;  
maskNg = numerics.SampleNumber_Mask_Y; 
waferNf = numerics.SampleNumber_Wafer_X; 
waferNg = numerics.SampleNumber_Wafer_Y; 

% 得到TCC像
% 利用奇异值分解(SVD)
[ U, S, ~ ] = svd(TCCMatrix_Stacked);
% socsNumber
if strcmpi(numerics.Hopkins_SettingType,'order')
    socsNumber = min(numerics.Hopkins_Order,size(U,2)-1);
elseif strcmpi(numerics.Hopkins_SettingType,'threshold')
    rateThreshold = numerics.Hopkins_Threshold; % 阈值设置
    singularVector = diag(S); % 奇异向量
    totalSingular = sum(singularVector);
    temp2 = 0;
    socsNumber = size(U,2)-1;
    for i = 1:length(singularVector)-1
        temp2 = temp2 + singularVector(i);
        if temp2>rateThreshold*totalSingular
            socsNumber = i;
            break;
        end
    end
else
    error('Error: TCCKernalSetting.method should be setNumber of setThreshold!');
end

TCCMatrix_Kernel = zeros(waferNg-1,waferNf-1,socsNumber);


if waferNf > maskNf
    rangeWaferNf = (waferNf-maskNf+2)/2 : (waferNf+maskNf-2)/2;
    rangeMaskNf = 1:maskNf-1;
else
    rangeWaferNf =1:waferNf-1;
    rangeMaskNf =  (maskNf-waferNf+2)/2 : (waferNf+maskNf-2)/2;
end

if waferNg > maskNg
    rangeWaferNg =  (waferNg-maskNg+2)/2 : (waferNg+maskNg-2)/2;
    rangeMaskNg = 1:maskNg-1;
else
    rangeWaferNg = 1:waferNg-1;
    rangeMaskNg = (maskNg-waferNg+2)/2 : (waferNg+maskNg-2)/2;
end

temp2 = zeros(maskNg-1,maskNf-1);
temp3 = zeros(waferNg-1,waferNf-1 );

for i = 1:socsNumber
    temp1 = reshape(U(:,i),FG_ValidSize(2),FG_ValidSize(1));% 1为Y方向，2为X方向
    temp2((maskNg - FG_ValidSize(2))/2+1 : (maskNg + FG_ValidSize(2))/2, ...
        (maskNf - FG_ValidSize(1))/2+1 : (maskNf + FG_ValidSize(1))/2) = temp1;
    temp3( rangeWaferNg, rangeWaferNf ) = temp2(rangeMaskNg,rangeMaskNf);
    TCCMatrix_Kernel(:,:,i) = fftshift(temp3);
end

diagS = diag(S);
diagS = reshape(diagS(1:socsNumber),1,1,[]);
TCCMatrix_Kernel = TCCMatrix_Kernel.*sqrt(diagS);

end