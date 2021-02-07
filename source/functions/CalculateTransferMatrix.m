function [M11,M12,M21,M22] = CalculateTransferMatrix(layers,sinThetaInc,indexInc,wavelength,polMode)
% 计算传输矩阵

% 没有增透膜和抗反膜时的TE传输矩阵
M11 = 1;
M12 = 0;
M21 = 0;
M22 = 1;

if ~isempty(layers) %有膜层
    eta0 = sqrt(Numerics.Mu0 / Numerics.Epsilon0);
    for iARC = 1:length(layers)
        ly = layers{iARC};
        indexARC = ly.IndexComplex; %need complex value
        thicknessARC = ly.Thickness;
        %sinThetaInner = sinThetaInc * indexInc / indexARC;
        cosThetaInner = sqrt(1-(sinThetaInc * indexInc / indexARC).^2);
        
        etaInner = eta0/indexARC; % 非磁介质近似
        switch polMode
            case 'TE'
                faiInner = cosThetaInner / etaInner;
            case 'TM'
                faiInner = -cosThetaInner * etaInner;
            otherwise
                error('Wrong polarization mode');
        end
        ndk0ct = 2*pi/wavelength*indexARC*thicknessARC*cosThetaInner;
        jsin_ndk0ct = 1j.*sin(ndk0ct);% can not replace by sqrt(1-sin2)
         
        m11 = cos(ndk0ct);
        m22 = m11;
        %         m12 = 1j ./ faiInner .* sin(ndk0ct);
        %         m21 = 1j .* faiInner .* sin(ndk0ct);
        m12 = jsin_ndk0ct./ faiInner;
        m21 = jsin_ndk0ct.* faiInner;
        
        ml11 = M11;
        ml12 = M12;
        ml21 = M21;
        ml22 = M22;
        
        %多层传输矩阵相乘
        M11 = ml11.*m11 + ml12.*m21;
        M12 = ml11.*m12 + ml12.*m22;
        M21 = ml21.*m11 + ml22.*m21;
        M22 = ml21.*m12 + ml22.*m22;
    end
end
end