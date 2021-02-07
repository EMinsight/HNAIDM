classdef Mask
    properties
        MaskType;
        MaskDefineType;%功能待定
        Orientation; %顺时针旋转
        Period_X;
        Period_Y;
        Nf;
        Ng;
        Background_Transmissivity;
        Background_Phase;
        bgTransmission;%兼容性
        bgPhase;%兼容性
        Feature;
        CutLine;
        
        % 频域定义
        Spectrum = 1
        Fx = 0
        Gy =0
    end
    
    methods
        function mk = Mask()
            mk.MaskType = '2D';
            mk.Period_X = 500;
            mk.Period_Y = 500;
            mk.Nf = 81;
            mk.Ng = 81;
            mk.Orientation = 0;
            mk.Background_Transmissivity = 0;
            mk.Background_Phase = 0;
            mk.MaskDefineType = '2D';
            mk.Feature = [];
        end
    end
    
    methods
        function  [ Hk, F, G, HkBlank ] = CalculateMaskSpectrum(mk,po,sr)
            NA = po.NA;
            Wavelength = sr.Wavelength;
            
            switch lower(mk.MaskType)
                case '1d'
                    [Hk,F, G] = mk.Mask1D2Spectrum(NA,Wavelength);
                case '1dpixel'
                    error('To be implement');
                    %  Hk = PixelMask1D2Spectrum(mk,NA,Wavelength);
                case '2d'
                    [Hk,F, G] = mk.Mask2D2Spectrum(NA,Wavelength);
                case '2dpixel'
                    [Hk,F, G] = mk.PixelMask2D2Spectrum(NA,Wavelength);
                case '3d'
                    % TODO: 三维掩模算法待完成
                case 'spectrum'
                    Hk = mk.Spectrum;
                    F = mk.Fx;
                    G = mk.Gy;
                    HkBlank = mk.Spectrum;
                otherwise
                    error('Error: There is not this type of mask!');
            end
        end
        
        function  [ Hk, F, G ] = CalculateMaskSpectrumV2(mk,po,sr)
            NA = po.NA;
            Wavelength = sr.Wavelength;
            
            switch lower(mk.MaskType)
                case '1d'
                    [Hk,F, G] = mk.Mask1D2Spectrum(NA,Wavelength);
                case '1dpixel'
                    error('To be implement');
                    %  Hk = PixelMask1D2Spectrum(mk,NA,Wavelength);
                case '2d'
                    [Hk,F, G] = mk.Mask2D2Spectrum(NA,Wavelength);
                case '2dpixel'
                    [Hk,F, G] = mk.PixelMask2D2Spectrum(NA,Wavelength);
                case '3d'
                    % TODO: 三维掩模算法待完成
                case 'spectrum'
                    Hk = mk.Spectrum;
                    F = mk.Fx;
                    G = mk.Gy;
                otherwise
                    error('Error: There is not this type of mask!');
            end
        end
        

        
        function mk = Convert2DiscreteMask(mk)
            X = linspace(-mk.Period_X/2,mk.Period_X/2,mk.Nf);
            Y = linspace(-mk.Period_Y/2,mk.Period_Y/2,mk.Ng);
            [X,Y] = meshgrid(X,Y);
            maskDiscreteData = zeros(size(X));
            inArea = zeros(size(X)); % 对所有掩模点的位置进行判断
            inArea = (inArea==1); % 这句话有什么意义？每次都是对整个矩阵先置零，怎么可能还有为1的位置呢？
            bgComplexAmplitude =  mk.Background_Transmissivity*exp(1j*mk.Background_Phase); % 掩模背景复振幅
            featureNumber = size(mk.Feature,2);
            for k = 1:featureNumber
                if strcmpi(mk.Feature(k).ShapeType,'r')    % 构成掩模的子图形为矩形,在后续编程时，可直接在循环中求出区域内的点，循环外一次性求复透过率
                    am = mk.Feature(k).ComplexAm(1) * exp(1j*mk.Feature(k).ComplexAm(2)); % 第k个子图形对应的复透过率
                    xVertexVector = [mk.Feature(k).BoundaryVertexX(1), mk.Feature(k).BoundaryVertexX(2), mk.Feature(k).BoundaryVertexX(2), mk.Feature(k).BoundaryVertexX(1)];
                    yVertexVector = [mk.Feature(k). BoundaryVertexY(1), mk.Feature(k). BoundaryVertexY(1), mk.Feature(k). BoundaryVertexY(2), mk.Feature(k). BoundaryVertexY(2)];
                    inRectangle = inpolygon(X,Y,xVertexVector,yVertexVector);
                    doubleCount = inArea & inRectangle; % 两者重合部分，表示计算了两次，需要在后面去掉一次;*********就是为了避免子图形之间有边是重合的************。
                    inRectangle(doubleCount) = 0; % 除去已经记载的区域内点，保证此时的inRectangle为每次新加入的点
                    inArea = inArea | inRectangle; % 存放已经做出判断的点
                    maskDiscreteData = maskDiscreteData + (am - bgComplexAmplitude).*double(inRectangle); %累加得到总的掩模图案
                elseif strcmpi(mk.Feature(k).ShapeType,'t')    % 构成掩模的子图形为三角形
                    am = mk.Feature(k).ComplexAm(1) * exp(1j*mk.Feature(k).ComplexAm(2));
                    xVertexVector = mk.Feature(k).BoundaryVertexX;
                    yVertexVector = mk.Feature(k). BoundaryVertexY;
                    inTriangle = inpolygon(X,Y,xVertexVector,yVertexVector);
                    doubleCount  = inArea & inTriangle;
                    inTriangle(doubleCount) = 0;
                    inArea = inArea | inTriangle;
                    maskDiscreteData = maskDiscreteData + (am - bgComplexAmplitude).*double(inTriangle);
                elseif strcmpi(mk.Feature(k).ShapeType,'p')    % 构成掩模的子图形为平行四边形 parallelogram
                    am = mk.Feature(k).ComplexAm(1) * exp(1j*mk.Feature(k).ComplexAm(2));
                    xVertexVector = mk.Feature(k).BoundaryVertexX;
                    yVertexVector = mk.Feature(k). BoundaryVertexY;
                    inParallelogram = inpolygon(X,Y,xVertexVector,yVertexVector);
                    doubleCount  = inArea & inParallelogram;
                    inParallelogram(doubleCount) = 0;
                    inArea = inArea | inParallelogram;
                    maskDiscreteData = maskDiscreteData + (am - bgComplexAmplitude).*double(inParallelogram);
                elseif strcmpi(mk.Feature(k).ShapeType,'c')    % 构成掩模的子图形为圆形
                    am = mk.Feature(k).ComplexAm(1) * exp(1j*mk.Feature(k).ComplexAm(2));
                    xCenterVector = mk.Feature(k).BoundaryVertexX(1);
                    yCenterVector = mk.Feature(k).BoundaryVertexX(2);
                    r = mk.Feature(k). BoundaryVertexY;
                    distance2 = (X - xCenterVector).^2 + (Y - yCenterVector).^2;
                    distance = sqrt(distance2);
                    inCircle = (distance - r)<1e-6;
                    doubleCount = inArea & inCircle;
                    inCircle(doubleCount) = 0;
                    inArea = inArea | inCircle;
                    maskDiscreteData = maskDiscreteData + (am - bgComplexAmplitude).*double(inCircle);
                end
            end
            maskDiscreteData = maskDiscreteData + bgComplexAmplitude; % 将多个子图形集成之后，再与背景相加，得到完整的掩模的复透过率
            
            mk.Feature = maskDiscreteData;
            mk.MaskType = '2DPixel';
        end
        
        function mk = Convert2PSM( mk, type, threshold)
            %             switch length( varargin )
            %
            %                 case 1
            %                     backgroundPhase = varargin{1};
            %                 case 2
            %
            %                 otherwise
            %             end
            if strcmpi(type,'binary')
                mk.Feature = abs(mk.Feature)>threshold(1);
            elseif strcmpi(type,'attpsm')
                mkMap = double(abs( mk.Feature) > threshold(1));
                mkMap( mkMap< threshold ) = mk.Background_Transmissivity*exp(1j*mk.Background_Phase);
                mk.Feature = mkMap;
            elseif strcmpi(type,'altpsm')
                mkMap = double( real(mk.Feature) > threshold(1));
                mkMap( real(mk.Feature) < threshold(2) ) = exp(1j*pi);
                mk.Feature = mkMap;
            end
        end
    end
    
    methods (Access = private)
        
        function [Spectrum,f,g] = Mask1D2Spectrum(mk,NA,Wavelength)
            
            normalized_Frequency = NA/Wavelength; % 进行傅里叶变换之前的系数表示
            bgComplexAmplitude = mk.Background_Transmissivity*exp(1j* mk.Background_Phase);
            normalized_Period_X = mk.Period_X/(Wavelength/NA); % 归一化采样间距
            f = (1/normalized_Period_X) * (-(mk.Nf-1)/2 : 1 : (mk.Nf-1)/2);
            g = 0;
            
            Spectrum = zeros(1,length(f));
            SubRectNum = length(mk.Feature);
            for i = 1:SubRectNum
                if strcmpi(mk.Feature(i).ShapeType,'r')
                    width = abs(mk.Feature(i).BoundaryVertexX(2) - mk.Feature(i).BoundaryVertexX(1)); % 图形宽度
                    
                    am = mk.Feature(i).ComplexAm(1) * exp(1j * mk.Feature(i).ComplexAm(2));
                    center = (mk.Feature(i).BoundaryVertexX(2) + mk.Feature(i).BoundaryVertexX(1))/2; % 图形的中心位置
                    temp = normalized_Frequency * (am - bgComplexAmplitude)*width.*exp(-1j*2*pi*center*normalized_Frequency.*f).*sinc(width*normalized_Frequency.*f);
                    Spectrum = Spectrum + temp;
                end
            end
      
            bgSpectrum = zeros(1,length(f)); % 背景频谱
            bgSpectrum(abs(f)<1e-6) = normalized_Period_X.* bgComplexAmplitude;
            Spectrum = Spectrum + bgSpectrum;
        end
        
        function [Hk,f,g] = Mask2D2Spectrum(mk,NA,Wavelength)
            bgComplexAmplitude=mk.Background_Transmissivity*exp(1j* mk.Background_Phase);
            normalized_Period_X = mk.Period_X/(Wavelength/NA); % 归一化采样间距
            normalized_Period_Y = mk.Period_Y/(Wavelength/NA); % 像方折射率已经放在NA中考虑进来了
            f = (1/normalized_Period_X) * (-(mk.Nf-1)/2 : 1 : (mk.Nf-1)/2);
            g = (1/normalized_Period_Y) * (-(mk.Ng-1)/2 : 1 : (mk.Ng-1)/2);
            normalized_Frequency=NA/Wavelength;
            
            Hk=zeros( length(g),length(f) );
            SubRectNum = length(mk.Feature);
            for ii=1:SubRectNum
                
                if strcmpi( mk.Feature(ii).ShapeType ,'r' )  %矩形图形
                    % 简化图形设计 只需定义对角点
                    xWidth=abs(mk.Feature(ii).BoundaryVertexX(2) - mk.Feature(ii).BoundaryVertexX(1));
                    yWidth=abs(mk.Feature(ii).BoundaryVertexY(2) - mk.Feature(ii).BoundaryVertexY(1)) ;
                    
                    a = mk.Feature(ii).ComplexAm(1)  *  exp( 1j*mk.Feature(ii).ComplexAm(2) );
                    
                    xCenter=( mk.Feature(ii).BoundaryVertexX(1) + mk.Feature(ii).BoundaryVertexX(2) ) /2;
                    yCenter=( mk.Feature(ii).BoundaryVertexY(1) + mk.Feature(ii).BoundaryVertexY(2) ) /2;
                    tempx=   exp(-1i*2*pi*xCenter*normalized_Frequency.*f).*sinc(xWidth*normalized_Frequency.*f)  ;
                    tempy=   exp(-1i*2*pi*yCenter*normalized_Frequency.*g).*sinc(yWidth*normalized_Frequency.*g)  ;
                    temp=(a-bgComplexAmplitude)*normalized_Frequency^2*xWidth*yWidth.* ( flipud(tempy') * tempx ) ;
                    Hk=Hk+temp;
                elseif strcmpi( mk.Feature(ii).ShapeType ,'c' ) %圆形图形
                    %besselj
                    a=mk.Feature(ii).ComplexAm(1)  *  exp(   1j*mk.Feature(ii).ComplexAm(2)   );
                    wr= mk.Feature(ii).BoundaryVertexY(1);
                    xCenter= mk.Feature(ii).BoundaryVertexX(1) ;
                    yCenter= mk.Feature(ii).BoundaryVertexX(2) ;
                    [f2,g2]=meshgrid(f,g);
                    rho = hypot(f2,g2);
                    rho = rho+1e-10;
                    tempx = exp(-1i*2*pi*xCenter*normalized_Frequency.*f);
                    tempy = exp(-1i*2*pi*yCenter*normalized_Frequency.*g);
                    temp2D= flipud(tempy') * tempx;
                    tempBessel=besselj(1,2*pi*wr*normalized_Frequency.*rho)./rho;
                    temp=(a-bgComplexAmplitude)*tempBessel.*temp2D;
                    Hk=Hk+temp;
                elseif strcmpi( mk.Feature(ii).ShapeType ,'t' ) %三角形图形
                    %注意归一化的问题
                    a=mk.Feature(ii).ComplexAm(1)  *  exp(   1j*mk.Feature(ii).ComplexAm(2)   );
                    
                    trx1=mk.Feature(ii).BoundaryVertexX(1);
                    trx2=mk.Feature(ii).BoundaryVertexX(2);
                    trx3=mk.Feature(ii).BoundaryVertexX(3);
                    
                    try1=mk.Feature(ii).BoundaryVertexY(1);
                    try2=mk.Feature(ii).BoundaryVertexY(2);
                    try3=mk.Feature(ii).BoundaryVertexY(3);
                    
                    trxv=[trx1,trx2,trx3];
                    tryv=[try1,try2,try3];
                    temp=triMask2spectrum(trxv,tryv,a,bgComplexAmplitude,normalized_Frequency,f,g);
                    
                    Hk=Hk+temp;
                    
                elseif strcmpi( mk.Feature(ii).ShapeType ,'p' )%单点
                    
                    a=mk.Feature(ii).ComplexAm(1)  *  exp(   1j*mk.Feature(ii).ComplexAm(2)   );
                    polyX=mk.Feature(ii).BoundaryVertexX(:);
                    polyY=mk.Feature(ii).BoundaryVertexY(:);
                    %         dt = DelaunayTri(polyX,polyY);
                    dt = delaunayTriangulation(polyX,polyY);
                    for iipoly=1:size(dt,1)
                        
                        trxv= polyX( dt(iipoly,:) ) ;
                        tryv= polyY( dt(iipoly,:) ) ;
                        graycenter=[ mean( trxv ),mean(  tryv ) ];
                        inone = inpolygon( graycenter(1),graycenter(2),polyX,polyY );
                        if inone
                            temp=triMask2spectrum(trxv,tryv,a,bgComplexAmplitude,normalized_Frequency,f,g);
                            Hk=Hk+temp;
                        end
                    end
                end
            end
            
  
            
            BkSpectrum = zeros(  length(g),length(f)  );
            BkSpectrum(  abs(g)<1e-9, abs(f)<1e-9 ) = (normalized_Period_X*normalized_Period_Y).*bgComplexAmplitude  ;%delt(f,g)=1/df/dg
            Hk=Hk+BkSpectrum;
            
        end
        
        function [Hk,f,g] = PixelMask2D2Spectrum(mk,NA,Wavelength)
            
            h= mk.Feature;
            [ PntNum_Y, PntNum_X ]=size(h);%PntNum is odd
            
            if mk.Nf ~= PntNum_X
                error('error: mask data size 2 must be equal to Mask.Nf ');
            end
            
            if mk.Ng ~= PntNum_Y
                error('error: mask data size 1 must be equal to Mask.Ng ');
            end
            
            period_X = mk.Period_X;
            period_Y = mk.Period_Y;
            Start_X = -mk.Period_X/2;
            Start_Y = -mk.Period_Y/2;
            
            normPitchX=period_X/(Wavelength/NA);
            normPitchY=period_Y/(Wavelength/NA);
            normSX=Start_X/(Wavelength/NA);
            normSY=Start_Y/(Wavelength/NA);
            
            f=(1/normPitchX)*(-(PntNum_X-1)/2:(PntNum_X-1)/2);
            g=(1/normPitchY)*(-(PntNum_Y-1)/2:(PntNum_Y-1)/2);
            
            normdx=normPitchX/(PntNum_X-1);
            normdy=normPitchY/(PntNum_Y-1);
            
            normDfx=1/normPitchX;
            normDfy=1/normPitchY;
            normfx=normDfx.*( -(PntNum_X-1)/2 :  (PntNum_X-1)/2 );
            normfy=normDfy.*( -(PntNum_Y-1)/2 :  (PntNum_Y-1)/2 );
            nDFTVectorX=0:PntNum_X-1;
            nDFTVectorY=0:PntNum_Y-1;
            
            %AFactor：PntNum_X行，PntNum_Y列
            AFactor1=exp( -1i*2*pi*normSY.*normfy(:) ) * exp( -1i*2*pi*normSX.*normfx );
            AFactor=(normdx*normdy).*AFactor1;
            
            hFactor=exp(1i*pi.*nDFTVectorY(:) )*exp(1i*pi.*nDFTVectorX );
            hh=h.*hFactor;
            
            h2D=hh(1:end-1,1:end-1);
            HYEnd=hh(end,1:end-1);
            HXEnd=hh(1:end-1,end);
            HXYEnd=hh(end,end);
            
            DFTH_2D= fft2(h2D);
            DFTH_2D(:,end+1)=DFTH_2D(:,1);
            DFTH_2D(end+1,:)=DFTH_2D(1,:);
            
            DFTHYEnd= fft(HYEnd);
            DFTHYEnd(end+1)=DFTHYEnd(1);
            DFTHM_12D=ones(PntNum_Y,1)*DFTHYEnd;
            
            DFTHXEnd= fft(HXEnd);
            DFTHXEnd(end+1)=DFTHXEnd(1);
            DFTHN_12D=DFTHXEnd(:)*ones(1,PntNum_X);
            
            Hk=AFactor.*( DFTH_2D+ DFTHM_12D+DFTHN_12D+HXYEnd );
            
            normalized_Frequency=NA/Wavelength;
            tempx=   sinc(mk.Period_X*normalized_Frequency.*f)  ;
            tempy=   sinc(mk.Period_Y*normalized_Frequency.*g)  ;
            HkBlank = normalized_Frequency^2 * mk.Period_X*mk.Period_Y.* (tempy' * tempx ) ;
            
        end
        
    end
    
    
    methods
        function mask = CreateSpaceMask(mask, lineCD, linePitch)
            mask.MaskType = '1D';
            mask.bgTransmission=0.0;
            mask.bgPhase=0;
            mask.Background_Transmissivity = 0;
            mask.Background_Phase = 0;
            mask.Period_X=linePitch;%原来为1000
            mask.Period_Y=linePitch;%原来为1000
            
            mask.Feature(1).ShapeType='r';% 点的坐标：先左上，再右下
            mask.Feature(1). BoundaryVertexX=[-lineCD/2, lineCD/2];%Boundary Vertex X&Y=[-400 -300];
            mask.Feature(1). BoundaryVertexY=[-linePitch/2, linePitch/2];%[400 200 ];
            mask.Feature(1). ComplexAm=[1 0];
            
            %设定 cutline
            mask.CutLine(1).X = [-linePitch/2, linePitch/2];
            mask.CutLine(1).Y = [0, 0];
        end
        
        function mask = CreateLineMask(mask, spaceCD, spacePitch)
            mask.MaskType = '1D';
            mask.bgTransmission=0.0;
            mask.bgPhase=0;
            mask.Background_Transmissivity = 1;
            mask.Background_Phase = 0;
            mask.Period_X = spacePitch;%原来为1000
            mask.Period_Y = spacePitch;%原来为1000
            
            mask.Feature(1).ShapeType='r';% 点的坐标：先左上，再右下
            mask.Feature(1). BoundaryVertexX=[-spaceCD/2, spaceCD/2];%Boundary Vertex X&Y=[-400 -300];
            mask.Feature(1). BoundaryVertexY=[-spacePitch/2, spacePitch/2];%[400 200 ];
            mask.Feature(1). ComplexAm=[0 0];
            
            %设定 cutline
            mask.CutLine(1).X = [-spacePitch/2, spacePitch/2];
            mask.CutLine(1).Y = [0, 0];
        end
    end
    
    methods (Static)
        
        function mask = CreateMask(maskType,varargin)
            % 创建初始Mask，内置Mask参数
            % 常见掩模图案
            mask = Mask();
            % 1D 掩模
            switch lower(maskType )
                case 'line'
                    switch length(varargin)
                        case 0
                            lineCD = 45;
                            linePitch = 90;
                        case 1
                            lineCD = varargin{1};
                            linePitch = lineCD*2;
                        case 2
                            lineCD = varargin{1};
                            linePitch = varargin{2};
                    end
                    mask = mask.CreateLineMask(lineCD,linePitch);
                    return;
                    
                case 'space'
                    switch length(varargin)
                        case 0
                            spaceCD = 45;
                            spacePitch = 90;
                        case 1
                            spaceCD = varargin{1};
                            spacePitch = spaceCD*2;
                        case 2
                            spaceCD = varargin{1};
                            spacePitch = varargin{2};
                    end
                    mask = mask.CreateSpaceMask( spaceCD, spacePitch );
                    return;
                    
                case 'space_alt'
                    switch length(varargin)
                        case 0
                            spaceCD = 45;
                            spacePitch = 90;
                        case 1
                            spaceCD = varargin{1};
                            spacePitch = spaceCD*2;
                        case 2
                            spaceCD = varargin{1};
                            spacePitch = varargin{2};
                    end
                    maskXL =  spacePitch*2;
                    maskYL =  spacePitch*2;
                    
                    mask.MaskType = '1D';
                    mask.Background_Transmissivity = 0;
                    mask.Background_Phase = 0;
                    mask.Period_X = maskXL;%原来为1000
                    mask.Period_Y = maskYL;%原来为1000
                    
                    mask.Feature(1).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(1). BoundaryVertexX=[-( spacePitch + spaceCD )/2, -( spacePitch - spaceCD )/2];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(1). BoundaryVertexY=[-maskYL/2, maskYL/2];%[400 200 ];
                    mask.Feature(1). ComplexAm=[1 0];
                    
                    mask.Feature(2).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(2). BoundaryVertexX=[( spacePitch - spaceCD )/2, ( spacePitch + spaceCD )/2];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(2). BoundaryVertexY=[-maskYL/2, maskYL/2];%[400 200 ];
                    mask.Feature(2). ComplexAm=[1 pi];
                    
                    %设定 cutline
                    mask.CutLine(1).X = [ 0, maskXL ];
                    mask.CutLine(1).Y = [ 0, 0 ];
                    return;
            end

            mask.bgTransmission=0.0;
            mask.bgPhase=0;
            mask.Background_Transmissivity = 0;
            mask.Background_Phase = 0;
            switch lower(maskType)                    
                case 'crossgate'%
                    mask.Period_X = 420;%原来为1000
                    mask.Period_Y = 420;%原来为1000
                    %rectangle
                    mask.Feature(1).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(1). BoundaryVertexX=[-180 -135];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(1). BoundaryVertexY=[-210 210];%[400 200 ];
                    mask.Feature(1). ComplexAm=[1 0];
                    % %
                    %rectangle
                    mask.Feature(2).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(2). BoundaryVertexX=[-75 -30];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(2). BoundaryVertexY=[-25 210];%[400 200 ];
                    mask.Feature(2). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(3).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(3). BoundaryVertexX=[-75 -30];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(3). BoundaryVertexY=[-210,-70 ];%[400 200 ];
                    mask.Feature(3). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(4).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(4). BoundaryVertexX=[30 75];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(4). BoundaryVertexY=[-210 -25];%[400 200 ];
                    mask.Feature(4). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(5).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(5). BoundaryVertexX=[-30 75];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(5). BoundaryVertexY=[-25 25];%[400 200 ];
                    mask.Feature(5). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(6).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(6). BoundaryVertexX=[30 75];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(6). BoundaryVertexY=[70 210];%[400 200 ];
                    mask.Feature(6). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(7).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(7). BoundaryVertexX=[135 180];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(7). BoundaryVertexY=[-210 210];%[400 200 ];
                    mask.Feature(7). ComplexAm=[1 0];
                    
                    % 设定cutline
                    mask. CutLine(1).X = [0,105];
                    mask. CutLine(1).Y = [100,100];
                    
                    mask. CutLine(2).X = [0,0];
                    mask. CutLine(2).Y = [-100,100];
                    
                    mask. CutLine(3).X = [-105,0];
                    mask. CutLine(3).Y = [100,100];
                    
                    mask. CutLine(4).X = [105,210];
                    mask. CutLine(4).Y = [0,0];
                    
                    
                case 'contact_holes' %periodic array of contact holes
                    mask.Period_X = 210;%原来为1000
                    mask.Period_Y = 210;%原来为1000
                    %rectangle
                    mask.Feature(1).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(1). BoundaryVertexX=[30 75];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(1). BoundaryVertexY=[30 75];%[400 200 ];
                    mask.Feature(1). ComplexAm=[1 0];
                    % %
                    %rectangle
                    mask.Feature(2).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(2). BoundaryVertexX=[-75 -30];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(2). BoundaryVertexY=[-75 -30];%[400 200 ];
                    mask.Feature(2). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(3).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(3). BoundaryVertexX=[30 75];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(3). BoundaryVertexY=[-75 -30];%[400 200 ];
                    mask.Feature(3). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(4).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(4). BoundaryVertexX=[-75 -30];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(4). BoundaryVertexY=[30 75];%[400 200 ];
                    mask.Feature(4). ComplexAm=[1 0];
                case 'line_space' %line_space
                    mask.Period_X = 720;%原来为1000
                    mask.Period_Y = 720;%原来为1000
                    %rectangle
                    mask.Feature(1).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(1). BoundaryVertexX=[-22.5 22.5];%Boundary Vertex X&Y=[-400
                    mask.Feature(1). BoundaryVertexY=[-300 300];%[400 200 ];
                    mask.Feature(1). ComplexAm=[1 0];
                    %rectangle
                    mask.Feature(2).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(2). BoundaryVertexX=[67.5 112.5];%Boundary Vertex X&Y=[-400
                    mask.Feature(2). BoundaryVertexY=[-300 300];%[400 200 ];
                    mask.Feature(2). ComplexAm=[1 0];
                    %rectangle
                    mask.Feature(3).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(3). BoundaryVertexX=[157.5 202.5];%Boundary Vertex X&Y=[-
                    mask.Feature(3). BoundaryVertexY=[-300 300];%[400 200 ];
                    mask.Feature(3). ComplexAm=[1 0];
                    %rectangle
                    mask.Feature(4).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(4). BoundaryVertexX=[-112.5 -67.5];%Boundary Vertex X&Y=[-
                    mask.Feature(4). BoundaryVertexY=[-300 300];%[400 200 ];
                    mask.Feature(4). ComplexAm=[1 0];
                    %rectangle
                    mask.Feature(5).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(5). BoundaryVertexX=[-202.5 -157.5];%Boundary Vertex X&Y=
                    mask.Feature(5). BoundaryVertexY=[-300 300];%[400 200 ];
                    mask.Feature(5). ComplexAm=[1 0];
                    
                case 'complex'
                    
                    mask.bgTransmission=0.0;
                    mask.bgPhase=0;
                    mask.Period_X = 1200;%原来为1000
                    mask.Period_Y = 1200;%原来为1000
                    %rectangle
                    %rectangle
                    mask.Feature(1).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(1). BoundaryVertexX=[-510 510];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(1). BoundaryVertexY=[427.5 472.5];%[400 200 ];
                    mask.Feature(1). ComplexAm=[1 0];
                    %%
                    %rectangle
                    mask.Feature(2).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(2). BoundaryVertexX=[240 465];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(2). BoundaryVertexY=[337.5 382.5];%[400 200 ];
                    mask.Feature(2). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(3).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(3). BoundaryVertexX=[240 345];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(3). BoundaryVertexY=[292.5 337.5];%[400 200 ];
                    mask.Feature(3). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(4).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(4). BoundaryVertexX=[120 345];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(4). BoundaryVertexY=[247.5 292.5];%[400 200 ];
                    mask.Feature(4). ComplexAm=[1 0];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                    %rectangle
                    mask.Feature(5).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(5). BoundaryVertexX=[120 345];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(5). BoundaryVertexY=[157.5 202.5];%[400 200 ];
                    mask.Feature(5). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(6).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(6). BoundaryVertexX=[240 345];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(6). BoundaryVertexY=[112.5 157.5];%[400 200 ];
                    mask.Feature(6). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(7).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(7). BoundaryVertexX=[240 465];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(7). BoundaryVertexY=[67.5 112.5];%[400 200 ];
                    mask.Feature(7). ComplexAm=[1 0];
                    %%
                    %rectangle
                    mask.Feature(8).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(8). BoundaryVertexX=[-345 -120];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(8). BoundaryVertexY=[247.5 292.5];%[400 200 ];
                    mask.Feature(8). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(9).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(9). BoundaryVertexX=[-345 -240];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(9). BoundaryVertexY=[292.5 337.5];%[400 200 ];
                    mask.Feature(9). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(10).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(10). BoundaryVertexX=[-465 -240];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(10). BoundaryVertexY=[337.5 382.5];%[400 200 ];
                    mask.Feature(10). ComplexAm=[1 0];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                    %rectangle
                    mask.Feature(11).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(11). BoundaryVertexX=[-465 -240];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(11). BoundaryVertexY=[67.5 112.5];%[400 200 ];
                    mask.Feature(11). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(12).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(12). BoundaryVertexX=[-345 -240];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(12). BoundaryVertexY=[112.5 157.5];%[400 200 ];
                    mask.Feature(12). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(13).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(13). BoundaryVertexX=[-345 -120];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(13). BoundaryVertexY=[157.5 202.5];%[400 200 ];
                    mask.Feature(13). ComplexAm=[1 0];
                    %%
                    mask.Feature(14).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(14). BoundaryVertexX=[-510 510];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(14). BoundaryVertexY=[-22.5 22.5];%[400 200 ];
                    mask.Feature(14). ComplexAm=[1 0];
                    %%
                    %rectangle
                    mask.Feature(15).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(15). BoundaryVertexX=[-345 -120];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(15). BoundaryVertexY=[-202.5 -157.5];%[400 200 ];
                    mask.Feature(15). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(16).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(16). BoundaryVertexX=[-345 -240];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(16). BoundaryVertexY=[-157.5 -112.5];%[400 200 ];
                    mask.Feature(16). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(17).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(17). BoundaryVertexX=[-465 -240];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(17). BoundaryVertexY=[-112.5 -67.5];%[400 200 ];
                    mask.Feature(17). ComplexAm=[1 0];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                    %rectangle
                    mask.Feature(18).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(18). BoundaryVertexX=[-465 -240];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(18). BoundaryVertexY=[-382.5 -337.5];%[400 200 ];
                    mask.Feature(18). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(19).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(19). BoundaryVertexX=[-345 -240];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(19). BoundaryVertexY=[-337.5 -292.5];%[400 200 ];
                    mask.Feature(19). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(20).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(20). BoundaryVertexX=[-345 -120];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(20). BoundaryVertexY=[-292.5 -247.5];%[400 200 ];
                    mask.Feature(20). ComplexAm=[1 0];
                    %%
                    %rectangle
                    mask.Feature(21).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(21). BoundaryVertexX=[240 465];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(21). BoundaryVertexY=[-112.5 -67.5];%[400 200 ];
                    mask.Feature(21). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(22).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(22). BoundaryVertexX=[240 345];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(22). BoundaryVertexY=[-157.5 -112.5];%[400 200 ];
                    mask.Feature(22). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(23).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(23). BoundaryVertexX=[120 345];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(23). BoundaryVertexY=[-202.5 -157.5];%[400 200 ];
                    mask.Feature(23). ComplexAm=[1 0];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                    %rectangle
                    mask.Feature(24).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(24). BoundaryVertexX=[120 345];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(24). BoundaryVertexY=[-292.5 -247.5];%[400 200 ];
                    mask.Feature(24). ComplexAm=[1 0];
                    %rectangle
                    mask.Feature(25).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(25). BoundaryVertexX=[240 345];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(25). BoundaryVertexY=[-337.5 -292.5];%[400 200 ];
                    mask.Feature(25). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(26).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(26). BoundaryVertexX=[240 465];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(26). BoundaryVertexY=[-382.5 -337.5];%[400 200 ];
                    mask.Feature(26). ComplexAm=[1 0];
                    %%
                    mask.Feature(27).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(27). BoundaryVertexX=[-510 510];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(27). BoundaryVertexY=[-472.5 -427.5];%[400 200 ];
                    mask.Feature(27). ComplexAm=[1 0];
                    
                    
                case 'mypattern'
                    mask.bgPhase=0;
                    mask.Period_X = 600;%原来为1000
                    mask.Period_Y = 600;%原来为1000
                    
                    %rectangle
                    mask.Feature(1).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(1). BoundaryVertexX=[ -200 200 ];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(1). BoundaryVertexY=[ 172 200 ];%[400 200 ];
                    mask.Feature(1). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(2).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(2). BoundaryVertexX=[172 200];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(2). BoundaryVertexY=[-200  172];%[400 200 ];
                    mask.Feature(2). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(3).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(3). BoundaryVertexX=[-200 -140];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(3). BoundaryVertexY=[ 90, 118 ];%[400 200 ];
                    mask.Feature(3). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(4).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(4). BoundaryVertexX=[ -80 -20];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(4). BoundaryVertexY=[ 90 118 ];%[400 200 ];
                    mask.Feature(4). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(5).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(5). BoundaryVertexX=[ -200 100 ];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(5). BoundaryVertexY=[ -200 -172 ];%[400 200 ];
                    mask.Feature(5). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(6).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(6). BoundaryVertexX=[-200 -172];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(6). BoundaryVertexY=[-172 0 ];%[400 200 ];
                    mask.Feature(6). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(7).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(7). BoundaryVertexX=[-132 -104];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(7). BoundaryVertexY=[-172 0 ];%[400 200 ];
                    mask.Feature(7). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(8).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(8). BoundaryVertexX=[-64 -36];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(8). BoundaryVertexY=[-172 0 ];%[400 200 ];
                    mask.Feature(8). ComplexAm=[1 0];
                    
                    %rectangle
                    mask.Feature(9).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(9). BoundaryVertexX=[72 100];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(9). BoundaryVertexY=[-172 -50 ];%[400 200 ];
                    mask.Feature(9). ComplexAm=[1 0];
                    
                case 'sram'
                    maskSize = 400;
                    mask.bgTransmission=0.0;
                    mask.bgPhase=0;
                    mask.Period_X=maskSize;%原来为1000
                    mask.Period_Y=maskSize;%原来为1000
                    
                    mask.Feature(1).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(1). BoundaryVertexX=[-510 510];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(1). BoundaryVertexY=[-472.5 -427.5];%[400 200 ];
                    mask.Feature(1). ComplexAm=[1 0];
                    
                    
                otherwise
                    error('error:This type of mask has not been included');
            end
            
        end
        
        function mask = CreateparameterizedMask(maskType,varargin)
            % 创建初始Mask，内置Mask参数
            % 常见掩模图案
            mask = Mask();
            % 1D 掩模
            switch lower(maskType )
                case 'line'
                    switch length(varargin)
                        case 0
                            lineCD = 45;
                            linePitch = 90;
                        case 1
                            lineCD = varargin{1};
                            linePitch = lineCD*2;
                        case 2
                            lineCD = varargin{1};
                            linePitch = varargin{2};
                    end
                    mask = mask.CreateLineMask(lineCD,linePitch);
                case 'space'
                    switch length(varargin)
                        case 0
                            spaceCD = 45;
                            spacePitch = 90;
                        case 1
                            spaceCD = varargin{1};
                            spacePitch = spaceCD*2;
                        case 2
                            spaceCD = varargin{1};
                            spacePitch = varargin{2};
                    end
                    mask = mask.CreateSpaceMask( spaceCD, spacePitch );
                case 'space_end_dense'
                    
                    switch length(varargin)
                        case 0
                            gapCD = 45;
                            spaceCD = 45;
                            spacePitch = 1000;
                        case 1
                            gapCD = varargin{1};
                            spaceCD = 45;
                            spacePitch = spaceCD*2;
                        case 2
                            gapCD = varargin{1};
                            spaceCD = varargin{2};
                            spacePitch = spaceCD*2;
                        case 3
                            gapCD = varargin{1};
                            spaceCD = varargin{2};
                            spacePitch = varargin{3};
                    end
                    
                    mask = Mask();
                    mask.bgTransmission = 0.0;
                    mask.bgPhase=0;
                    mask.Background_Transmissivity = 0;
                    mask.Background_Phase = 0;
                    mask.Period_X = spacePitch;%默认90
                    mask.Period_Y = spacePitch;%默认90
                    mask.Feature(1).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(1). BoundaryVertexX=[-mask.Period_X/2, -gapCD/2];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(1). BoundaryVertexY=[-spaceCD/2, spaceCD/2];%[400 200 ];
                    mask.Feature(1). ComplexAm=[1 0];
                    
                    mask.Feature(2).ShapeType='r';% 点的坐标：先左上，再右下
                    mask.Feature(2). BoundaryVertexX=[gapCD/2, mask.Period_X/2];%Boundary Vertex X&Y=[-400 -300];
                    mask.Feature(2). BoundaryVertexY=[-spaceCD/2, spaceCD/2];%[400 200 ];
                    mask.Feature(2). ComplexAm=[1 0];
                case 'space_end_t'
                    
            end
        end
        function temp=triMask2spectrum(trxv,tryv,ComplexAm,BgComplexAm,normFre,f,g)
            
            trx1=trxv(1);
            trx2=trxv(2);
            trx3=trxv(3);
            
            try1=tryv(1);
            try2=tryv(2);
            try3=tryv(3);
            
            [ff,gg]=meshgrid(f,g);
            
            newf2=(trx2-trx1)*normFre.*ff+(try2-try1)*normFre.*gg;
            newg2=(trx3-trx1)*normFre.*ff+(try3-try1)*normFre.*gg;
            
            detTri=normFre*normFre*( (trx2-trx1)*(try3-try1)-(trx3-trx1)*(try2-try1) );
            
            triPhase=exp(-1i*2*pi.*(trx1*normFre.*ff+try1*normFre.*gg) );
            
            Hfg=unit_tri_fft(newf2 , newg2);
            
            temp=detTri.*triPhase.*Hfg;
            temp=(ComplexAm-BgComplexAm).*temp;
            
        end
    end
end