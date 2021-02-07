classdef ImagingView
    
    methods (Static)
        
        % 绘制掩模衍射谱
        function PlotMaskSpectrum(mask,source,projector)
            %% 掩模衍射级次绘图

            %Nf,Ng：归一化频率的倍数，对应最大衍射级次
            %Nf→直接计算出来的像面点数：决定直接算出来的间隔→与仿真精度有关
            %Nf有最小值：要求f的范围至少为-2~2,保证在光瞳偏移范围内
            %采用零级谱进行归一化

%             Nf = 1001;
%             Ng = 1001;
%             mask.Nf = Nf;
%             mask.Ng = Ng;

            [ Spectrum, fx, fy] = mask.CalculateMaskSpectrum(projector,source);
            %%%%%%&&&&&%%%%%%
            if  strcmpi( mask.MaskType,'1D') || strcmpi( mask.MaskType,'1Dpixel') 
                view_fx = fx;
                valid_fx = abs(view_fx)<2;
                view_fx = view_fx(valid_fx);
                view_fy = zeros(size(view_fx));
                absSpectrum = abs(Spectrum);
                view_Spectrum = absSpectrum(valid_fx);
                view_Spectrum = 2* view_Spectrum/max(view_Spectrum(:));
                quiver3(view_fx,view_fy,zeros(size(view_fx)),zeros(size(view_fx)),zeros(size(view_fx)),view_Spectrum,0);
                
            elseif strcmpi( mask.MaskType,'2D') || strcmpi( mask.MaskType,'2Dpixel') 
                [view_fx,view_fy] = meshgrid(fx,fy);
                
                valid_fx = abs(view_fx)<2;
                valid_fy = abs(view_fy)<2;
                validfxy = valid_fx & valid_fy;
                
                view_fx2 = view_fx(validfxy);
                view_fy2 = view_fy(validfxy);
                absSpectrum = abs(Spectrum);
                view_Spectrum = absSpectrum(validfxy);
                view_Spectrum = 2* view_Spectrum/max(view_Spectrum(:));
                
                quiver3(view_fx2,view_fy2,zeros(size(view_fx2)),zeros(size(view_fx2)),zeros(size(view_fx2)),view_Spectrum,0);
            else
                
            end
            rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'EdgeColor','r','LineWidth',1.5);
            axis equal
            box on
            xlabel('Frequency X');
            ylabel('Frequency Y');
            zlabel('Amplitude');
            title('Mask Spectrum');
            set(gca,'fontsize',13);
            set(gca,'FontName','Times New Roman');
            
        end
        
        % 绘制投影物镜光瞳图
        function PlotPupilPattern(mask,source,projector)
            %% 绘制衍射光瞳图像
            %******************************************
            %Pupil Filter进一步的支持和完善
            %******************************************
            %             Wavelength=source.Wavelength;
            %             NA=projector.NA;
            %             Pitch_X=mask.Period_X;
            %             Pitch_Y=mask.Period_Y;
            %             normPitch_X=Pitch_X/(Wavelength/NA);
            %             normPitch_Y=Pitch_Y/(Wavelength/NA);
            
            %%%%%%&&&&&%%%%%%
            [~,Source_Matrix] = Calc_SourceSimple( source );
            %%%%%%&&&&&%%%%%%
            
            %Nf,Ng：归一化频率的倍数，对应最大衍射级次
            %Nf→直接计算出来的像面点数：决定直接算出来的间隔→与仿真精度有关
            %Nf有最小值：要求f的范围至少为-2~2
            % normPitch：归一化的掩模周期
            % (1/normPitch)：衍射频谱的采样频率

            if strcmpi(mask.MaskType,'2D')||strcmpi(mask.MaskType,'2DPixel')
                N=5;
            elseif strcmpi(mask.MaskType,'1D')||strcmpi(mask.MaskType,'1DPixel')
                N=1;
            else
                N=1;
            end
            %             f=(1/normPitch_X)  *  (  -(Nf-1)/2:  1:   (Nf-1)/2   );
            %             f=(1/normPitch_Y)  *  (  -(Ng-1)/2:  1:   (Ng-1)/2   );
            
            %             masktransformer.fRange=[ f(1),f(end) ];
            %             masktransformer.gRange=[ f(1),f(end)  ];
            
            %%%%%%&&&&&%%%%%%
            %             spectrum = Mask2Spectrum(  masktransformer  );
            [ spectrum, f, g] = mask.CalculateMaskSpectrum(projector ,source); % 只取函数返回值的第一项
            %%%%%%%%%%%%%%%%%
            
            view_size = 51;
            view_fx = linspace(-1,1,view_size);
            view_gy = view_fx;
            [ source_fx,source_gy ] = meshgrid(view_fx,view_gy);
            [ source_fx2,source_gy2 ] = meshgrid(linspace(-1,1,size(Source_Matrix,1)));
            Source_temp = interp2(source_fx2,source_gy2,Source_Matrix,source_fx,source_gy,'spline');
            diffractionPattern = zeros(view_size);
            
            
            % absSpectrum = abs(Spectrum);
            % quiver3(view_fx,view_fy,zeros(size(view_fx)),zeros(size(view_fx)),zeros(size(view_fx)),absSpectrum);
            % rectangle('Position',[-1,-1,2,2],'Curvature',[1,1])
            % axis([-2,2,-2,2,0,0.15]);
            
            [ calc_fx,calc_gy ] = meshgrid(f,g);
            
            valid_fx = abs(calc_fx)<2;
            valid_gy = abs(calc_gy)<2;
            validfxy = valid_fx & valid_gy;
            
            view_fx2 = calc_fx(validfxy);
            view_gy2 = calc_gy(validfxy);
            if size(spectrum,1) == 1
                spectrum = repmat(spectrum,N,1);
            end
            view_Spectrum = spectrum(validfxy);
            
            view_fx2 = reshape(view_fx2,[],1);
            view_gy2 = reshape(view_gy2,[],1);
            view_Spectrum = reshape(view_Spectrum,[],1);
            
            for i = 1:length(view_fx2)
                temp_fx = source_fx + view_fx2(i);
                temp_fy = source_gy + view_gy2(i);
                
                tempdata = interp2(temp_fx,temp_fy,Source_temp ,source_fx,source_gy)* abs(view_Spectrum(i));
                tempdata(isnan(tempdata)) = 0;
                diffractionPattern = diffractionPattern + tempdata;
            end
            %             r2 = source_fx.^2+source_gy.^2;
            %             diffractionPattern(r2>1) = nan;
            pcolor(view_fx,view_gy,abs(diffractionPattern)/max(diffractionPattern(:)));
            set(gca,'YDir','normal');
            shading interp;
            colormap hot;
            rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'EdgeColor','r','LineWidth',1.5);
            set(gca,'FontName','Times New Roman');
            set(gca,'fontsize',14);
            colorbar;
            axis square;
            axis xy;
            xlabel('Pupil X');
            ylabel('Pupil Y');
            
            
        end
        
        % 绘制掩模图
        function PlotMask(mk,varargin)
            % 绘制掩模图像
            %Draw Mask
            %MATLAB:mathmatics->interpolation and computational geometry
            %绘图： fill (或patch)
            %
            
            switch length(varargin)
                case 0
                    drawSet.CHSV_H=0.8333;
                    drawSet.fontname='Times New Roman';
                    drawSet.fontsize=13;
                    drawSet.title='Mask in one period'; %'Mask in one period';
                    drawSet.xlabel='X / nm';
                    drawSet.ylabel='Y / nm';
                otherwise
                    drawSet = varargin{1};
            end
            
            
            Mask_Start_X=-mk.Period_X/2;
            Mask_Start_Y=-mk.Period_Y /2;
            maskEdge_X=[Mask_Start_X,   Mask_Start_X+mk.Period_X, ...
                Mask_Start_X+mk.Period_X , Mask_Start_X ];
            maskEdge_Y=[Mask_Start_Y,   Mask_Start_Y,...
                Mask_Start_Y+mk.Period_Y,     Mask_Start_Y+mk.Period_Y];
            
            CHSV_H = drawSet.CHSV_H; %色调 Hue~1
            CHSV_S = 1-abs(mod(mk.Background_Phase,2*pi)/pi -1); %饱和度: 0~1
            if mk.Background_Phase	> eps
                CHSV_V = 0.5;
            else
                CHSV_V = mk.Background_Transmissivity; %强度0~1
            end
            CBg = hsv2rgb([CHSV_H,CHSV_S,CHSV_V]);
            fill(maskEdge_X,maskEdge_Y,CBg);
            
            set(gca,'fontsize',drawSet.fontsize);
            set(gca,'fontName',drawSet.fontname);
            hold on;
            
            set(gca,'XLim',[-mk.Period_X/2 mk.Period_X/2 ] );
            set(gca,'YLim',[-mk.Period_Y/2 mk.Period_Y/2 ] );
            %背景完毕
            if strcmpi(mk.MaskType,'spectrum')
                error('Unsupported Mask Type !');
            elseif strcmpi(mk.MaskType,'2DPixel')
                [r,c]=size(mk.Feature);
                mx=linspace(-mk.Period_X/2,mk.Period_X/2,c);
                my=linspace(-mk.Period_Y/2,mk.Period_Y/2,r);
            
                mapHSV = zeros(mk.Ng,mk.Nf,3);
                mapHSV(:,:,1) = drawSet.CHSV_H.*ones(size(mk.Feature)); %色调 Hue~1
                mapHSV(:,:,2) = 1-abs(mod(angle(mk.Feature),2*pi)/pi -1); %饱和度: 0~1
                tempV = abs(mk.Feature);
                %                 tempV( logical((tempV>0).*(tempV<1)) ) = 0.5;%调试掩模需要注释
                mapHSV(:,:,3) = tempV;
                
                mapRGB = hsv2rgb(mapHSV);
                
                imagesc(mx,my,mapRGB);
                %                 colormap gray;
                %axis xy;
            else
                FeatureNum=size(mk.Feature,2);%特征图形个数
                for ii=1:FeatureNum
                    if strcmp(mk.Feature(ii).ShapeType,'r')
                        xv=[ mk.Feature(ii). BoundaryVertexX(1) ,  mk.Feature(ii). BoundaryVertexX(2)  ,...
                            mk.Feature(ii). BoundaryVertexX(2) ,  mk.Feature(ii). BoundaryVertexX(1) ];
                        yv=[ mk.Feature(ii). BoundaryVertexY(1) ,  mk.Feature(ii). BoundaryVertexY(1)  ,...
                            mk.Feature(ii). BoundaryVertexY(2) ,  mk.Feature(ii). BoundaryVertexY(2)  ];
                        
                        CHSV_H=drawSet.CHSV_H; %色调 Hue~1
                        CHSV_V=mk.Feature(ii).ComplexAm(1); %强度0~1
                        CHSV_S= 1- abs(mod( mk.Feature(ii).ComplexAm(2),2*pi )/pi-1) ; %饱和度: 0~1
                        Cm=hsv2rgb([CHSV_H,CHSV_S,CHSV_V]);
                        
                        fill(xv,yv,Cm,'EdgeColor',Cm);%用白色填充
                    elseif  strcmp(mk.Feature(ii).ShapeType,'t')
                        xv= mk.Feature(ii). BoundaryVertexX ;
                        yv= mk.Feature(ii). BoundaryVertexY  ;
                        
                        CHSV_H=drawSet.CHSV_H; %色调 Hue~1
                        CHSV_V=mk.Feature(ii).ComplexAm(1); %强度0~1
                        CHSV_S=1- abs(mod( mk.Feature(ii).ComplexAm(2),2*pi )/pi-1) ; %饱和度: 0~1
                        Cm=hsv2rgb([CHSV_H,CHSV_S,CHSV_V]);
                        
                        fill(xv,yv,Cm,'EdgeColor',Cm);
                    elseif  strcmp(mk.Feature(ii).ShapeType,'p')
                        xv= mk.Feature(ii). BoundaryVertexX ;
                        yv= mk.Feature(ii). BoundaryVertexY  ;
                        
                        CHSV_H=drawSet.CHSV_H; %色调 Hue~1
                        CHSV_V=mk.Feature(ii).ComplexAm(1); %强度0~1
                        CHSV_S=1- abs(mod( mk.Feature(ii).ComplexAm(2),2*pi )/pi-1) ; %饱和度: 0~1
                        Cm=hsv2rgb([CHSV_H,CHSV_S,CHSV_V]);
                        
                        fill(xv,yv,Cm,'EdgeColor',Cm);
                    elseif  strcmp(mk.Feature(ii).ShapeType,'c')
                        x1= mk.Feature(ii). BoundaryVertexX(1) ;
                        y1= mk.Feature(ii). BoundaryVertexX(2) ;
                        r= mk.Feature(ii). BoundaryVertexY  ;
                        
                        CHSV_H=drawSet.CHSV_H; %色调 Hue~1
                        CHSV_V=mk.Feature(ii).ComplexAm(1); %强度0~1
                        CHSV_S=1- abs(mod( mk.Feature(ii).ComplexAm(2),2*pi )/pi-1) ; %饱和度: 0~1
                        Cm=hsv2rgb([CHSV_H,CHSV_S,CHSV_V]);
                        
                        rectangle('Position',[x1-r,y1-r,2*r,2*r],'Curvature',[1,1],'FaceColor', Cm,'EdgeColor',Cm);
                    end
                end
            end
            % 画Cutline
            try
                if  ~isempty( mk.CutLine )
                    for i = 1:length(mk.CutLine)
                        line(mk.CutLine(i).X,mk.CutLine(i).Y,'Color','red','LineWidth',1);
                        plot(mk.CutLine(i).X,mk.CutLine(i).Y,'*','Color','red','MarkerSize',10);
                    end
                end
            catch
                
            end
            
            xlabel(drawSet.xlabel);
            ylabel(drawSet.ylabel);
            title(drawSet.title);
            set(gca,'FontName','Times New Roman');
            set(gca,'fontsize',14);
            axis equal;
            axis tight;
            axis([Mask_Start_X,-Mask_Start_X,Mask_Start_Y,-Mask_Start_Y])
        end
        
        % 照明光瞳绘图
        function PlotSource(source,view_size)
            flag_interp = false;
            switch nargin
                case 2
                    if strcmpi(source.Shape,'pixel')|| strcmpi(source.Shape,'point')
                        flag_interp = true;
                    else
                        source.PntNum = view_size;
                    end
            end
            
            [~,sourceMatrix,fx,fy]= source.Calc_SourceSimple();%
            
            if flag_interp
                view_fx = linspace( -1, 1, view_size);
                view_fy = view_fx;
                [ view_fxm, view_fym ] = meshgrid(view_fx, view_fy);
                view_sourceMatrix = interp2( fx, fy, sourceMatrix, view_fxm, view_fym);
            else
                view_fx = linspace( -1, 1, source.PntNum);
                view_fy = view_fx;
                [ view_fxm, view_fym ] = meshgrid(view_fx, view_fy);
                view_sourceMatrix = sourceMatrix;
            end
            

            r2m = view_fxm.^2 +view_fym.^2;
            view_sourceMatrix(r2m>1+eps) = 0;
            
             % 绘制光源
            imagesc( view_fx, view_fy, view_sourceMatrix );
            set(gca,'YDir','normal');
            shading flat;
            colormap hot;
            caxis([0,1]);
            colorbar;
            axis xy;
            axis square;
            xlabel('Pupil X');
            ylabel('Pupil Y');
            rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'EdgeColor','r','LineWidth',1.5);
            set(gca,'FontName','Times New Roman');
            set(gca,'fontsize',14);
            axis([-1,1,-1,1]);
            title('Illuminiation Source');
        end
 
        function PlotSource_MMA_PSF(source,view_size)

            srTemp = source;
            srTemp.MMA_Number = 1;
            srTemp.MMA_Coordinate = 'polar';
            srTemp.MMA_Rho = 0;
            srTemp.MMA_Theta = 0;
            srTemp.MMA_SourceConstructType = 'full';
            
            switch nargin
                case 2
                    srTemp.PntNum = view_size;
            end
            
            [~,temp,~,~] = srTemp.Calc_SourceAll;
            fx = linspace( -1, 1, source.PntNum);
            fy = fx;
            imagesc( fx, fy, temp);
            rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'EdgeColor','r','LineWidth',1.5);
            set(gca,'YDir','normal');
            shading interp;
            colormap hot;
            colorbar;
            axis xy; axis square;
            xlabel('Pupil X'); ylabel('Pupil Y');
            title('MMA PSF');
            set(gca,'FontName','Times New Roman');
            set(gca,'fontsize',14);
        end
        
        function PlotSource_MMA_Center(source)
             
            switch lower(source.MMA_Coordinate)
                case 'polar'
                    [ mmax, mmay ] = pol2cart(source.MMA_Theta,source.MMA_Rho);
                case 'cartesian'
                    mmax = source.MMA_X;
                    mmay = source.MMA_Y;
                otherwise
                    error('Unsupported Coordinate');
            end
            switch lower(source.MMA_SourceConstructType)
                case  'quarter'
                    mmaxFull = [ mmax,-mmax,-mmax, mmax];
                    mmayFull = [ mmay,mmay,-mmay, -mmay];
                case 'full'
                    mmaxFull = mmax;
                    mmayFull = mmay;
            end
            % 调试代码
            %                 xx = linspace(-1,1,101);
            %                 [ xx, yy ] = meshgrid( xx );
            %                 pcolor(xx,yy,Source_Partial);hold
            
            %     Source_Matrix(core_Xr(i),core_Yr(i)) = 1;
            
            scatter(mmaxFull,mmayFull,4,'filled');
            rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'EdgeColor','r','LineWidth',1.5);
            axis xy;
            axis square;
            xlabel('Pupil X');
            ylabel('Pupil Y');
            title('MMA Spot Center');
            set(gca,'FontName','Times New Roman');
            set(gca,'fontsize',14);
            box on;
            
        end
        
        % 绘制偏振分布
        function PlotPolarization(source)
           view_number_Polarization = 12;
           view_size = 301;

            %光源绘图
            [~,sourceMatrix, fx, fy]= source.Calc_SourceSimple();%
%             [ fx, fy ] = meshgrid(linspace(-1,1, source.PntNum));
            [ view_fx, view_fy ] = meshgrid(linspace( -1, 1, view_size));
            view_sourceMatrix = interp2( fx, fy, sourceMatrix, view_fx, view_fy);

            contour( view_fx, view_fy, view_sourceMatrix);
%             colormap hot;
             % 偏振数据计算
             [Px,Py] = meshgrid(linspace(-1,1,view_number_Polarization));
            Px = reshape(Px,1,[]);
            Py = reshape(Py,1,[]);
            [ theta, rho ] = cart2pol(Px, Py);
            
            Px(rho>=1) = [];
            Py(rho>=1) = [];
            theta(rho>=1) = [];
            rho(rho>=1) = [];
            
             [ PolarizedX, PolarizedY ]  = source.Calc_PolarizationMap(theta,rho);
            
            if length(PolarizedX)==1
                PolarizedX = PolarizedX*ones(size(Px));
                PolarizedY = PolarizedY*ones(size(Py));
            end
            % 偏振绘图
            hold on
            colorblue = [0,113,188]/255;
            colororange = [216,82,24]/255;
            if strcmpi( source.PolarizationType, 'unpolarized' )
                quiver(Px,Py,PolarizedX,PolarizedY,0.25,'Color',colorblue);%% 蓝色
                quiver(Px,Py,-PolarizedX,-PolarizedY,0.25,'Color', colorblue);
                quiver(Px,Py,-PolarizedX,PolarizedY,0.25,'Color',colororange);%橙色
                quiver(Px,Py,PolarizedX,-PolarizedY,0.25,'Color',colororange);

            else
                quiver(Px,Py,PolarizedX,PolarizedY,0.25,'Color', colorblue);
                quiver(Px,Py,-PolarizedX,-PolarizedY,0.25,'Color', colorblue);
                
%                                 quiver(Px,Py,PolarizedX,PolarizedY,0.25);
%                                 quiver(Px,Py,-PolarizedX,-PolarizedY,0.25);
            end
            
            rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'EdgeColor','r','LineWidth',1.5);
            axis square
            xlim([-1 1]);
            ylim([-1 1]);
            xlabel('Pupil X');
            ylabel('Pupil Y');
            set(gca,'FontName','Times New Roman');
            set(gca,'fontsize',14);
            title('Source Polarization');
            
        end
        
        % 绘制投影物镜波像差
        function PlotAberration(projector)
              view_size = 301;
              
              view_x = linspace(-1,1,view_size);
              view_y = view_x;
              [ view_x, view_y ]=meshgrid(view_x,view_y);
              [ theta, rho ] = cart2pol(view_x, view_y);
               pupil = projector.CalculateAberration( rho, theta, 0 );
               pupil(rho>1) = nan;
               pcolor( view_x, view_y, pupil );
               axis square; xlim([-1 1]); ylim([-1 1]); shading flat; colorbar;
               xlabel('Pupil X');ylabel('Pupil Y');title('Pupil Map');
               set(gca,'FontName','Times New Roman');
               set(gca,'fontsize',14);
        end
        
        % 绘制空间像
        function PlotAerialImage(lithoImage)
            if strcmpi( lithoImage.ImageType,'1D')
                if length( lithoImage.ImageZ)<2
                    plot( lithoImage.ImageX, lithoImage.Intensity,'LineWidth',1 );
                    xlabel('X/nm');
                    ylabel('Intensity');
                else
                    imagesc( lithoImage.ImageX, lithoImage.ImageZ, lithoImage.Intensity );
                    set(gca,'YDir','normal');
                    colorbar;
                    xlabel('X/nm');
                    ylabel('Z/nm');
                end
            elseif strcmpi( lithoImage.ImageType,'2D')
                if length( lithoImage.ImageZ ) >1
                    gridNum = 1001;
                    ai_xs = lithoImage.ImageX;
                    ai_ys = lithoImage.ImageY;
                    [ai_xm,ai_ym] = meshgrid(ai_xs, ai_ys);
                    view_xs = linspace(ai_xs(1),ai_xs(end),gridNum);
                    view_ys = linspace(ai_ys(1),ai_ys(end),gridNum);
                    [ view_xm, view_ym] = meshgrid(view_xs, view_ys);
                    view_intensity = interp2(ai_xm,ai_ym,lithoImage.Intensity,view_xm,view_ym);
                    imagesc( view_xs, view_ys,view_intensity );
                    %                     imagesc( lithoImage.ImageX, lithoImage.ImageY, lithoImage.Intensity );
                    set(gca,'YDir','normal');
                    colorbar;
                    xlabel('X/nm');
                    ylabel('Y/nm');

                else
                    gridNum = 1001;
                    ai_xs = lithoImage.ImageX;
                    ai_ys = lithoImage.ImageY;
                    [ai_xm,ai_ym] = meshgrid(ai_xs, ai_ys);
                    view_xs = linspace(ai_xs(1),ai_xs(end),gridNum);
                    view_ys = linspace(ai_ys(1),ai_ys(end),gridNum);
                    [ view_xm, view_ym] = meshgrid(view_xs, view_ys);
                    view_intensity = interp2(ai_xm,ai_ym,lithoImage.Intensity,view_xm,view_ym);
                    imagesc( view_xs, view_ys,view_intensity );
                    %                     imagesc( lithoImage.ImageX, lithoImage.ImageY, lithoImage.Intensity );
                    set(gca,'YDir','normal');
                    colorbar;
                    xlabel('X/nm');
                    ylabel('Y/nm');
                end
            end
            title('Aerial Image')
            set(gca,'FontName','Times New Roman');
            set(gca,'fontsize',14);
            
        end
        
        % 绘制光刻胶像
        function PlotResistImage(lithoImage)
            if strcmpi( lithoImage.ImageType,'1D')
                if length( lithoImage.ImageZ)<2
                    plot( lithoImage.ImageX, lithoImage.Intensity,'LineWidth',1 );
                    xlabel('X/nm');
                    ylabel('Intensity');
                else
                    imagesc( lithoImage.ImageX, lithoImage.ImageZ, lithoImage.Intensity );
                    set(gca,'YDir','normal');
                    colorbar;
                    xlabel('X/nm');
                    ylabel('Z/nm');
                    axis equal
                    axis tight
                end
            elseif strcmpi( lithoImage.ImageType,'2D')
                if length( lithoImage.ImageZ ) >1
                    gridNum = 201;
                    ai_xs = lithoImage.ImageX;
                    ai_ys = lithoImage.ImageY;
                    ai_zs = lithoImage.ImageZ;
                    
                    [ai_xm,ai_ym,ai_zm] = meshgrid(ai_xs, ai_ys, ai_zs);
                    view_xs = linspace(ai_xs(1),ai_xs(end),gridNum);
                    view_ys = linspace(ai_ys(1),ai_ys(end),gridNum);
                    view_zs = linspace(ai_zs(1),ai_zs(end),gridNum);
                    [ view_xm, view_ym,view_zm] = meshgrid(view_xs, view_ys,view_zs);
                    
                    xslice = 0;
                    yslice = 0;
                    zslice = ai_zs(end);
                    view_intensity = interp3(ai_xm,ai_ym,ai_zm,lithoImage.Intensity,view_xm,view_ym,view_zm,'cubic');
                    slice(view_xm,view_ym,view_zm,view_intensity,xslice,yslice,zslice);
                    shading flat
                    colorbar;
                    xlabel('X/nm');
                    ylabel('Y/nm');
                    zlabel('Z/nm');
                    title('Resist Image')
                    set(gca,'FontName','Times New Roman');
                    set(gca,'fontsize',14);
                else
                    imagesc( lithoImage.ImageX, lithoImage.ImageY, lithoImage.Intensity );
                    set(gca,'YDir','normal');
                    colorbar;
                    xlabel('X/nm');
                    ylabel('Y/nm');
                    axis equal
                    axis tight
                end
            end
            title('Resist Image')
            set(gca,'FontName','Times New Roman');
            set(gca,'fontsize',14);
        end
        
        % 绘制光刻胶轮廓
        function PlotResistContour(lithoImage,nLevel)
            if strcmpi( lithoImage.ImageType,'2D')
                if length( lithoImage.Focus_Range)<2
                    contour ( lithoImage.ImageX, lithoImage.ImageY, lithoImage.Intensity,nLevel,'LineWidth',1.5);
                    
                    %                     imagesc( lithoImage.ImageX, lithoImage.ImageY, lithoImage.Intensity );
                    xlabel('X/nm');
                    ylabel('Y/nm');
                else
                    error( 'Unsupported Data' );
                end
            else
                error( 'This function only support 2D mask' );
            end
            box on
            title('Resist Contour');
            set(gca,'FontName','Times New Roman');
            set(gca,'fontsize',14);
        end
        
        % 绘制光刻胶堆栈
        function PlotFilmStack(filmStack)
            for iARC = 1:length(filmStack.Layers)
                
            end
        end
    end
end