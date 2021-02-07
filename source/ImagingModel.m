classdef ImagingModel < handle
    
    properties
        Numerics
        Source
        Mask
        Projector
        LithoImageProtype
        FilmStack
        Metrology
        Receipe
    end
    
    methods        
        function model = ImagingModel()
            model.Numerics = Numerics;
            model.Source = Source;
            model.Mask = Mask.CreateMask('line_space');
            model.Projector = ProjectionObjective;
            model.FilmStack = FilmStack;
            model.LithoImageProtype = ImageData;
            model.Receipe = Receipe;
        end
        
        function  ali = CalculateAerialImage( model )            
            sr = model.Source;
            mk = model.Mask;
            po = model.Projector;
            rp = model.Receipe;
            nm = model.Numerics;
            lip = model.LithoImageProtype;
            if strcmpi( model.Mask.MaskType, '1D' )
                ali = Calculate1DAerialImage( sr, mk, po, lip, rp, nm );
            elseif strcmpi( model.Mask.MaskType, '2D' ) || strcmpi( model.Mask.MaskType, '2DPixel' )
                if strcmpi(nm.ImageCalculationMethod,'abbe')
                    ali = Calculate2DAerialImage( sr, mk, po, lip, rp, nm );
                elseif strcmpi(nm.ImageCalculationMethod,'hopkins')
                    ali = lip;
                    [ TCCMatrix_Stacked,FG_ValidSize ] = Calculate2DTCCMatrix(sr,mk,po,rp,nm);
                    TCCMatrix_Kernel = DecomposeTCC_SOCS( TCCMatrix_Stacked,FG_ValidSize,nm);
                    ali.Intensity = CalculateAerialImage_SOCS(mk,TCCMatrix_Kernel,sr,po,nm);
                    ali.ImageX = linspace(-mk.Period_X/2, mk.Period_X/2, nm.SampleNumber_Wafer_X);
                    ali.ImageY = linspace(-mk.Period_Y/2, mk.Period_Y/2, nm.SampleNumber_Wafer_Y);
                    ali.ImageZ = rp.FocusRange;
                    ali.ImageType = '2d';
                    ali.SimulationType = 'aerial';
                else
                    error('Unsupported Calculation Method');
                end
            else
                error('Unsupported Mask');
            end
            
            if nm.Normalization_Intensity
                ni = CalculateNormalImage( sr, mk, po, rp, nm);
                ali.Intensity = ali.Intensity / ni;
            end
                
        end
        
        function rli = CalculateResistImage( model )
            sr = model.Source;
            mk = model.Mask;
            po = model.Projector;
            fs = model.FilmStack;
            rp = model.Receipe;
            nm = model.Numerics;
            lip = model.LithoImageProtype;
            if strcmpi( model.Mask.MaskType, '1D' )
                rli = Calculate1DResistImage( sr, mk, po, fs, lip, rp, nm);
            elseif strcmpi( model.Mask.MaskType, '2D' ) || strcmpi( model.Mask.MaskType, '2DPixel' )
                if strcmpi(nm.ImageCalculationMethod,'abbe')
                rli = Calculate2DResistImage( sr, mk, po, fs, lip, rp, nm);
                elseif strcmpi(nm.ImageCalculationMethod,'hopkins')
                    rli = lip;
                    [ TCCMatrix_Stacked,FG_ValidSize ] = Calculate2DResistTCCMatrix(sr,mk,po,fs,rp,nm);
                    TCCMatrix_Kernel = DecomposeTCC_SOCS( TCCMatrix_Stacked,FG_ValidSize,nm);
                    rli.Intensity = CalculateAerialImage_SOCS(mk,TCCMatrix_Kernel,sr,po,nm);
                    rli.ImageX = linspace(-mk.Period_X/2, mk.Period_X/2, nm.SampleNumber_Wafer_X);
                    rli.ImageY = linspace(-mk.Period_Y/2, mk.Period_Y/2, nm.SampleNumber_Wafer_Y);
                    rli.ImageZ = rp.FocusRange;
                    rli.ImageType = '2d';
                    rli.SimulationType = 'aerial';
                else
                   error('Unsupported Calculation Method'); 
                end
            else
                error('Unsupported Mask');
            end
            if nm.Normalization_Intensity
                ni = CalculateNormalImage( sr, mk, po, rp, nm);
                rli.Intensity = rli.Intensity / ni;
            end
        end
        
        function eli = CalculateExposedLatentImage(model,rli)
            if nargin == 1
                rli = model.CalculateResistImage();
            elseif nargin == 2
            
            else
                error('Unsupported Input');
            end
            fs = model.FilmStack;
            rp = model.Receipe;
            nm = model.Numerics;
            eli = CalculateExposedLatentImage(rli,fs,rp,nm);
        end
    end
    
end