function NIRSPreProcess(run,par,foldername,filename,override)

    checkdir(foldername)
    finalfilelob = [filename,'.lob'];
    finalfilemat = [filename,'.mat'];
    flagfile = [foldername,'flag.mat'];
    
    okflag = true;
    if exist(finalfilelob,'file')==2 && override == false
        okflag = false;
    elseif exist(flagfile,'file')==2
        load (flagfile,'flag')
        if flag == 1 && override == false
            okflag = false;
        end
            
    end
    
    if okflag
        try
            % Adding important Variables
            run.SD.MeasListAct = ones(size(run.d,2),1);
            run.SD.SCI = run.ComputeScalpCouplingIndex(par.SCIFrequency);
            run.SD.parameters = par;
            run = run.MarkBadChannels(par);

            % > Optical Density Conversion
            dOD = run.Convert2OD;
            run.SD.dOD_original = dOD;

            %pfolder = fullfile(foldername,'Plots');
            %fourierpchanel(run,dOD,'Fourier by channel',pfolder);

            % > Corrections
            % >> Spline

            dOD_spline = run.SplineCorrection(dOD, par.Spline, []);
            run.SD.dOD_spline = dOD_spline;

            %splotstitle = strcat('Spline',replace(num2str(par.Spline),'.',','));
            %chanels_corrections_plots(dOD,dOD_spline,pfolder,splotstitle)

            % >> Wavelet

            dOD_wavelet = run.WaveletCorrection(dOD_spline, par.Wavelet);
            run.SD.dOD_wavelet = dOD_wavelet;

            %wplotstitle = [splotstitle, ' Wavelet', replace(num2str(par.Wavelet),'.',',')];
            %chanels_corrections_plots(dOD_spline,dOD_wavelet,pfolder,wplotstitle)

            % > Concentration Conversion

            dc = run.Convert2Conc(dOD_wavelet,par.Convert2Conc);
            run.SD.dc_original = dc;

            % > Filters
            % >> Bandpass
            bp_filtered = run.BPFilter(dc,par.FrequencyRange); 

            bp_filtered = bp_filtered(par.BordersRemoval:end-par.BordersRemoval,:,:);
            run.SD.bp_filtered = bp_filtered;

            run.t = run.t(par.BordersRemoval:end-par.BordersRemoval,:,:);

            % >> PCA
            run.dc = bp_filtered;
            pca_filtered = run.PCAFilter(par);
            run.dc = pca_filtered;

            %figtitle = 'Post Bandpass, PC-1';
            %chanels_concet_plots(bp_filtered,pca_filtered,pfolder,figtitle)


            flag=false;
            runstruct = struct(run);

            save(finalfilelob,'run')
            save(finalfilemat,'runstruct')
            save(flagfile,'flag')
        catch
            flag=true;
            save(flagfile,'flag')
        end
    end
end
                
                



