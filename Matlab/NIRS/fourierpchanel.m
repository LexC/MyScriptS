function fourierpchanel(run,dOD,figtitle,figdir,override)
    close all
    newfigdir = fullfile(figdir,strcat(figtitle,'.png'));
    
    if exist(newfigdir,'file')==2 && override == false
        okflag = false;
    else
        okflag = true;
    end
    
    if okflag
        checkdir(figdir)

        dODw1 = dOD(:,1:size(dOD,2)/2);
        dODw2 = dOD(:,size(dOD,2)/2+1:end);

        f = figure(1);
        f.Position=[0,0,1920,1080];

        sgtitle(figtitle)

        for k = 1:64
            ax(k) = subplot(8,8,k);
        end
        frec = 1/(run.t(2,1) - run.t(1,1)); % frequencia de mostragem

        for i = 1:64
            subplot(ax(i));
            if i<=size(dODw1,2)
                [pxx,wx1] = pwelch(dODw1(:,i),[],[],[],frec);
                plot(wx1,10*log10(pxx),'- b')
                hold on
                [pxx,wx1] = pwelch(dODw2(:,i),[],[],[],frec);
                plot(wx1,10*log10(pxx),'- k')
                hold off
            else
                plot([0,0],[0,0],'- b')
            end
            title(sprintf('C:%d',i));
            grid
            xlim([0,2])
        end


        l = legend(strcat(num2str(run.SD.Lambda(1)),' nm'),strcat(num2str(run.SD.Lambda(2)),' nm'),...
            'FontSize',12);
        l.Position(1:2)=[0.851,0.925];
        saveas(f,newfigdir,'png')
        close all
    end
end