function chanels_corrections_plots(dOD,new_dOD,figdir,figtitle,override)
    close all
    checkdir(figdir)
    
    aux = replace(figtitle,' ','_');
    newfigdir = fullfile(figdir,[aux,'_',num2str(ceil(size(dOD,2)/16)),'.png']);
    
    if exist(newfigdir,'file')==2 && override == false
        okflag = false;
    else
        okflag = true;
    end
    
    if okflag
        for i=1:ceil(size(dOD,2)/16)
            f = figure(i);
            f.Position=[0,0,1920,1080];

            sgtitle(figtitle)
            for j = 1:16
                ax(j) = subplot(4,4,j);
            end

            for j = 1:16
                u = j+(i-1)*16;
                subplot(ax(j));
                if u<=size(dOD,2)
                    plot(dOD(:,u),'k')
                    hold on
                    plot(new_dOD(:,u),'r')
                    xlim([0 size(dOD,1)])
                else
                    plot([0,0],[0,0],'- b')
                end
                hold off
                title(sprintf('C:%d',u));
                grid
            end
            aux = replace(figtitle,' ','_');
            newfigdir = fullfile(figdir,[aux,'_',num2str(i),'.png']);

            l = legend('Original','New',...
                'FontSize',12);
            l.Position(1:2)=[0.851,0.925];

            saveas(gcf,newfigdir,'png')

        end
    end
    close all
end