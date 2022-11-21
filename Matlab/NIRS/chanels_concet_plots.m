function chanels_concet_plots(conc,new_conc,figdir,figtitle,override)
    close all
    checkdir(figdir)
    
    i=ceil(size(conc,2)/16);
    newfigdir = replace(replace(figtitle,',',''),' ','_');
    newfigdir = fullfile(figdir,[newfigdir,'_',num2str(i),'.png']);
    
    if exist(newfigdir,'file')==2 && override == false
        okflag = false;
    else
        okflag = true;
    end
    
    if okflag
        for i=1:ceil(size(conc,2)/16)
            f=figure(i);
            f.Position=[0,0,1920,1080];
            sgtitle(figtitle)
            for j = 1:16
                ax(j) = subplot(4,4,j);
            end

            for j = 1:16
                u = j+(i-1)*16;
                subplot(ax(j));
                if u<=size(conc,2)
                    plot(conc(:,u,3),'--k')
                    hold on
                    plot(new_conc(:,u,3),'k')
                    plot(conc(:,u,1),'--r')
                    plot(conc(:,u,2),'--b')
                    plot(new_conc(:,u,1),'r')
                    plot(new_conc(:,u,2),'b')

                else
                    plot([0,0],[0,0],'- b')
                end
                hold off
                title(sprintf('C:%d',u));
                grid
                xlim([0 size(conc,1)])
            end

            newfigdir = replace(replace(figtitle,',',''),' ','_');
            newfigdir = fullfile(figdir,[newfigdir,'_',num2str(i),'.png']);

            l = legend('Original','New',...
                'FontSize',12);
            l.Position(1:2)=[0.851,0.925];
            saveas(gcf,newfigdir,'png')
        end
    end
    close all
end