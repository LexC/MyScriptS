rho=0.05:.05:.95;
x=10;y=15;z=2;

%% 
u=1;w=1;
for i=1:(x+y+z)
    %%
    if i<=x
        for j=1:size(CorrMat(i).graphs.K,2)
            k{i}(j)=CorrMat(i).graphs.K{j};
            c{i}(j)=CorrMat(i).graphs.C{j};
            p{i}(j)=CorrMat(i).graphs.D{j};
        end
        filename{i}='Run';
    elseif i<=x+y
        for j=1:size(CorrMatMean(w).graphs.K,2)
            k{i}(j)=CorrMatMean(w).graphs.K{j};
            c{i}(j)=CorrMatMean(w).graphs.C{j};
            p{i}(j)=CorrMatMean(w).graphs.D{j};
        end
        w=w+1;
        filename{i}='Subject';
    else
        for j=1:size(CorrMatGroupMean(u).graphs.K,2)
            k{i}(j)=CorrMatGroupMean(u).graphs.K{j};
            c{i}(j)=CorrMatGroupMean(u).graphs.C{j};
            p{i}(j)=CorrMatGroupMean(u).graphs.D{j};
        end
        u=u+1;
        filename{i}='Group';
    end

end

%%
u=3;
for i=1:size(k,2)
    %%
    switch u
        case 1
            figure (1)
                plot(rho,k{i},'-o')
                title(strcat(filename{i},' - \rho vs <K>'))
                xlabel('\rho')
                ylabel('<K>')
                grid on
                axis([0 1 0 120])
%                 hold on
            pause(1)
        case 2
            figure (2)
                plot(rho,p{i},'-o')
                title(strcat(filename{i},' - \rho vs <Shortest Path>'))
                xlabel('\rho')
                ylabel('<K>')
                grid on
                axis([0 0.6 1 5])
%                 hold on
            pause(1)
        case 3
            figure (3)
                plot(rho,c{i},'-o')
                title(strcat(filename{i},' - \rho vs <Clustering>'))
                xlabel('\rho')
                ylabel('<K>')
                grid on
                hold on
                axis([0 1 0 1])
            pause(1)
    end
end