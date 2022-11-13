function [ stru ] = graphparameters( map )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %% -------> Variables <-------
    u=0;
    g{19}=[];adjmap=g;k=g;K=g;StdK=g;p=g;d=g;D=g;StdD=g;c=g;C=g;StdC=g;
    %%
    
    for i=0.05:0.05:0.95
        fprintf('...')
        %% -------> Making Adjacency Matrices <-------
        u=u+1;
        aux=map;
        aux(map>=i)=1;
        aux(map<i)=0;
        adjmap{u}=aux;

        
        %% -------> Graphs Parameters <-------
        
        % graph(g{u})
        g{u}=graph(aux,'omitselfloops');
        
        %% Degree(k{u})
        k{u}=degree(g{u});
        
        % Average Degree(K{u})
        K{u}=mean(k{u});
        StdK{u}=std(k{u});
       

        %% Shortest Path(p{u}) and distance(d{u})
        p{u}{size(aux,1),size(aux,2)}=[];
        d{u}=zeros(size(aux));
        for x=1:size(aux,1)
            for y=1:size(aux,2)
                if x>y
                    [p{u}{x,y},d{u}(x,y)]=shortestpath(g{u},x,y);
                end
            end
        end
        
        % Average Path Length (D{u})
        D{u}=mean(d{u}(d{u}>0));
        StdD{u}=std(d{u}(d{u}>0));
        %% Clustering Coefficient(c{u})
        
        c{u}=zeros(size(k{u}));
        
        for x1=1:size(k{u},1)
            % Defining neighbors of x1
            s1=aux(:,x1);
            t1=find(s1==1);
            t1(t1==x1)=[];
            
            if isempty(t1)~=1
                for x2=1:size(t1,1)

                    % Defining links between neighbors of x1
                    s2=aux(:,t1(x2));
                    t2=find(s2==1);
                    t2(t2==t1(x2))=[];
                    neibs=logical((ismember(t2,t1)-1)*(-1));
                    t2(neibs)=[];

                    % Computing pairs
                    if x2==1 
                        pairs=[t1(x2)*ones(size(t2)),t2];
                    elseif isempty(pairs)==1
                        pairs=[t1(x2)*ones(size(t2)),t2];
                    else
                        pairs=[pairs;t1(x2)*ones(size(t2)),t2];
                    end
                end

                % Cleaning repeated pairs
                pairs=sort(pairs,2);
                pairs=unique(pairs,'rows');
            else
                pairs=[];
            end
            
            if isempty(pairs)==1
                c{u}(x1,1)=0;
            else
                c{u}(x1,1)=2*size(pairs,1)/(k{u}(x1)*(k{u}(x1)-1));
            end

        end
        
        % Average Clustering Coefficient(C{u})
        C{u}=mean(c{u});
        StdC{u}=mean(c{u});
        
    end
    disp('.')
    stru.AdjMat=adjmap;
    stru.graph=g;
    stru.k=k;
    stru.K=K;
    stru.StdK=StdK;
    stru.p=p;
    stru.d=d;
    stru.D=D;
    stru.StdD=StdD;
    stru.c=c;
    stru.C=C;
    stru.StdC=StdC;
end

