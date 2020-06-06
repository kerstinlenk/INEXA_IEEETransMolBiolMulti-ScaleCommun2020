%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m: Graph theory analysis for simulation results from INEXA model
% (see publication with doi: 10.3389/fncom.2019.00092
% author: Barbara Genocchi & Michael T. Barros & Kerstin Lenk 
% date: 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT: 
% change the path where the results are stored $results_path$.
% maximum distance between astrocytes $d$,
% the number of simulation runs $run$

%% OUTPUT:
% Figures are saved in run1 for each distances, respectively
% Graph_analysis.mat with L, k and N_act is saved in run10 for each distances, respectively

clear all
close all


 d = [70, 80, 90, 100, 110, 120]; %change here

for link = 1:length(d)
    
    N_ast_run = [];
    
    N_neu_run = [];
    
    for run = 1: 10 %change here
        %change result_path
        results_path = ['C:\linkRadius_d_' num2str(d(link)) '_reducNeuCon0.55_c_0.02_y_0.70_run' num2str(run) '/'];
        cd(results_path)
        
        %% astro
        AstroConnections = readtable('AstrocyteConnections.csv');
        AstroConnections = AstroConnections(1:end,1:end-1);
        AstroConnections = table2array(AstroConnections);
        AstroData = readtable('AstroData_0.0200_0.7000_0.7000.csv');
        AstroData = table2array(AstroData(1:end,1:end-1));
        %needed for coloured graph
        AstroInfo = readtable('AstrocyteNetworkTopology.csv');
        AstroInfo = AstroInfo(1:end,1:end-1);
        AstroInfo = table2array(AstroInfo);
        
        tinMS = 5*60*1000; 
               
        if ~((isempty(AstroData))||((length(AstroData(:,1)))<14))
            activity = AstroData(14:14+length(AstroConnections)-1, 1:end-1);
        end
        
        %% N_act number of astro active at least once in the simulation
        size_act = size(activity);
        
        ActiveAstro = 0;
        find_astro = zeros(size_act(1),1);
        aa = [];
        
        for j = 1 : size_act(1)
            
            aa(j,1) = sum(activity(j,:));
        end
        
        astro_act = 0;
        for ll = 1: length(aa)
            b = isequal(aa(ll,1),0);
            
            if b == 0
                astro_act = astro_act +1;
            end
        end
        
        N_act_a = astro_act;
        
        N_ast_run(run,1) = N_act_a;
        
                
        if run ==1 %same topology no need to evaluate each run
            
        %% mean shortest path
        
        Adj_mat_ast = AstroConnections;
        
        Astro_graph = graph(Adj_mat_ast);
        
        L_astro_min = [];
        
        for bb = 1: length(AstroConnections)
            for cc = 1: length(AstroConnections)
                
                if bb~=cc %no self connections
                
                [L_path, min_path] = shortestpath(Astro_graph,bb,cc);
                
                L_astro_min = [L_astro_min,min_path];
                end
            end
        end
        
        L_astro_mean = mean(L_astro_min);
        L_astro_std = std(L_astro_min);
        
             
        %% degree from graph
        
        k_ast = degree(Astro_graph);
        k_ast_mean = mean(k_ast);
        k_ast_std = std(k_ast);
        
        %% coloured graph
        
        if  run ==1
        f1 = figure(1);
        
        colormap autumn
        deg = degree(Astro_graph);
        nSizes = 2*sqrt(deg-min(deg)+0.2);
        nColors = deg;
        plot(Astro_graph,'XData',AstroInfo(1,:),'YData',AstroInfo(2,:),'ZData',AstroInfo(3,:),'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.1,'NodeLabel',{},'LineWidth',2)
        colorbar
        title_astro_col = ['Graph-coloured-map-astro-d' num2str(d(link)) ''];
        title(title_astro_col)
        title_astro_fig_col = [title_astro_col, '.fig'];
        
        [M,I] = max(k_ast);
         savefig(f1, title_astro_fig_col)
         close(f1)
        end
        
        %normal figure
        
        f2 = figure(2);
        plot(Astro_graph)
        title_astro = ['Graph-map-astro-d' num2str(d(link)) ''];
        title(title_astro)
        
        title_astro_fig = [title_astro, '.fig'];
        savefig(f2, title_astro_fig)
        close(f2)
        end
               
        %% neurons
        NeuronConnections = readtable('BaseNetwork.csv');
        NeuronConnections = table2array(NeuronConnections(1:end,1:end-1));
            
        %% L mean shortest path - unweighted
        if run ==1
        NeuronConnections_unw = zeros(length(NeuronConnections),length(NeuronConnections));
        for xx = 1: length(NeuronConnections)
            for jj = 1 : length(NeuronConnections)
                if NeuronConnections(xx,jj)~=0  || NeuronConnections(jj,xx)~=0
                    NeuronConnections_unw(xx,jj) = 1;
                    NeuronConnections_unw(jj,xx) = 1;
                end
            end
        end
        
        Adj_mat_neu = NeuronConnections_unw;
        
        Neu_graph_unw = graph(Adj_mat_neu);
        
        L_neu_min = [];
        
        for bb = 1: length(NeuronConnections)
            for cc = 1: length(NeuronConnections)
                
                if bb~=cc %no self connections
                    [L_path, min_path] = shortestpath(Neu_graph_unw,bb,cc);

                    L_neu_min = [L_neu_min, min_path];
                end
            end
        end
        
        L_neu_mean = mean(L_neu_min);
        L_neu_std = std(L_neu_min);
        
        %% degree from graph
        
        k_neu = degree(Neu_graph_unw);
        k_neu_mean = mean(k_neu);
        k_neu_std = std(k_neu);

        %% N_act_neurons
        
        Neuron_activity = readtable('timestamps_0.0200_0.7000_0.7000.csv');
        Neuron_activity = table2array(Neuron_activity(2:end,1:end-1));
        
        [stamps NumberOfNeuron] = size(Neuron_activity);
        N_act_neu = 0;
        
        for tt = 1: NumberOfNeuron
            if Neuron_activity(1,tt) ~= 0
                N_act_neu = N_act_neu + 1;
            end
        end
        
        N_neu = N_act_neu;    
        
        %% coloured graph 
        f5 = figure(5);
        plot(Neu_graph_unw)
        title_neu = ['Graph-map-neurons-d' num2str(d(link)) ''];
        title(title_neu)
        
        title_neu_fig = [title_neu, '.fig'];
        savefig(f5, title_neu_fig)
        close(f5)
        end

    end
    
    N_ast_run_mean = mean(N_ast_run);
    N_ast_run_std = std(N_ast_run);
        
    save_file = [results_path, 'Graph_analysis.mat'];
    save(save_file, 'L_astro_mean', 'k_ast_mean', 'k_ast_std', 'N_ast_run', 'N_ast_run_mean', 'N_ast_run_std', 'N_neu', 'L_neu_mean', 'L_neu_std', 'k_neu_mean', 'k_neu', 'k_neu_std', 'k_ast_std', 'L_neu_std', 'L_astro_std');
end
