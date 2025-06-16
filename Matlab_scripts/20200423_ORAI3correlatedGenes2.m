% Loads & organizes data

load('03.14.19_workspace_immgenRNAseq_full.mat')
load('03.14.19_workspace_haemosphereMouse_RNA_full.mat')

% Process Immgen
values = cell2mat(immgen_rna([2:end],[2:end]));
cells = immgen_rna(1,[2:end]);
genes = immgen_rna([2:end],1);

immgenMouse = {};
immgenMouse.values = values;
immgenMouse.cells = cells;
immgenMouse.genes = genes;

% Process Haemopedia Mouse
values = cell2mat(haemosphereMouse_rna([2:end],[2:end]));
cells = haemosphereMouse_rna(1,[2:end]);
genes = haemosphereMouse_rna([2:end],1);

haemopediaMouse = {};
haemopediaMouse.values = values;
haemopediaMouse.cells = cells;
haemopediaMouse.genes = genes;

clear immgen_rna haemosphereMouse_rna values cells genes

%% Acquires correlated genes
%
% Input:
%   celltypes               X-by-2 list of cells (columns represent
%                           different datasets)
%

dataSets = {'immgenMouse';
    'haemopediaMouse'};
celltypes = cellQ;
geneQ = 'ORAI3';
exprThresh = 10;
rhoPctle = 0.1;

for i = [1:length(dataSets(:,1))]
    dataSet = dataSets{i,1};
    eval(['values = ',dataSet,'.values;']);
    eval(['cells = ',dataSet,'.cells;']);
    eval(['genes = ',dataSet,'.genes;']);
    
    celltype = celltypes(:,i);
    exclude = find(cellfun('isempty',celltype));
    celltype(exclude,:) = [];

    % Converts cell-names to char. list
    for j = [1:length(cells)]
        cells{1,j} = char(string(cells{1,j}));
    end
    for j = [1:length(genes(:,1))]
        genes{j,1} = char(string(genes{j,1}));
    end
    
    % Acquires appropriate cells
    celltype = unique(celltype);
    include = [];
    for j = [1:length(celltype(:,1))]
        cell = celltype{j,1};
        loc = find(strcmp(lower(cell),lower(cells)));
        if not(isempty(loc))
            include = [include, loc];
        end
    end
    values = values(:,include);
    cells = cells(1,include);
    
    % Acquires 'geneQ'
    loc = find(strcmp(lower(genes),lower(geneQ)));
    if not(isempty(loc))
        val1 = values(loc,:);
    end
    values(loc,:) = [];
    genes(loc,:) = [];
    
    % Acquires correlated genes (& removes those under an expression threshold)
    genesOut = {};
    valOut = [];
    for j = [1:length(genes(:,1))]
        gene = genes{j,1};
        val2 = values(j,:);
        
        % Removes genes falling below a certain threshold
        meanVal = mean(val2);
        if meanVal < exprThresh
            continue
        end
        
        % Determines correlation
        rho = corr(val1',val2');
        
        % Appends information
        genesOut = [genesOut; gene];
        valOut = [valOut; rho];
    end
    
    temp = [genesOut, num2cell(valOut)];
    temp = sortrows(temp,2,'desc');
    topNum = ceil(length(genesOut)*rhoPctle);
    temp = temp([1:topNum],:);
    eval(['temp',char(string(i)),' = temp;']);
end

% Performs overlap analysis
allGenes = [temp1(:,1);temp2(:,1)];
uniqueGenes = unique(allGenes);
genesOut = {};
valOut = [];
for i = [1:length(uniqueGenes(:,1))]
    gene = uniqueGenes{i,1};
    loc1 = find(strcmp(lower(temp1(:,1)),lower(gene)));
    loc2 = find(strcmp(lower(temp2(:,1)),lower(gene)));
    if not(isempty(loc1)) & not(isempty(loc2))
        val1 = temp1{loc1(1),2};
        val2 = temp2{loc2(1),2};
        meanVal = mean([val1 val2]);
        
        genesOut = [genesOut; gene];
        valOut = [valOut; meanVal];
    end
end
output = [genesOut, num2cell(valOut)];
output = sortrows(output,2,'desc');

clear allGenes cell cells celltype celltypes dataSet dataSets exclude exprThresh
clear gene geneQ genes genesOut i include j loc loc1 loc2 meanVal rho
clear rhoPctle rhoThresh temp temp1 temp2 topNum uniqueGenes val1 val2
clear valOut values
    
%%

directory1 = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\05.20.18_ORAI3 Analyses\Data\2020.04.23_ORAI3correlations2\datasetCorrelations';
directory2 = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\05.20.18_ORAI3 Analyses\Data\2020.04.23_ORAI3correlations2\ARCHS4Correlations';

directories = {'directory1';
    'directory2'};


for i = [1:length(directories(:,1))]
    eval(['directory = ',directories{i,1},';']);
    ls = dir(directory);
    
    % Sifts through files
    pathwayOut = {};
    pathwayFamilyOut = {};
    for j = [1:length(ls(:,1))]
        filename = ls(j).name;
        pathwayFamily = strsplit(filename,'_');
        pathwayFamily = pathwayFamily{1};
        loc = strfind(filename,'.txt');
        if not(isempty(loc))
            filename = [directory,'\',filename];
        else
            continue
        end
        
        % Imports data
        data = readtable(filename);
        data = table2cell(data);
        
        % Filter out pathways with p-values < 0.05
        for k = [1:length(data(:,1))]
            if data{k,3} < 0.05
                pathwayOut = [pathwayOut; data(k,:)];
                pathwayFamilyOut = [pathwayFamilyOut; pathwayFamily];
            end
        end
    end
    
    pathwayOut = [pathwayFamilyOut, pathwayOut];
    pathwayOut = sortrows(pathwayOut,9,'desc');
    eval(['pathwayOut',char(string(i)),' = pathwayOut;']);
end

% Overlap analysis
allPathways = [pathwayOut1(:,2); pathwayOut2(:,2)];
uniquePathways = unique(allPathways);
pathOut = {};
famOut = {};
fracOut = {};
pOut = [];
fdrOut = {};
combinedOut = [];
genesOut = [];
for i = [1:length(uniquePathways(:,1))]
    pathway = uniquePathways{i,1};
    loc1 = find(strcmp(lower(pathwayOut1(:,2)),lower(pathway)));
    loc2 = find(strcmp(lower(pathwayOut2(:,2)),lower(pathway)));
    if not(isempty(loc1)) & not(isempty(loc2))
        val1 = pathwayOut1(loc1(1),:);
        val2 = pathwayOut2(loc2(1),:);
        
        path = pathway;
        fam = val1{1};
        frac = [val1{3},'|',val2{3}];
        pval = pfast([val1{4},val2{4}]);
%         pval = mean([val1{4},val2{4}]);
        fdr = [char(string(val1{5})),'|',char(string(val2{5}))];
        combined = mean([val1{9},val2{9}]);
        genes = [val1{10},'|',val2{10}];
        
        pathOut = [pathOut; pathway];
        famOut = [famOut; fam];
        fracOut = [fracOut; frac];
        pOut = [pOut; pval];
        fdrOut = [fdrOut; fdr];
        combinedOut = [combinedOut; combined];
        genesOut = [genesOut; {genes}];
    end
end
output = [pathOut,famOut,fracOut,num2cell(pOut),fdrOut,num2cell(combinedOut),genesOut];
output = sortrows(output,6,'desc');

clear allPathways combined combinedOut data directories directory directory1
clear directory2 fam famOut fdr fdrOut filename frac fracOut genes genesOut
clear i j k loc loc1 loc2 ls path pathOut pathway pathwayFamily pathwayFamilyOut
clear pathwayOut pathwayOut1 pathwayOut2 pOut pval uniquePathways val1 val2
        
%% Bar figures

data = input;
colorMap = cmap;
cmapRange = [0 7];

combined = cell2mat(data(:,6));
combined = log(combined);

pval = cell2mat(data(:,4));
pval = -log10(pval);
pval = [cmapRange(2); pval; cmapRange(1)];
pval = ceil(rescale(pval,1,length(cmap(:,1))));
pval = pval([2:end-1]);



figure
hold
for i = [1:length(combined(:,1))]
    ymid = length(combined) + 1 - i;
    combinedScore = combined(i,1);
    xval = [0 combinedScore combinedScore 0];
    yval = [ymid-0.4 ymid-0.4 ymid+0.4 ymid+0.4];
    
    pInd = pval(i,1);
    pColor = colorMap(pInd,:);
    
    f = fill(xval,yval,pColor);
end
set(gca,'Visible','off');    
    
clear cmapRange colorMap combined combinedScore data f i pColor pInd
clear pval xval xmid yval


        
        
        
        
        
        















