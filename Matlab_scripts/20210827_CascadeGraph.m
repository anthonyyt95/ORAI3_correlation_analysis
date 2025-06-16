%% Creates cascading bar graph

data = input;
qGenes = {'ORAI3';'ORAI1';'ORAI2'};
lineCol = [0.2 0.2 0.2;
    0.5 0.5 0.5;
    0.8 0.8 0.8];

% Organize data
genes = data([2:end],1);
values = cell2mat(data([2:end],[2:end]));
cells = data(1,[2:end]);

for i = [1:length(genes(:,1))]
    val = genes{i,1};
    val = char(val);
    genes{i,1} = val;
end

% Acquires ORAI genes
valOut = [];
for i = [1:length(qGenes(:,1))]
    gene = qGenes{i,1};
    loc = find(strcmp(lower(genes),lower(gene)));
    valOut = [valOut; values(loc,:)];
end

% Acquires ORAI3/sum(ORAI) metric
oraiMetric = [];
for i = [1:length(valOut(1,:))]
    val = valOut(:,i);
    metric = val./sum(val);
    oraiMetric = [oraiMetric, metric];
end

% Orders cells according to ORAI3
valOut = [cells; num2cell(oraiMetric)];
valOut = valOut';
valOutOrdered = sortrows(valOut,2,'desc');

% Creates cascading bargraph
figure
hold
xval = 1;
for i = [1:length(valOutOrdered(:,1))]
    val = cell2mat(valOutOrdered(i,[2:end]));
    
    lastVal = 0;
    for j = [1:length(val)]
        yval = lastVal + val(j);
        p = plot([xval,xval],[lastVal,yval]);
        p.LineWidth = 4;
        p.Color = lineCol(j,:);
        
        lastVal = yval;
    end
    xval = xval + 1;
end
set(gca,'TickDir','out');
set(gca,'XTick',[]);
set(gca,'YTickLabels',[]);
set(gcf,'Position',[500 500 1000 300]);
ylim([0,1]);
xlim([0,length(valOut(:,1))]);

exportgraphics(gcf,'figure.pdf');

clear cells data gene genes i loc metric oraiMetric qGenes val valOut 
clear valOutOrdered values j lastVal lineCol p xval yval

    
    



