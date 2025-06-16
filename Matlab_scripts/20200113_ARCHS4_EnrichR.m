%% Acquires all pathways in a file containing gene X

filename = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\05.20.18_ORAI3 Analyses\Data\2020.01.11_EnrichRPathways\2020.01.13_EnrichRpathways.txt';
gene = 'ORAI3';

data = importdata(filename);

genePresent = [];
output = {};
step = 1;
for i =[1:length(data(:,1))]
    i
    val = data{i,1};
    ind = strfind(val,'"');
    val(ind) = [];
    val = strsplit(val,',');
    val(end) = [];
    
    pway = val{1};
    descriptor = val{2};
    val = val([3:end]);
    
    loc = find(strcmp(lower(val), lower(gene)));
    if not(isempty(loc))
        output{step,1} = descriptor;
        output{step,2} = pway;
        output{step,3} = strjoin(val,',');
        step = step + 1;
    end
end
  
clear data descriptor filename gene genePresent i ind loc pway step val

%% Searches for pathways containing ORAI3

data = output;

distMatrix = NaN(length(data(:,1)),length(data(:,1)));
for i =[1:length(data(:,1))]
    i
    val1 = data(i,:);
    pway = [val1{1,1},': ',val1{1,2}];
    val1 = strsplit(val1{1,3},',');
    
    % Contributes to distance matrix
    for j = [i:length(data(:,1))]
        val2 = data(j,:);
        val2 = strsplit(val2{1,3},',');
        
        l1 = length(val1);
        l2 = length(val2);
        if l1 < l2
            len = l1;
        elseif l2 < l1
            len = l2;
        else
            len = l1;
        end
        
        index = length(intersect(val1,val2))/len;
        distMatrix(j,i) = index;
    end
        
end

clear data descriptor dirName filename gene i j list loc loc1 loc2
clear pway step val k tempVal

%% Performs clustering anaysis based on distance matrix

distMatrixFile = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\05.20.18_ORAI3 Analyses\Matlab\2020.01.13_EnrichR_ORAI3pathways.csv';
pathwayFile = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\05.20.18_ORAI3 Analyses\Data\2020.01.11_EnrichRPathways\2020.01.13_EnrichRpathways.txt';

pathways = importdata(pathwayFile);
distMat = importdata(distMatrixFile);






%% EnrichR pathway file cleanup

filename = 'C:\Users\antho\Downloads\Achilles_fitness_decrease.txt';
data = importdata(filename,',');

output = {};
txtFile = '';
for i = [1:length(data(:,1))]
    val = data{i,1};
    val = regexprep(val,'\t',',');
    loc1 = strfind(val,',,');
    if not(isempty(loc1))
        pway = val([1:loc1-1]);
        val([1:loc1+1]) = [];
    end
    
    line = [pway, val];
    output{i,1} = pway;
    output{i,2} = val;
    
    loc1 = strfind(filename,'\');
    loc1 = loc1(length(loc1));
    loc2 = strfind(filename,'.');
    descriptor = filename([loc1+1:loc2-1]);
    
    if i == 1
        line = [pway,',',descriptor,',',val];
        txtFile = [txtFile, line];
    else
        line = ['?',pway,',',descriptor,',',val];
        txtFile = [txtFile, line];
    end
end
txtFile = regexprep(txtFile,'?','\n');
txtFile = regexprep(txtFile,',','\t');
dlmwrite([descriptor,'.gmt'],txtFile,'delimiter','');



%% Create comma-delimited text file from above output

data = output_annot;

txtFile = '';
for i = [1:length(data(:,1))]
    pway = data{i,1};
    descriptor = data{i,2};
    list = data{i,3};
    i
    if i == 1
        line = [pway,',',descriptor,',',list];
    else
        line = ['>',pway,',',descriptor,',',list];
    end
    txtFile = [txtFile,line];
end
txtFile = regexprep(txtFile,',','\t');
txtFile = regexprep(txtFile,'>','\n');
dlmwrite('2020.01.13_EnrichRPathways.txt',txtFile,'delimiter','');


