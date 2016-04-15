% To normalize RACIPE predicted gene expression data

function Data_normalization(name, numG, maxstates, varargin)

if nargin == 3
    meanV = 0;
    stdV  = 0;
elseif nargin == 5
    meanV = varargin{1};
    stdV  = varargin{2};
end

fsname = [name, '_zscore.dat'];

fidout = fopen(fsname, 'w');
fprintf(fidout, '%s', 'name');
for i = 1:numG 
    fprintf(fidout, '\t%s', ['Gene' num2str(i)]);
end
fprintf(fidout, '\n');

len = zeros(1, numG);

datacomb = [];

for i = 1:maxstates
    tmpname = [name, '_', num2str(i), '.dat'];
    tmpdata = importdata(tmpname);
    tmpdata = rearragendata(tmpdata(:,4:end), i, numG);

    len(i) = size(tmpdata, 1);
    
    datacomb = [datacomb; tmpdata];
end

tmpdata2 = center_data(datacomb(:, 3:(2+numG)), numG, meanV, stdV);

tmp = 1;

for i = 1:maxstates
    for j = tmp:sum(len(1:i))
        fprintf(fidout, '%s', [num2str(i), '_', num2str(datacomb(j,1)), '_', num2str(datacomb(j,2)) ]);
        for h = 1:numG
            fprintf(fidout, '\t%f', tmpdata2(j, h));
        end
        fprintf(fidout, '\n');
    end
    tmp = tmp + len(i);
end

fclose(fidout);
        
end

function out = rearragendata(data, cnt, numG) 
    r = size(data, 1);
    out = zeros(cnt*r ,2+numG);

    for i = 1:r
        for j = 1:cnt
            out((i-1)*cnt + j, 1) = i;
            out((i-1)*cnt + j, 2) = j;
            out((i-1)*cnt + j, 3:(2+numG)) = data(i,(1 + (j-1)*numG):(numG + (j-1)*numG));
        end
    end
end

function out = center_data(data, numG, meanV, stdV)
    out = zeros(size(data));
    
    if meanV == 0
        for i = 1:numG
            out(:, i) = zscore(data(:, i));
        end
    else
        for i = 1:numG
            tmpout   = data(:,i) - meanV(i);
            out(:,i) = tmpout./stdV(i); 
        end
    end
end
