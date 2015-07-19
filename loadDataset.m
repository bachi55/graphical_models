function dataTable = loadDataset (modelName, modelDir) 
    if nargin < 2
        modelDir = '/home/bach/Documents/studies/ss15/graphische-modelle/projekt2/data/';
%         modelDir = '/home/va54tuz/Documents/ss15/graphische-modelle/projekt2/data/';
%         modelDir = 'M:\Documents\ss15\graphische-modelle\projekt2\data\';
    end % if

    dataTable = readtable (sprintf ('%s/%s.csv', modelDir, modelName));
end % function