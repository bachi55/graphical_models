addpath ('library/QUIC/');
addpath ('GUI/');

dataTable = loadDataset ('trees');

g = Gaussian2 (dataTable);
g.disp ();
g.plot (dataTable);

g.marg ({'Girth', 'Volume'});
g.disp ();

g = Gaussian (dataTable);
g.cond ({'Volume'}, 5);
g.disp ();

sg = SparseGaussian (dataTable, 0.95);
sg.disp();
sg.plot();

dataTable = loadDataset ('iris');
dataTable = dataTable (:, 1:4);
mg = MixtureOfGaussians (dataTable, ...
    struct ('components', 3, 'maxIt', 500, 'lhTol', 1e-6, 'verbose', 1));
mg2 = mg; 
mg2.marg ({'Sepal_Length', 'Petal_Length'});

mg.plot (dataTable);
mg.aggregation();
mg.plotCondRVPair ({'Sepal_Length', 'Petal_Length'}, {'Sepal_Width', 3.5});

mg2 = mg;
mg2.cond ({'Sepal_Width'}, 3.5);
mg2.plot (dataTable);

mg.disp();
mg.marg ({'Sepal_Width', 'Petal_Length'});

dataTable = loadDataset ('trees2');
data = dataTable (:, [1, 2, 3]);
mg = MixtureOfGaussians (data, ...
    struct ('components', 3, 'maxIt', 500, 'lhTol', 1e-6, 'verbose', 1));
mg.plot (data);
figure;
mg.plot (data, table2array (dataTable(:, 5)));
mg.aggregation();

dataTable = loadDataset ('trees2');
dataTable = dataTable (:, [1:3, 5]);
cg = ConditionallyGaussian (dataTable, ...
    struct ('categorial', {{'pl'}}, ...
            'continuous', {{'crown_diameter', 'tree_height', 'trunk_diameter'}}), ...
    struct ('verbose', 1));

dataTable = loadDataset ('australian-crabs');
dataTable = dataTable (:, [1, 2, 4, 5, 6, 7, 8]);
cg = ConditionallyGaussian (dataTable, ...
    struct ('categorial', {{'species', 'sex'}}, ...
            'continuous', {{'FL', 'RW', 'CL', 'CW', 'BD'}}), ...
    struct ('verbose', 1));

cg.marg ({'species', 'FL'});
