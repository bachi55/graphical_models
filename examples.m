%% prepare workspace
addpath ('QUIC/');
addpath ('GUI/');

%% Gaussian
dataTable = loadDataset ('trees');

g = Gaussian2 (dataTable);
g.plot (dataTable);

g.marg ({'Girth', 'Volume'});
g.disp;
g.cond ({'Volume'}, 5);
g.disp;

%% Sparse Gaussian
sg = SparseGaussian (dataTable, 0.95);
sg.disp;
sg.plot;

%% Mixture Of Gaussians
dataTable = loadDataset ('iris');
dataTable = dataTable (:, 1:4);
mg = MixtureOfGaussians (dataTable, ...
    struct ('components', 3, 'maxIt', 500, 'lhTol', 1e-6, 'verbose', 1));

% plot without data
mg.plot;
% plot with data
mg.plot (dataTable)

% use the inspector
% NOTE: Always change the variables in order: RV1 -> RV2 -> RVCond
mg.inspect;

%% Conditional Gaussians
dataTable = loadDataset ('trees2');
dataTable = dataTable (:, [1:3, 5]);
cg = ConditionallyGaussian (dataTable, ...
    struct ('categorical', {{'pl'}}, ...
            'continuous',  {{'crown_diameter', 'tree_height', 'trunk_diameter'}}), ...
    struct ('verbose', 1));

dataTable = loadDataset ('australian-crabs');
dataTable = dataTable (:, [1, 2, 4, 5, 6, 7, 8]);
cg = ConditionallyGaussian (dataTable, ...
    struct ('categorical', {{'species', 'sex'}}, ...
            'continuous',  {{'FL', 'RW', 'CL', 'CW', 'BD'}}), ...
    struct ('verbose', 1));