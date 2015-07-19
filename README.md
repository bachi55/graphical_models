Graphical-Models
================

Projected created during the course Graphical-Models in the summer-term 2015. 

Description of the project:
---------------------------
Aim of the project was to implement so called *queries* for different 
probability models. 

Implementation strategy:
------------------------
The different probability models are implemented as Matlab classes. This 
provides an efficient way to re-use methods like *marginalize*, 
*conditioning*, *plotting*, etc. for different models.

E.g.: Class-tree for a Gaussian:

        ProbabilityDensityFunction
        
                |
                V
                
        MixturesOfGaussians --> SparseMixturesOfGaussians (A)
        
                |
                V
                
        Gaussian
        
                |
                V
                
        SparseGaussian (B)
               
(A) not implemented yet, (B) implemented, but not integrated in the 
class hierarchy

Probability Models:
-------------------
(already support *marg*, *cond*, *aggr*, and *plot*)
* Gaussian (FILE: Gaussian2.m) 
* MixturesOfGaussians (FILE: MixturesOfGaussians.m)
* SparseGaussians (FILE: SparseGaussian.m)

(only support parameter estimation)
* ConditionalGaussians

Status of project:
------------------
Everything based on *MixturesOfGaussians* is ready. The model can be 
estimated, and queries can be performed. The *inspect* method allows a 
simple evaluation of the model and shows, that the queries are working. 

The *ConditionalGaussians* are not ready now. Till now only the model 
estimation is implemented. No queries can be performed and the 
'inspector' does not work for them till now. 

Example applications:
---------------------
First change into the project folder. And prepare your workspace:
```
%% prepare workspace
addpath ('library/QUIC/');
addpath ('GUI/');

```

### Gaussian
Than lets estimate and plot a simple Gaussian:
```matlab
dataTable = loadDataset ('trees');
g = Gaussian2 (dataTable);
g.plot (dataTable);
```

Ok, now lets marginalize and conditioning a bit:
```
g.marg ({'Girth', 'Volume'});
g.cond ({'Volume'}, 5);
g.disp;
```

Please not that the implementation of *marg* marginalizes the random 
variables **not** given. 

### Mixture of Gaussians
But a single Gaussian can be boring:
```
dataTable = loadDataset ('iris');
dataTable = dataTable (:, 1:4);
mg = MixtureOfGaussians (dataTable, ...
    struct ('components', 3, 'maxIt', 500, 'lhTol', 1e-6, 'verbose', 1));
mg.plot (dataTable);
```

Also here we could marginalize and conditioning, but we will use now the 
*inspect* function to make a bit more fancy:
```
mg.inspect;
```

This simple tool allows us to chose two random variables we want to keep 
after marginalization and one to condition on. The slide-bar at the 
bottom changes the value of the conditioned variable. Also some 
adjustments for the plot can be done. 

### Sparse Gaussian
If we use the L1-regularization to estimate a spare precision matrix 
(the inverse of the sigma-matrix), we obtain the adjacency matrix for 
the Hammersley-Clifford Graph. This graph consists of the random 
variables as nodes and there conditioned independence as edges. Lets 
inspect such a graph: 
```
dataTable = loadDataset ('trees');
sg = SparseGaussian (dataTable, 0.95);
sg.plot;
```

### Sparse Mixture of Gaussians
This class is not implemented yet, but would be straightforward. For 
every component we would need to estimate the L1-regularized precision 
matrix separated.

### Conditional Gaussian (CG)
Unfortunately the implementation for the Conditional Gaussian is not 
ready. But at least we can estimate a model already: 
```
dataTable = loadDataset ('australian-crabs');
% we do not need the index categorical variable
dataTable = dataTable (:, [1, 2, 4, 5, 6, 7, 8]);
cg = ConditionallyGaussian (dataTable, ...
    struct ('categorical', {{'species', 'sex'}}, ...
            'continuous',  {{'FL', 'RW', 'CL', 'CW', 'BD'}}), ...
    struct ('verbose', 1));
```

Interesting data-sets:
----------------------
### trees2
```
dataTable = loadDataset ('trees2');
data = dataTable (:, [1, 2, 3]);
mg = MixtureOfGaussians (data, ...
    struct ('components', 3, 'maxIt', 500, 'lhTol', 1e-6, 'verbose', 1));
```
This data-set contains about 2000 measurements from trees. For every 
tree there are four continuous variables given:
* crown-diameter
* tree-height
* trunk-diameter
* basal-area
Furthermore there is categorical variable given indication the trees' 
species:
* pl

The Mixture of Gaussians can't "catch" out these tree species, since 
they're to mixed into each other. To see this different plots:
```
mg.plot (data);
figure;
mg.plot (data, table2array (dataTable(:, 5)));
```
The first one colors the tree instances according to the most likely 
mixture-component and the second one according to there labels. 

A Conditional Gaussian could include the species as categorical variable 
and using *inspect* could show, how different tree properties "behave" 
for different species.

TODO:
-----
(high priority)
* Add 'marg', 'cond', 'plot' and 'aggr' to the Conditional Gaussians
* Extend the 'inspect'or to handle Conditional Gaussians
* How to plot Conditional Gaussians?
* Add SparseGaussian to the class-hierarchy

(middle priority)
* Exhaustive testing of the classes. Especially for "border-cases". 
* Evaluate the performance of the implementation 

(low priority)
* Create an SQL-like language to talk to the models.

