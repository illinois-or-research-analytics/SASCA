# ABM

## How to run
The general command is given below.
```
abm --edgelist ${INPUT_EDGELIST} --nodelist ${INPUT_NODELIST} --out-degree-bag ${OUTDEGREE_BAG} --recency-probabilities ${RECENCY_PROBABILITIES} --alpha ${ALPHA} --growth-rate ${GROWTH_RATE} --num-cycles ${NUM_CYCLES} --output-file ${OUTPUT_FILE} --log-file ${OUTPUT_LOG} --num-processors ${NUM_THREADS} --log-level ${LOG_LEVEL}
```

Input edgelist is a two column csv with source and target integer node ids on each row.

Input nodelist is a two column csv with integer node ids in the first column and the corresponding integer publication year in the second column.

Outdegree bag is a two column csv with integer columns. The first column is not read. The second column is a sample from an outdegree distribution. The program will pick one line at random and use the second column as a new node's outdegree.

Recency probabilities is a two column csv with integer value n and floating point value k in each row. It encodes the users belief about the probability (k) of a publication being cited given the age (n) of the publication.

Alpha can be left to be -1 for a random model or a constant value such as 0.5.

Growth rate is the exponent for the exponential growth formula. It controls how fast the graph will grow.

Num cycles is the time in the exponential growth formula.

## Interface
### `void ReadOutDegreeBag();`
This function reads one of the inputs to the program. The file read specifies a possible out-degree per line of file. When an agent is initilaized, a random line is chosen from this file and used as the out-degree.

### `void ReadRecencyProbabilities();`
This function reads one of the inputs to the program. The file read specifies the probability of citing a node from the past where the first column specifies the number of years in the past while the second column specifies the probability of citing a node from that number of years ago.


### `std::map<int, int> BuildContinuousNodeMapping(Graph* graph);`
This function is a standard mapper function that takes in a graph pointer and outputs a mapping from node ids in the graph to continuous integer ids from 0 to graph size.

### `std::map<int, int> ReverseMapping(std::map<int, int> mapping);`
This function simply takes in a map of int to int and reverses the keys and values. This is used in conjunction with `BuildContinuousNodeMapping()` and should never run into key collisions as long as the continious mapping was built correctly.

### `std::vector<int> GetComplement(Graph* graph, const std::vector<int>& base_vec, const std::map<int, int>& reverse_continuous_node_mapping)'
This function takes in a vector of node ids and gets the node ids in the graph that are not in the input base vector.

### `std::vector<int> GetNeighborhood(Graph* graph, const std::vector<int>& generator_nodes, const std::map<int, int>& reverse_continuous_node_mapping)`
This function takes the input generator node vector and returns the k-hop neighborhood of all the generator nodes.

### `std::vector<int> GetGeneratorNodes(Graph* graph, const std::map<int, int>& reverse_continuous_node_mapping)`
This function picks generator nodes.

### `int GetFinalGraphSize(Graph* graph);`
The ABM model has a set growth rate and an initial size. Therefore, it is possible to calculate the final graph size once all of the input parameters are read in. This function is responsible for calculating the exact final graph size which is then used to allocate enough memory space for score arrays and other data structures needed.

### `void FillInDegreeArr(Graph* graph, const std::map<int, int>& continuous_node_mapping, int* in_degree_arr);`
This function takes an empty array (`in_degere_arr`) that has already been allocated and populates it with the in-degrees of each node according to the node mapping indices.

### `InitializeFitness(Graph* graph);`
This function calls three other functions that are templated. The three functions are responsible for initializing fitness related attributes for each node in the graph.

### `void FillFitnessArr(Graph* graph, const std::map<int, int>& continuous_node_mapping, int current_year, int* fitness_arr);`
Fitness is a dynamically calculated value where each node has three phases. During th e first phase, a node has the minimum possible fitness value which is 1 in our current simulator. After the first phase is over, the fitness value assigned to a node immediately climbs to the maximum value assigned. This value is assigned per node from another function. T

### `void FillRecencyArr(Graph* graph, const std::map<int, int>& continuous_node_mapping, int current_year, double* recency_arr);`
This function takes an already allocated `recency_arr` and populates it with probabilities for each node where it's the normalized probabilites for the current year - existing node year obtained from one of the input files.

### `void PopulateWeightArrs(double* pa_weight_arr, double* rec_weight_arr, double* fit_weight_arr, int len);`
This function creates randomized agent weights that sum to 1 to each of the nodes in the graph. The size of the graph given by `len` which is also the length of the weight arrays passed in.

### `PopulateAlphaArr(double* alpha_arr, int len);`
This function populates `alpha_arr` with either random doubles from 0 to 1 or a set value passed in from the program arguments.

### `int GetMaxYear(Graph* graph);`
This function is a utility function that is called once in order to identify the last year in the seed graph.

### `int GetMaxNode(Graph* graph);`
This function is a utility function that is called once in order to identify the highest node id in the seed graph.

### `void PopulateOutDegreeArr(int* out_degree_arr, int len);`
This function populates randomly draws an out degree from the input out degree file for the seed graph.

### `void CalculateScores(int* src_arr, double* dst_arr, int len);`
This function takes each element from `src_arr` and raises it to the power of `gamma` and adds 1 before storing in `dst_arr`.

### `int MakeCitations(int* citations, double* score_arr, double* random_weight_arr, int len, int num_citations);`
This function takes the `score_arr` as weights for weighted random sampling without replacement and stores `num_citations` randomly drawn indices to `citations`.

### `void AssignFitnessLagDuration(Graph* graph, const T& container)`
This function assigns a random positive integer to be the length of fitness lag duration for each node in the graph.

### `void AssignFitnessPeakDuration(Graph* graph, const T& container)`
This function assigns a random positive integer to be the length of fitness peak duration for each node in the graph.

### `void AssignPeakFitnessValues(Graph* graph, const T& container)`
This function assigns a random positive integer to be the value of peak fitness value for each node in the graph.

## Formulas and Algorithms
### Fitness values
Fitness values are now dynamically calculated based on the current year. Every node in the graph has an initial "Lag" duration during which the fitness value for a node stays at the minimum possible value. Afterwards, the fitness value immediately climbs to the peak fitness value. This value is sustained for some number of years until it starts decaying.

### Weighted random sampling without replacement
We are continuing to explore the random sampling algorithm part.
