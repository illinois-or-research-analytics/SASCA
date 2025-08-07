# SASCA: Scalable Agent-based Simulator for Citation Analysis

## One time setup
run the `setup.sh` script to locally install [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and [PCG](https://www.pcg-random.org/).

## How to run
The general command is given below.
```console
abm --edgelist ${INPUT_EDGELIST} --nodelist ${INPUT_NODELIST} --out-degree-bag ${OUTDEGREE_BAG} --recency-probabilities ${RECENCY_PROBABILITIES} --alpha ${ALPHA_PROPORTION} --growth-rate ${GROWTH_RATE} --fully-random-citations ${FULLY_RANDOM_CITATIONS_PROPORTION} --same-year-proportion ${SAME_YEAR_PROPORTION} --num-cycles ${NUM_CYCLES} --output-file ${OUTPUT_EDGELIST} --auxiliary-information-file ${OUTPUT_AUX_FILE} --log-file ${OUTPUT_LOG} --num-processors ${NUM_THREADS} --log-level ${LOG_LEVEL}
```
- Input edgelist is a two column csv with source and target integer node ids on each row.
- Input nodelist is a two column csv with integer node ids in the first column and the corresponding integer publication year in the second column.
- Outdegree bag is a two column csv with integer columns. The first column is ignored. The second column defines the out-degree distribution. The simulator will effectively pick one line at random and use the second column as a new agent's outdegree.
- Recency probabilities is a two column csv with integer value `n` and integer value `k` in each row. Engineered to be empirically driven by recording the number of citations (`k`) that are made to publications from `n` years ago. for a real-world citation network. More generally, it encodes the users belief about the probability of a publication being cited given the age of the publication.
- Alpha determines the proportion of citations that are made to the 1-hop node of the generator node relative to the number of 2-hop citations. It can be left to be -1 for a random model or a constant value such as 0.5.
- Growth rate is the exponent for the exponential growth formula. It controls how fast the network grows.
- Fully random citations is a proportion, typically set to 0.05, which determines the proportion of citations that are made uniformly random to the entire network.
- Same year proportion, typically set to 0.12, determines the number of new agents that make a single citation to other new agents in the same year.
- Num cycles is the time in the exponential growth formula. Combined with the growth rate and input nodelist size, num cycles determinsitcally determines the number of nodes in the final network.
- Output file is the two column CSV edgelist output representing the final network.
- Auxiliary information file is a multi-column CSV output that contains metadata about each node such as its assigned out-degree, alpha, agent weights, or the final in-degree and out-degree.
- Log file has the logs, the most useful of which is the size of network at the beginning of each cycle and time elapsed.
- Num processors directly controls the number of new agents that are processed in parellel for each year.
- Log level is an integer with values 0, 1, and 2 corresponding to silent, info, and verbose modes.
