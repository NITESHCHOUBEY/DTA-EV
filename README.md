# Capacitated Equilibrium for Dynamic Traffic Assignment
This code computes an approximate `dynamic capacitated equilibrium` for multicommodity networks with resource
constrained paths between given source-destination pairs. A fixed-point algorithm is implemented that
starts from a feasible path-based flow and iteratively converges to an approximate dynamic equilibrium by discretizing
the continuous time-scale.

# Installation (Linux)
* Install Python 3.8 (if not already available) using `sudo apt install python3.8`
* Install the [NetworkX](https://networkx.org/) package using `pip3 install networkx`

# Preparing Input Files
* Network

The code reads the network using a row-wise description of its edges in a space-separated format as follows:

`nodeFrom nodeTo edgeCapacity edgeTravelTime edgeResourceConsumption edgePrice`

* Commodities

Different commodities are described as rows of a simple space-separated text file in the following format:

`nodeFrom nodeTo resourceBudget priceBudget time0 time1 time2 ... flowValue0_1 flowValue1_2 ...`

The input flow is assumed to be piecewise constant and is described using the time-breakpoints `time0 time1 ...`
and the values of the inflow between consecutive breakpoints (for example, `flowValue0_1` indicates the inflow
in the time interval `[time0, time1)`.


# Running the Code
The code can be run using the script `main.py` with parameters passed as space separated arguments as follows:

`python3 main.py <networkFile> <commodityFile> <instanceName> <timeHorizon> <iterationLimit> <timeLimit> <precision> <alpha0> <timeStep> <priceToTime> <numThreads>`

where,

`networkFile`: a plain text file containing the attributes of edges per line,

`commodityFile`: a plain text file describing each commodity per line,

`instanceName`: a name to indicate the instance,

`timeHorizon`: an approximate time up to which the last agent to expected to reach the sink,

`iterationLimit`: the maximum number of iterations to be run by the fixed-point algorithm,

`timeLimit`: the maximum (wall clock) time limit for the code after which it will be terminated,

`precision`: the change in the L-1 norm of path flows desired in the path inflows returned in consecutive iterations,

`alpha0`: starting value of the step-size parameter of the fixed-point algorithm,

`timeStep`: discretization time-step of the fixed-point algorithm,

`priceToTime`: a scalar to converts price to equivalent travel-time units,

`numThreads`: the number of threads to be used in the code (default 1).


# A Toy Example

Consider the following network.

![Toy Example](examples/toyExamples/evExample6.png)

The networkFile for reading this network is [evExample6EdgesWR.txt](examples/toyExamples/evExample6EdgesWR.txt). An example of a commodityFile for this network is <a href="examples/toyExamples/evExample6Comm.txt" target="_blank"> evExample6Comm.txt</a>. The lines starting with '#' in the files are ignored by the code. Also, the code can read the values of edge attributes as rationals provided in the form `p/q` where `p` and `q` are integers.

To find an approximate equilibrium for the above mentioned inputs, we can run the following command:

`python3 main.py ./examples/toyExamples/evExample6EdgesWR.txt ./examples/toyExamples/evExample6Comm.txt example6 50 2 3600 0.001 0.75 0.25 0 1`

'python3 main.py ./examples/nguyen/nguyenEdges.txt ./examples/nguyen/nguyenComm.txt example6 50 40000 7200 0.01 0.005 0.25 0 1'


