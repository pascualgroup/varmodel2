# varmodel2

Ed Baskerville

Last update: 26 October 2017

This code is an implementation of a model of malaria var gene evolution within an epidemic simulation.
Malaria strains are represented as unordered sets of var genes, which are in turn composed of abstract loci.
A number of alleles can appear at each locus, and the allelic composition of a gene across loci governs immune dynamics in the host.
Individual hosts are infected by strains, and infections can be transmitted between hosts.
Each infection expresses a single var gene at a time, and the sequence of expressions is explicitly represented in the simulation.
The simulation also includes immigration of new strains into the population, recombination during transmission and during an infection, and mutation.

The simulation is modeled as a sequence of discrete events (state changes) that happen in continuous time.

## Building and running: cheat sheet

1. Create a parameters file in Python (`params.py`) format, or generate one (or many) in JSON format (`params.json`).

2. Build the code using the specified parameters into run directory `rundir`:

```sh
path/to/varmodel2/build.py -p params.py -d rundir
```

3. Run the model:

```sh
cd rundir
bin/varmodel2
```

If you want to specify a different random seed without recompiling the model, just do:

```sh
bin/varmodel 1234
```

## Building and running the model: details

### Each run lives in a separate working directory

The code is designed to be used with a fresh working directory for each run that contains output as well as the compiled model, source code, and generated source code.
The file `parameters.hpp` is generated before compilation from a file containing parameter values in either Python or JSON format.
Storing a copy of the code with each run ensures that you know exactly what version of the code and what parameter values were used to generate a set of outputs.

One exception: if you do replicates with the same parameter values, it is possible to do multiple runs with different random seeds using the same compiled model.

### Setting up parameters

To do an individual run, first set up a parameter file as a simple Python file (module).
An example is provided in `parameters-example.py`.

This does not require you to know much Python; parameter values are simply assigned to variable names.
The parameter names should match variables defined in `src/parameters.hpp.template`, and the values should match the appropriate type: e.g., floating-point numbers for `double`; integers for `int64_t` and `uint64_t`; lists for `array`; etc.

When you build the model, parameter values are inserted directly into `parameters.hpp.template` to create `parameters.hpp`.

E.g., since `parameters.hpp.template` contains the line

```cpp
const double {P_IMMIGRATION_INCLUDES_NEW_GENES};
```

a parameters file, e.g., `myparameters.py`, must contain a line like

```python
P_IMMIGRATION_INCLUDES_NEW_GENES = 0.5
```

which will result in the following line in the generated file `parameters.hpp`:

```python
const double P_IMMIGRATION_INCLUDES_NEW_GENES = 0.5;
```

Strings are passed verbatim—without quotation marks—which enables the use of both enumerated (`enum`) types and strings.

E.g., if you want to pass an `enum` value, use a regular Python string:

```python
SELECTION_MODE = 'SPECIFIC_IMMUNITY'
```

which yields

```cpp
const SelectionMode SELECTION_MODE = SPECIFIC_IMMUNITY;
```

and if you want to pass a string value, add double-quote (`"`) marks inside the string:

```python
SAMPLE_DB_FILENAME = '"output_samples.sqlite"'
```

which yields

```cpp
const std::string SAMPLE_DB_FILENAME = "output_samples.sqlite";
```

to specify a string, with inner quotation marks.

You can also write parameter files in JSON format, e.g.:

```javascript
{
    // ...
    "SAMPLE_DB_FILENAME" : "\"output_samples.sqlite\"",
    "SELECTION_MODE" : "SPECIFIC_IMMUNITY",
    "P_IMMIGRATION_INCLUDES_NEW_GENES" : 0.5,
    // ...
}
```

This is most useful as an output format for parameter sweep scripts that generate many parameter files.

### Building the model

To build the code with the specified parameters, do:

```sh
path/to/varmodel2/build.py -p params.py -d rundir
```

This will do the following things:

* Generate the file `rundir/generated/parameters.hpp` by inserting values from `params.py`
* Generate `rundir/generated/*Manager.*` files, which contain code to manage objects of different types and save/load them to checkpoint databases
* If all goes well, compile the model into `rundir/bin/varmodel2`.

### Running the model

You can now run the model via:

```sh
cd rundir
./bin/varmodel2
```

which will generate an output database in SQLite format based on the parameter `SAMPLE_DB_FILENAME`, and will periodically write checkpoints to `CHECKPOINT_SAVE_FILENAME` if checkpointing is on.

## Model specification

### Overview of structure

In this model, hosts carry infections of different strains of the malaria parasite.

Strains are composed of a collection of `N_GENES_PER_STRAIN` var genes.
Strain identity is defined by this collection independent of order.
Although unlikely, the same gene may occur multiple times in a strain.

Genes are composed of `N_LOCI` loci.
At the outset, each locus `i` has one of `N_ALLELES_INITIAL[i]` possible values, indexed from `0` to `N_ALLELES_INITIAL[i] - 1`.
Mutation events create new alleles, so the number of distinct alleles at each locus can increase over time.

At any time, hosts may be infected multiple times by the same or different strains.
Each infection defines a sequence of expression of var genes by a strain.
The order of expression is randomized distinctly for each infection.

Hosts have an immune history, the details of which depend on the immune selection model being used.

### A discrete-event, continuous-time model

The simulation consists of a sequence of discrete events occurring in continuous time.
Conceptually, the state of a model consists of the collection of hosts, some of whom contain infections, along with a single *event queue* consisting of events that will occur in the future at specified times.
The simulation progresses by looking at the event on the queue with the lowest time, advancing the clock to that time, and then executing the state change corresponding to that event.
State changes will typically modify both host/infection state as well as the set of future events.

In Pythonic pseudo-code, the entire simulation looks approximately like this:

```python
now = 0.0
queue = initialize_event_queue()
while queue.next_event_time() <= T_END:
    event_time, event = queue.pop_next_event()
    now = event_time
    execute_event(event)
```

In fact, in order to avoid language features that might be confusing to users unfamiliar with C++, this version of the code implements the event queue as *multiple* event queues, one for each type of event.
First, the code looks across all event queues to find the one with the lowest next event time.
Then, it executes just that single event.

Event queues are implemented as indexed priority heaps, as in the [next-reaction method by Gibson and Bruck](https://scholar.google.com/scholar?cluster=16469525020041545503) (see "Event queue implementation", below).

Most events in the simulation are modeled as Poisson processes, so the times associated with those events on the queue are probabilistic realizations: specifically, times are drawn from an exponential distribution.
If a state change causes the underlying rate of the event's Poisson process to change, the time will be re-drawn using the new rate, following the memoryless property of Poisson processes.

In this version of the code, each exponential draw and re-draw is implemented explicitly.
(This is less automatic than the previous version of the code, but it is somewhat easier to see exactly what is going on.)

### Overview of simulation

1. Initialize the population hosts and initial infections.
2. Add all initial events to the event queue.
3. Repeatedly execute events until the simulation is done.

These are the different events that can occur:

* *Biting events* represent the transmission of strains between two hosts.
* *Immigration events* represent the appearance of a new strain from outside the population.
* *Immunity loss events* represent the loss of previously gained immunity.
* *Infection transition events* represent a transition in an infection between the expression of two different var genes.
* *Mutation events* mutate a gene within an infection.
* *Recombination events* recombine genes within an infection.
* *Death events* represent the death of a host, and the birth of a new host in its place.

### Initialization

### Biting events

### Immigration events

### Mutation events

### Recombination events

### Death events

### Immune selection modes

### Expression dynamics

## Event queue implementation

Event queues are implemented using an *indexed priority heap*, an efficient data structure for retrieving the next item in a collection according to some property (key), in this case, the next time of an event, and for adding, removing, and updating items in the collection.

A *heap* is a binary tree, typically implemented using an array using the following indexing scheme:

```
           0
     1           2
  3     4     5     6
07 08 09 10 11 12 13 14
```

Each event in the tree has the property that its time is less than its two children.
E.g., the event at index 6 has time less than the events at indices 13 and 14;
the event at index 2 has time less than the events at indices 5 and 6;
and the event at index 0 has time less than the events at indices 1 and 2.

This means that, if the heap is in a consistent state, the event with the lowest time is always stored at index 0.

The heap is *indexed* using a hash table mapping event IDs to the location of the event in the heap.
This makes it possible to find an event whose time needs updating in the tree in constant time.

When an event time is changed, the heap can be made consistent again by moving the event to the appropriate height in the tree by swapping parents and children.
This process is limited by the height of the tree, meaning that update operations are logarithmic in the number of events; this is why the heap is so efficient.

## History and overview of changes

This code is a new, simpler implementation of the malaria var gene evolution model by Qixin He, in which var genes are composed of loci with varying alleles.
That model code was based on a model by Yael Artzy-Randrup, implemented by Ed Baskerville, in which var genes are simply distinguished from each other by identity.

The main changes from the previous implementation are as follows:

* Some details of model behavior have changed (see section below).
* Parameters are compile-time constants instead of being dynamically loaded from JSON at runtime, which added unnecessary complexity to the C++ code.
A Python script is used to generate the `parameters.hpp` constants file from either JSON or from a Python module.
* The abstraction layer for database output was removed.
Output tables (for sampling) are written using raw calls to the SQLite API.
* You can save and load checkpoints&mdash;complete on-disk representations of the simulation state.
The checkpointing code is automatically generated by a Python script that parses C++ data structures.
Checkpoints are SQLite databases following a simple object relational model (ORM) for the simulation state, and are suitable as a basis for more detailed analysis than periodic sampling.
* C++ language features are used only where they clearly make the code better without impairing readability for people with limited familiarity with C++.

## Code generation for parameters and database output

## Debugging and testing

Debugging output is produced via the `PRINT_DEBUG` macro, whose first argument is an integer debugging level.

The `PRINT_DEBUG_LEVEL` parameter controls what debugging output is actually produced.

Note that turning on all debugging output will substantially slow down the simulation, not to mention.
(Because `PRINT_DEBUG_LEVEL` is a compile-time constant, the compiler will eliminate all unused debugging statements from the code.)

Before each simulation, a small set of unit tests are run on parts of the code to ensure they are working correctly.

Additionally, every `VERIFICATION_PERIOD` in simulation time units, data structures are checked for integrity.

## Behavioral changes from previous model

All changes listed here are relative to [VarModel repo commit f66d253](https://github.com/pascualgroup/VarModel/commit/f66d253176960a539db9628c1a7aeaa7fa4ab6f1).

### Elimination of inactive state between gene expressions

There is no longer an "inactive" state between gene expressions.
Instead, state transitions go directly from the liver (waiting) stage to the first expression, then the second expression, and so on.

### Dead code removal

TODO: list unused functionality that was removed

### Bug fix: off-by-one error in ectopic recombination

In the previous version of the code, there was an off-by-one error in choosing the crossover point.

The crossover point was never chosen to be just before the last locus; rather, the farthest along it could be was before the second-to-last locus.
This has been fixed in the new version of the code.

### Bug fix in recombination during transmission

In choosing the two strains to be recombined, the recombination probability was set assuming that a distinct pair of strains would be selected, but in fact the same pair could be selected for recombination.

Now, the code simply chooses two random strains.
If they are the same, no recombination occurs; if they are different, recombination occurs.
