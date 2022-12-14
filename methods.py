
from qiskit import QuantumCircuit, Aer, execute
import numpy as np
import math
import networkx as nx
import matplotlib.pyplot as plt
from itertools import product as iterprod
from qiskit.tools.visualization import circuit_drawer

import json
import os

loc = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(loc, 'data', 'angles_regular_graphs.json'), 'r') as json_file:
    fixed_angles_data = json.load(json_file)
    
    

def get_size_dist(counts, sizes, optimal_size):
    """ For given measurement outcomes, i.e. combinations of counts and sizes, return counts corresponding to each cut size.
    Args:
        counts (list): List of integers, denoting number of times each cut was measured 
        sizes (list): List of integers. Cut size corresponding to each measured cut
        optimal_size (int) : Max cut size for the graph
    Returns:
        full_counts_list (list) : Number of times each size was obtained
        full_size_list (list) : List of all possible sizes (0,1,...,optimal_size)
    """
    unique_sizes = list(set(sizes))
    unique_counts = [0] * len(unique_sizes)
    
    for i, size in enumerate(unique_sizes):
        # corresp_counts = [counts[ind] for ind,s in enumerate(sizes) if s == size]
        corresp_counts = [c for (c,s) in zip(counts, sizes) if s == size]
        unique_counts[i] = sum(corresp_counts)
    
    # Make sure that the scores are in non-decreasing order
    s_and_c_list = [[a,b] for (a,b) in zip(unique_sizes, unique_counts)] #size and cut list
    s_and_c_list = sorted(s_and_c_list, key = lambda x : x[0]) # sort according to sizes
    unique_sizes = [x[0] for x in s_and_c_list]
    unique_counts = [x[1] for x in s_and_c_list]
    
    full_size_list = list(range(optimal_size + 1))
    full_counts_list = [unique_counts[unique_sizes.index(s)] if s in unique_sizes else 0 for s in full_size_list]
    
    return full_counts_list, full_size_list
    
def compute_energy_expectation(counts, sizes, **kwargs):
    """
    Compute the mean of cut sizes
    This approximates the expectation value of energy

    Parameters
    ----------
    counts : ndarray of ints
        measured counts corresponding to cuts
    sizes : ndarray of ints
        cut sizes (i.e. number of edges crossing the cut)
    **kwargs : optional arguments
        will be ignored

    Returns
    -------
    float : energy expectation
    """
    # Convert counts and sizes to ndarrays, if they are lists
    counts, sizes = np.array(counts), np.array(sizes)

    return np.sum(counts * sizes) / np.sum(counts)

def compute_cvar(counts, sizes, alpha = 0.1, **kwargs):
    """
    Obtains the Conditional Value at Risk or CVaR for samples measured at the end of the variational circuit.
    Reference: Barkoutsos, P. K., Nannicini, G., Robert, A., Tavernelli, I. & Woerner, S. Improving Variational Quantum Optimization using CVaR. Quantum 4, 256 (2020).

    Parameters
    ----------
    counts : ndarray of ints
        measured counts corresponding to cuts
    sizes : ndarray of ints
        cut sizes (i.e. number of edges crossing the cut)
    alpha : float, optional
        Confidence interval value for CVaR. The default is 0.1.

    Returns
    -------
    float
        CVaR value

    """
    # Convert counts and sizes to ndarrays, if they are lists
    counts, sizes = np.array(counts), np.array(sizes)
    
    # Sort the negative of the cut sizes in a non-decreasing order.
    # Sort counts in the same order as sizes, so that i^th element of each correspond to each other
    sort_inds = np.argsort(-sizes)
    sizes = sizes[sort_inds]
    counts = counts[sort_inds]

    # Choose only the top num_avgd = ceil(alpha * num_shots) cuts. These will be averaged over.
    num_avgd = math.ceil(alpha * np.sum(counts))

    # Compute cvar
    cvar_sum = 0
    counts_so_far = 0
    for c, s in zip(counts, sizes):
        if counts_so_far + c >= num_avgd:
            cts_to_consider = num_avgd - counts_so_far
            cvar_sum += cts_to_consider * s
            break
        else:
            counts_so_far += c
            cvar_sum += c * s

    return cvar_sum / num_avgd

def compute_gibbs(counts, sizes, eta = 0.5, **kwargs):
    """
    Compute the Gibbs objective function for given measurements

    Parameters
    ----------
    counts : ndarray of ints
        measured counts corresponding to cuts
    sizes : ndarray of ints
        cut sizes (i.e. number of edges crossing the cut)
    eta : float, optional
        Inverse Temperature
    Returns
    -------
    float
        Gibbs objective function value

    """
    # Convert counts and sizes to ndarrays, if they are lists
    counts, sizes = np.array(counts), np.array(sizes)
    ls = max(sizes)#largest size
    shifted_sizes = sizes - ls
    
    # gibbs = - np.log( np.sum(counts * np.exp(eta * sizes)) / np.sum(counts))
    gibbs = - eta * ls - np.log(np.sum (counts / np.sum(counts) * np.exp(eta * shifted_sizes)))
    return -gibbs


def compute_best_cut_from_measured(counts, sizes, **kwargs):
    """
    From the measured cuts, return the size of the largest cut
    """
    return np.max(sizes)
    
class qaoa_anstaz:
    def __init__(self, nodes, edges, betas, gammas, num_shots = 5000, eta = 0.5, alpha = 0.1, optimal_cut = [], optimal_size = -1):
        """
        Initialize a class for a QAOA ansatz circuit. Create corresponding circuit as well.

        Args:
            nodes (int): number of nodes in the graph
            edges (list of tuples): list of (undirected) edges of the graph. Each element is a two-tuple of node labels, which should be between 0 and (nodes-1) 
            betas (list of floats): List of ?? angles. Number of rounds is equal to the length of this list.
            gammas (list of floats): List of ?? angles. Number of rounds is equal to the length of this list.
            num_shots (int): Number of times the qaoa ansatz circuit is to be measured
            eta (float): Number between 0 and 1. Parameter for Gibbs' objective function
            alpha (float): Number between 0 and 1. Parameter for CVaR function.
            optimal_cut (list): List of length=nodes. Each element is 0 or 1. This list specifies one (out of possibl multiple) answers for the maxcut problem. If not specified, it is computed using a brute force approach while initialization.
            optimal_size (int): The max cut size (i.e. the solution to the maxcut problem). If not specified, it is computed using a brute force approach while initialization.
        """
        self.nodes = nodes
        self.edges = edges
        self.betas = betas
        self.gammas = gammas
        self.num_shots = num_shots
        self.eta = eta
        self.alpha = alpha
        self.nxgraph = nx.from_edgelist(edges)
        
        assert len(self.betas) == len(self.gammas), "Betas and Gammas lists should be of the same length"
        self.rounds = len(self.betas) 
        
        if len(optimal_cut) == 0 or optimal_size == -1:
            # Compute max cut using a brute force approach.
            self.brute_force_compute_maxcut()
        else: 
            # Provided values will be assumed to be correct.
            self.optimal_cut = optimal_cut
            self.optimal_size = optimal_size
        
        self.create_ansatz_circuit()
        
        self._bestcut_ratio = None
        self._gibbs_ratio = None
        self._cvar_ratio = None
        self._approx_ratio = None
        
        self._svs_bestcut_ratio = None
        self._svs_gibbs_ratio = None
        self._svs_cvar_ratio = None
        self._svs_approx_ratio = None
        
        # Do a random uniform sampling of strings and obtain corresponding distribution of cut sizes
        self.random_sampling_dist()
        
        
    @property
    def gibbs_ratio(self):
        if self._gibbs_ratio is None:
            gibbs = compute_gibbs(self.counts, self.sizes, eta = self.eta)
            self._gibbs_ratio = gibbs / self.optimal_size / self.eta
        return self._gibbs_ratio
            
    @property
    def bestcut_ratio(self):
        if self._bestcut_ratio is None:
            best_measured_size = compute_best_cut_from_measured(self.counts, self.sizes)
            self._bestcut_ratio = best_measured_size / self.optimal_size
        return self._bestcut_ratio
    
    @property
    def approx_ratio(self):
        if self._approx_ratio is None:
            energy_expectation = compute_energy_expectation(self.counts, self.sizes)
            self._approx_ratio = energy_expectation / self.optimal_size
        return self._approx_ratio
    
    @property
    def cvar_ratio(self):
        if self._cvar_ratio is None:
            cvar = compute_cvar(self.counts, self.sizes, alpha = self.alpha)
            self._cvar_ratio = cvar / self.optimal_size
        return self._cvar_ratio
    
    @property
    def svs_gibbs_ratio(self):
        if self._svs_gibbs_ratio is None:
            gibbs = compute_gibbs(self.svs_counts, self.svs_sizes, eta = self.eta)
            self._svs_gibbs_ratio = gibbs / self.optimal_size / self.eta
        return self._svs_gibbs_ratio
            
    @property
    def svs_bestcut_ratio(self):
        if self._svs_bestcut_ratio is None:
            best_measured_size = compute_best_cut_from_measured(self.svs_counts, self.svs_sizes)
            self._svs_bestcut_ratio = best_measured_size / self.optimal_size
        return self._svs_bestcut_ratio
    
    @property
    def svs_approx_ratio(self):
        if self._svs_approx_ratio is None:
            energy_expectation = compute_energy_expectation(self.svs_counts, self.svs_sizes)
            self._svs_approx_ratio = energy_expectation / self.optimal_size
        return self._svs_approx_ratio
    
    @property
    def svs_cvar_ratio(self):
        if self._svs_cvar_ratio is None:
            cvar = compute_cvar(self.svs_counts, self.svs_sizes, alpha = self.alpha)
            self._svs_cvar_ratio = cvar / self.optimal_size
        return self._svs_cvar_ratio
        
    def svs_simulate(self):
        """Implement a state vector simulation
        """
        # Create QAOA circuit
        svs_circuit = self.circuit.copy()
        svs_circuit.remove_final_measurements()
        # Obtain probabilities using state vector simulation 
        sv = Aer.get_backend('statevector_simulator')
        svs_counts = execute(svs_circuit, sv).result().get_counts()
        
        self.svs_cuts = list(svs_counts.keys())
        self.svs_counts = list(svs_counts.values())
        self.svs_sizes = [self.eval_cut(solution[::-1]) for solution in self.svs_cuts] # Reverse each cut passed ot eval_cut, since qiskit uses the little-endian format
        
        self.svs_dist_counts, self.svs_dist_sizes = get_size_dist(self.svs_counts, self.svs_sizes, self.optimal_size)
            
    def random_sampling_dist(self):
        # Obtain num_shots number of uniform random samples between 0 and 2 ** nodes
        unif_cuts = np.random.randint(2 ** self.nodes, size=self.num_shots).tolist()
        unif_cuts_uniq = list(set(unif_cuts))

        # Get counts corresponding to each sampled int/cut
        unif_counts = [unif_cuts.count(cut) for cut in unif_cuts_uniq]
        unif_cuts = list(set(unif_cuts))

        def int_to_bs(numb):
            # Function for converting from an integer to list of 0's and 1's
            strr = format(numb, "b") #convert to binary
            strr = '0' * (self.nodes - len(strr)) + strr
            strr = list(strr)
            return strr

        unif_cuts = [int_to_bs(i) for i in unif_cuts]
        unif_sizes = [self.eval_cut(cut) for cut in unif_cuts]

        # Also get the corresponding distribution of cut sizes
        self.unif_dist_counts, self.unif_dist_sizes = get_size_dist(unif_counts, unif_sizes, self.optimal_size)
            
    def draw_graph(self):
        "Draw graph"
        nx.draw(self.nxgraph)
    
    def draw_ansatz(self):
        circuitdiagram = circuit_drawer(self.circuit)
        print(circuitdiagram)
    
    def brute_force_compute_maxcut(self):
        """Compute the max cut size and (one of possibly multiple) max cut.
        Compute classical using a brute force approach, going through all possible cuts.
        """
        list_0_1s = [[0]] + [[0,1]] * (self.nodes-1)
        
        cut_list = list(iterprod(*list_0_1s))
        cut_size_list = [self.eval_cut(cut) for cut in cut_list]
        self.optimal_size = max(cut_size_list)
        self.optimal_cut = cut_list[cut_size_list.index(self.optimal_size)] # Only one optimal cut will be stored
        
    def eval_cut(self, cut):
        """Compute the cut size corresponding to a given solution

        Args:
            cut : List of 0s and 1s, denoting a cut

        Returns:
            int: cut size
        """
        cut_size = len([1 for (i,j) in self.edges if cut[i]!=cut[j]])
        return cut_size
        
    def execute_circuit(self):
        """
        Implement a noiseless simulation of the ansatz circuit.
        Compute the measured cuts, along with the number of times they were measured
        Creates:
            cuts (list): List of bitstrings, denoting measured cuts
            counts (list): List of integers, denoting number of times each cut was measured 
            sizes (list): List of integers. Cut size corresponding to each measured cut
        """
        # Implement a noiseless simulation using the qasm_simulator 
        backend = Aer.get_backend('qasm_simulator')
        cuts_counts_dict = backend.run(self.circuit, shots=self.num_shots).result().get_counts() # Dictionary of ('bitstring', number of times it was measured) pairs
        
        self.cuts = list(cuts_counts_dict.keys())
        self.counts = list(cuts_counts_dict.values())
        self.sizes = [self.eval_cut(solution[::-1]) for solution in self.cuts] # Reverse each cut passed ot eval_cut, since qiskit uses the little-endian format
        
        self.dist_counts, self.dist_sizes = get_size_dist(self.counts, self.sizes, self.optimal_size)
        
        
    def plot_cutsize_dist(self):
        fig, axs = plt.subplots(1, 1)

        suptitle = "Empirical Distribution of Cut Sizes\n Graph Size={}".format(self.nodes)
        plt.title(suptitle)

        axs.plot(np.array(self.dist_sizes) / self.optimal_size, np.array(self.dist_counts) / self.num_shots, marker='o',ls='-', c='k', ms=2, mec='k', mew=0.4, lw=1, label=f"Circuit Sampling")

        # Also plot the distribution obtained from uniform random sampling
        axs.plot(np.array(self.unif_dist_sizes) / self.optimal_size, np.array(self.unif_dist_counts) / self.num_shots,
             marker='o', ms=1, mec = 'k',mew=0.2, lw=10,alpha=0.5,
             ls = '-', label = "Uniform Random Sampling", c = "pink")  # " degree={deg}") # lw=1,

        # # Plot vertical lines corresponding to the various metrics
        # plotted_metric_values = []
        # for metric in ['approx_ratio', 'cvar_ratio', 'bestcut_ratio', 'gibbs_ratio']:
        #     lw=1; ls='solid'
        #     curmetricval = 
        #     if curmetricval in plotted_metric_values:
        #         # for lines that will coincide, assign different styles to distinguish them
        #         lw=1.5; ls='dashed'
        #     plotted_metric_values.append(curmetricval)
        #     axs.axvline(x=curmetricval, color=curdict['color'], label=curdict['label'], lw=lw, ls=ls)
            

        axs.set_ylabel('Fraction of Total Counts')
        axs.set_xlabel(r'$\frac{\mathrm{Cut\ Size}}{\mathrm{Max\ Cut\ Size}}$')
        axs.grid()
        axs.set_xlim(left=-0.02, right=1.02)
        axs.legend(loc='upper left')

        fig.tight_layout()

    

    def create_ansatz_circuit(self):
        """
        Create the ansatz circuit (determined by the ?? and ?? parameters)
        """
        self.circuit = QuantumCircuit(self.nodes)
        for i in range(self.nodes):
            # Apply the Hadamard gate on each qubit
            self.circuit.h(i)

        for roundIndex in range(self.rounds):
            # For each round, apply RZZ gates followed by RX gates
            beta = self.betas[roundIndex]
            gamma = self.gammas[roundIndex]

            for i,j in self.edges:
                # exp(i gamma/2 Z_i Z_j)
                self.circuit.rzz(- gamma, i, j)
                # circuit.rzz( -gamma, i, j)
            for i in range(self.nodes):
                # exp(-i beta X_i)
                self.circuit.rx(2 * beta, i)
                
        self.circuit.measure_all()
        




