
from qiskit import QuantumCircuit, Aer, execute
import numpy as np
import math
import networkx as nx
import matplotlib.pyplot as plt
from itertools import product as iterprod

def get_size_dist(counts, sizes):
    """ For given measurement outcomes, i.e. combinations of counts and sizes, return counts corresponding to each cut size.
    Args:
        counts (list): List of integers, denoting number of times each cut was measured 
        sizes (list): List of integers. Cut size corresponding to each measured cut
    Returns:
        unique_counts (list) : Number of times each size was obtained
        unique_sizes (list) : List of all sizes
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
    return unique_counts, unique_sizes
    
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
    return gibbs


def compute_best_cut_from_measured(counts, sizes, **kwargs):
    """
    From the measured cuts, return the size of the largest cut
    """
    return np.max(sizes)
    
class qaoa_anstaz:
    def __init__(self, nodes, edges, betas, gammas, num_shots = 5000, eta = 0.5, optimal_cut = [], optimal_size = -1):
        """
        Initialize a class for a QAOA ansatz circuit. Create corresponding circuit as well.

        Args:
            nodes (int): number of nodes in the graph
            edges (list of tuples): list of (undirected) edges of the graph. Each element is a two-tuple of node labels, which should be between 0 and (nodes-1) 
            betas (list of floats): List of β angles. Number of rounds is equal to the length of this list.
            gammas (list of floats): List of γ angles. Number of rounds is equal to the length of this list.
            num_shots (int): Number of times the qaoa ansatz circuit is to be measured
        """
        self.nodes = nodes
        self.edges = edges
        self.betas = betas
        self.gammas = gammas
        self.num_shots = num_shots
        self.eta = eta
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
    
    def draw_graph(self):
        "Draw graph"
        nx.draw(self.nxgraph)
    
    def draw_ansatz(self):
        self.circuit.draw(output='mpl')
    
    
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
        cuts_counts_dict = backend.run(self.circuit, nshots=self.num_shots).result().get_counts() # Dictionary of ('bitstring', number of times it was measured) pairs
        
        self.cuts = list(cuts_counts_dict.keys())
        self.counts = list(cuts_counts_dict.values())
        self.sizes = [self.eval_cut(solution[::-1]) for solution in self.cuts] # Reverse each cut passed ot eval_cut, since qiskit uses the little-endian format
        
        self.unique_counts, self.unique_sizes = get_size_dist(self.counts, self.sizes)
        
        self.compute_all_metrics()
        
    def compute_all_metrics(self):
        self.best_measured_size = compute_best_cut_from_measured(self.counts, self.sizes)
        self.gibbs = compute_gibbs(self.counts, self.sizes, self.eta)
        self.cvar = compute_cvar(self.counts, self.sizes, self.alpha)
        self.energy_expectation = compute_energy_expectation(self.counts, self.sizes)
        
        self.bestcut_ratio = self.best_measured_size / self.optimal_size
        

    def create_ansatz_circuit(self):
        """
        Create the ansatz circuit (determined by the β and γ parameters)
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
        




