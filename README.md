### Quantum state and channel discrimination

This is a WIP, the code and this README is still incomplete/not optimal. 

Files:
- non_optimal_discrimination.py: Distinguishes depolarizing channel from identity for a qubit and 
qutrit system, along with a depolarizing factor. Doesn't use a choi jamiolkowski matrix, hence 
non-optimal.

- qcd_cj_matrix.py: Uses CJ matrix, and also computes the graph for different p values and probabilities.
Adapted code from https://github.com/vsiddhu/SDP-Quantum-OR/blob/main/Notebook%203%20-%20Quantum%20Channel%20Discrimination.ipynb
