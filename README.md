# Error-Profile-Analysis
This software visualizes the error vectors associated with decoding failures of QLDPC codes, facilitating code construction and decoder optimization.

## Inputs
The program requires two input files:
- Parity check matrix file: A '.txt' file containing the parity check matrix of the QLDPC code (nonbinary format, with 0, 1, 2, 3 representing I, X, Y, and Z, respectively; entries should be separated by spaces).
  ```
  H=load('H_3_7_m3.txt');
  ```
- Error vector file: A `.txt` file containing the error vectors associated with decoding failures obtained from the decoder (nonbinary format, with 0, 1, 2, 3 for I, X, Y, and Z; entries separated by spaces).
  ```
  E = load('nb_deg_3_7_m3_0.04_500.txt');
  ```

## Outputs
For each error vector in the input file, the program plots the neighborhood of the erroneous qubits, which includes:
- Erroneous qubits: Represented by circle nodes, with X, Z, and Y errors highlighted in red, blue, and green, respectively.
- Check nodes (stabilizers): Adjacent to the erroneous qubits, represented by square nodes. X and Z stabilizers are highlighted in red and blue, respectively.
