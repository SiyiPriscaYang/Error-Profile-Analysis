# Error-Profile-Analysis
This software plots the error vectors associated with decoding failures of QLDPC codes, helping visualize the errors and faciliates code construction and decoder optimization.

## Inputs
The file has two inputs:
- A .txt file containing the parity check matrix of the QLDPC code (nonbinary, 0, 1, 2, 3 for I, X, Y, Z, use spaces to separate entries)
  ```
  H=load('H_3_7_m3.txt');
  ```
- A .txt file containing the error vectors associated with decoding failures obtained from the decoder (nonbinary, 0, 1, 2, 3 for I, X, Y, Z, use spaces to separate entries)
  ```
  E = load('nb_deg_3_7_m3_0.04_500.txt');
  ```
  
## Outputs
