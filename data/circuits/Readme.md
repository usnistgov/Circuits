Each file contains the list of representative Boolean functions with an AND optimal implementation. The information is presented in the following way:

 - First line is the algebraic normal form of the Boolean function whose implementation is provided.
 - Second line is the topology of the circuit.
 - Third line is the *output mask* (written in decimal but interpreted in binary) that shows which AND gate outputs are going to be added to form the final output. The least significant bit denotes the AND gate A0, the next bit A1, and so on.
 - The following *2k* lines (*k* is the number of AND gates in the circuit) denotes the linear functions to be provided to the AND gates. The first two linear functions will be used by A0, etc.
 - The last line is an affine function named *final mask* that is going to be added to the output in order to obtain the desired function.

Example:

```
x1x2+x1x2x3x4
(A0 L : L) (A1 L : 0) (A2 L : 0 1)
4
x1+x4
x4
x3+x4
x4
x2
x1+x4
0
```

This example implements the Boolean function `x1x2+x1x2x3x4`. The topology suggests that the function requires three AND gates. The *output mask* is 4, which is `100` in binary, meaning that only the output of A2 will be used. The following six linear functions will be fed into the three AND gates; two linear functions per AND gate. The *final mask* in this case is `0`, so the output is just going to be the output of the last AND gate A2.

This representation can be transformed to a more descriptive format as follows:

```
T = (A0 L : L) (A1 L : 0) (A2 L : 0 1)
a0 = (x4+x1) * (x4)                      # a0 = x4+x1x4
a1 = (x4+x3) * (a0+x4)                   # a1 = x1x4+x1x3x4
a2 = (x2) * (a0+a1+x4+x1)                # a2 = x1x2+x1x2x3x4
f  = a2
```


