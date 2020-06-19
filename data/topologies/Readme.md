A topology in the context of circuits represents the relations of AND gates between each other. It represents which AND gate takes input from which other AND gates. Each file here contains a list of topologies with a certain number of AND gates. A topology is specified in a single line such as:

`(A0 L : L) (A1 L : 0) (A2 0 : 1)`

Each AND gate is specified between paranthesis with a gate index starting from 0, followed by two sets of inputs separated by a colon. The inputs denote other AND gate indices. If a specific input does not take an AND gate input, it is denoted by **L**, which stands for a linear (or affine) function. For the other inputs that specify indices, the **L** input is implicitly assumed.

A topology definition also satisfies the following conditions:

 - An AND gate can only use previously defined gate indices, there is no forward reference. 
 - The first gate is always of the form `(A0 L : L)` because there is no other AND gate to be used as input.
 
 For more information, see:
 
 -  Ç. Çalık, M. Sönmez Turan, and R. Peralta, The multiplicative complexity of 6-variable Boolean functions. Cryptography and Communications. Special Issue on Boolean Functions and Their Applications, pp. 1–15, 2018. [doi:10.1007/s12095-018-0297-2](https://doi.org/10.1007/s12095-018-0297-2). Also at <ia.cr/2018/002>