This folder contains data that has been used and produced in the circuit complexity research.


## Descriptions of folders

---
- `ANF' Files with Algebraic Normal Form (ANF) of various polynomials over GF(2)
- `circuits` List of circuits for producing the representatives of Boolean functions with up to 6-variables.
- `classes` The representatives of affine equivalence classes.
- `mc_dim` The representatives of affine equivalence classes grouped w.r.t. the multiplicative complexity and dimension.
- `slp` A different representation of the circuit information in folder `circuits`.
- `topologies` List of topologies up to 5 AND gates.

The file formats are described inside the folders. 


## n6_group_data.txt

The file `n6_group_data.txt` is used as precomputed data to speed up the method of identifying the classes of 6-variable Boolean functions. The file contains the 150357 representative Boolean functions in 30883 groups, where all the files in the same group have the same signature (algebraic degree, distribution of absolute values in the Walsh spectrum and autocorrelation spectrum). Each group is identified by a line declaring its group number, the number of functions in the group, and additional information for equivalence matching if there are more than one functions in the group. This line is followed by the truth tables of the functions in hexadecimal format.

The following excerpt from the file descibes three groups:

```
group 3 size 1
8888888888888888

group 4 size 3 ind 2 
16e6da2abc4c7080
4ed236aae4789c00
b2e42e78ca9c5600

group 5 size 2 lc 2 
be281248d4e87888
eed2b47896aacc00
```

`group 3` has only one function, so no additional information is required for class identification. All functions having the signature of `group 3` are from the same equivalence class.

For groups having more than one Boolean functions, a hint to a class identification method is provided. We want to find in the most efficient way which of the three functions a given function is affine equivalent to, since we cannot distinguish it based on the signaure; there are more than one candidates. This hint is calculated based on the benchmarks of different identification methods and parameters and the fastest one is chosen.

`group 4 has` three functions, so  an identification hint is required, which is `ind 2`, meaning that the *indicators* of all three functions up to distance 2 will be computed and used in the identification process.

`group 5 has` two functions. Similar to the previous case, but a *local connections* distinguisher (`lc 2`) is used this time.

The details of those two methods can be found in the following thesis:

- Joanne Elizabeth Fuller, **Analysis of Affine Equivalent Boolean Functions for Cryptography**,  PhD thesis, Queensland University of Technology, 2003.

