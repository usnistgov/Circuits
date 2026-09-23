# AES Circuits

Boolean circuits (straight-line programs) related to the *Advanced Encryption Standard* (AES).

This page is maintained by the NIST [Circuit Complexity program](https://csrc.nist.gov/projects/circuit-complexity), for educational and research purposes.
Refer to NIST [ACVP](https://github.com/usnistgov/ACVP) for validation of AES implementations aimed for production.

**Highlights:** In the various tables, selected columns highlight the lowest (best) value in bold, underlined.

<details open>
<summary><h2>Index of subfolders</h2></summary>

- **Building blocks:** [aes-sbox](aes-sbox/README.md), [aes-invsbox](aes-invsbox/README.md), [aes-mixcols](aes-mixcols/README.md), [aes-invmixcols](aes-invmixcols/README.md).
- **Folded circuits:** See [aes-fold1](aes-fold1/README.md), covering `keyexp`, `cipher`, `encipher`, `invcipher`, and `decipher` (with regard to AES-128, AES-192, and AES-256), without unfolding the building blocks.
- **Flat circuits:** [KeyExpansion](aes-keyexp/README.md), [Cipher](aes-cipher/README.md), [Encipher](aes-encipher/README.md), [InvCipher](aes-invcipher/README.md), [Decipher](aes-decipher/README.md).


</details>
<details open>
<summary><h2>Building blocks: [inv]sbox and [inv]mixcols</h2></summary>

**Building blocks:**
- `sbox` (8-bit to 8-bit): implements `SBox()` from FIPS 197.
- `invsbox` (8-bit to 8-bit): implements `InvSBox()` from FIPS 197.
- `mixcols` (32-bit to 32-bit): implements `MixColumns()` from FIPS 197.
- `invmixcols` (32-bit to 32-bit): implements `InvMixColumns()` from FIPS 197.

Recently added circuits were obtained with optimization techniques developed with assistance from AI. With additional computation, most of these circuits can likely be improved in some way. The circuits can be easily verified as correct based on input/output behavior.

<details open>
<summary><h3>AES S-box</h3></summary>

Example circuits for the AES S-box (8-bit to 8-bit function), for both Forward and Inverse directions.

#### S-box Forward with #AND ≤ 34

| File | #AND<br>(*A*) | AND<br>depth | #Gates<br>(*G*) | Gate<br>depth | XX<br>(*x*+*x'*) | #XOR<br>(*x*) | #XNOR<br>(*x'*) |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ.txt](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [<ins><strong>29</strong></ins>](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | 5 | 139 | 35 | 110 | 90 | 20 |
| [circ.txt](./aes-sbox/aes-sbox-a29-ad5-g184-gd16-xx155-26.circ.txt) | [<ins><strong>29</strong></ins>](./aes-sbox/aes-sbox-a29-ad5-g184-gd16-xx155-26.circ.txt) | 5 | 184 | 16 | 155 | 129 | 26 |
| [circ.txt](./aes-sbox/aes-sbox-a30-ad4-g147-gd36-xx117-13.circ.txt) | 30 | [<ins><strong>4</strong></ins>](./aes-sbox/aes-sbox-a30-ad4-g147-gd36-xx117-13.circ.txt) | 147 | 36 | 117 | 104 | 13 |
| [circ.txt](./aes-sbox/aes-sbox-a30-ad4-g205-gd15-xx175-54.circ.txt) | 30 | [<ins><strong>4</strong></ins>](./aes-sbox/aes-sbox-a30-ad4-g205-gd15-xx175-54.circ.txt) | 205 | [<ins><strong>15</strong></ins>](./aes-sbox/aes-sbox-a30-ad4-g205-gd15-xx175-54.circ.txt) | 175 | 121 | 54 |
| [circ.txt](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 32 | 5 | [<ins><strong>110</strong></ins>](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 23 | 78 | 75 | 3 |
| [circ.txt](./aes-sbox/aes-sbox-a34-ad4-g110-gd22-xx76-3.circ.txt) | 34 | [<ins><strong>4</strong></ins>](./aes-sbox/aes-sbox-a34-ad4-g110-gd22-xx76-3.circ.txt) | [<ins><strong>110</strong></ins>](./aes-sbox/aes-sbox-a34-ad4-g110-gd22-xx76-3.circ.txt) | 22 | [<ins><strong>76</strong></ins>](./aes-sbox/aes-sbox-a34-ad4-g110-gd22-xx76-3.circ.txt) | 73 | 3 |
| [circ.txt](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 34 | [<ins><strong>4</strong></ins>](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 128 | [<ins><strong>15</strong></ins>](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 94 | 90 | 4 |

Note: The circuits with #AND < 32 were obtained from transformations applied to the [AES S-Box circuit](https://github.com/umizame/S-box_29-AND/blob/master/circuits/aes-sbox-fwd-g228-a29-d35-ad6.slp) from [@umizame](https://github.com/umizame/S-box_29-AND), with 29 AND, 195 XOR, 4 NOT, depth 35, and AND-depth 6.

#### S-box Forward with #AND > 34

| File | #AND<br>(*A*) | AND<br>depth | #Gates<br>(*G*) | Gate<br>depth | XX<br>(*x*+*x'*) | #XOR<br>(*x*) | #XNOR<br>(*x'*) |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ.txt](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [<ins><strong>36</strong></ins>](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 4 | [<ins><strong>138</strong></ins>](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 14 | [<ins><strong>102</strong></ins>](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 80 | 22 |
| [circ.txt](./aes-sbox/aes-sbox-a39-ad4-g148-gd13-xx109-24.circ.txt) | 39 | 4 | 148 | 13 | 109 | 85 | 24 |
| [circ.txt](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 45 | 4 | 151 | [<ins><strong>12</strong></ins>](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 106 | 80 | 26 |
| [circ.txt](./aes-sbox/aes-sbox-a46-ad3-g250-gd15-xx204-4.circ.txt) | 46 | [<ins><strong>3</strong></ins>](./aes-sbox/aes-sbox-a46-ad3-g250-gd15-xx204-4.circ.txt) | 250 | 15 | 204 | 200 | 4 |
| [circ.txt](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 47 | [<ins><strong>3</strong></ins>](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 225 | 15 | 178 | 174 | 4 |

</details>
<details open>
<summary><h3>AES Inverse S-box</h3></summary>

#### S-box Inverse with #AND ≤ 34

| File | #AND<br>(*A*) | AND<br>depth | #Gates<br>(*G*) | Gate<br>depth | XX<br>(*x*+*x'*) | #XOR<br>(*x*) | #XNOR<br>(*x'*) |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ.txt](./aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [<ins><strong>29</strong></ins>](./aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | 5 | 145 | 32 | 116 | 87 | 29 |
| [circ.txt](./aes-invsbox/aes-invsbox-a29-ad5-g186-gd16-xx157-43.circ.txt) | [<ins><strong>29</strong></ins>](./aes-invsbox/aes-invsbox-a29-ad5-g186-gd16-xx157-43.circ.txt) | 5 | 186 | 16 | 157 | 114 | 43 |
| [circ.txt](./aes-invsbox/aes-invsbox-a30-ad4-g150-gd30-xx120-33.circ.txt) | 30 | [<ins><strong>4</strong></ins>](./aes-invsbox/aes-invsbox-a30-ad4-g150-gd30-xx120-33.circ.txt) | 150 | 30 | 120 | 87 | 33 |
| [circ.txt](./aes-invsbox/aes-invsbox-a30-ad4-g203-gd15-xx173-52.circ.txt) | 30 | [<ins><strong>4</strong></ins>](./aes-invsbox/aes-invsbox-a30-ad4-g203-gd15-xx173-52.circ.txt) | 203 | [<ins><strong>15</strong></ins>](./aes-invsbox/aes-invsbox-a30-ad4-g203-gd15-xx173-52.circ.txt) | 173 | 121 | 52 |
| [circ.txt](./aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | 32 | 5 | [<ins><strong>112</strong></ins>](./aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | 27 | [<ins><strong>80</strong></ins>](./aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | 71 | 9 |
| [circ.txt](./aes-invsbox/aes-invsbox-a34-ad4-g114-gd27-xx80-11.circ.txt) | 34 | [<ins><strong>4</strong></ins>](./aes-invsbox/aes-invsbox-a34-ad4-g114-gd27-xx80-11.circ.txt) | 114 | 27 | [<ins><strong>80</strong></ins>](./aes-invsbox/aes-invsbox-a34-ad4-g114-gd27-xx80-11.circ.txt) | 69 | 11 |
| [circ.txt](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | 34 | [<ins><strong>4</strong></ins>](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | 134 | [<ins><strong>15</strong></ins>](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | 100 | 82 | 18 |

Note: The circuits with #AND < 32 were obtained from transformations applied to the [AES S-Box circuit](https://github.com/umizame/S-box_29-AND/blob/master/circuits/aes-sbox-fwd-g228-a29-d35-ad6.slp) from [@umizame](https://github.com/umizame/S-box_29-AND), with 29 AND, 195 XOR, 4 NOT, depth 35, and AND-depth 6.

#### S-box Inverse with #AND > 34

| File | #AND<br>(*A*) | AND<br>depth | #Gates<br>(*G*) | Gate<br>depth | XX<br>(*x*+*x'*) | #XOR<br>(*x*) | #XNOR<br>(*x'*) |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ.txt](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | 36 | 4 | 147 | 14 | 111 | 87 | 24 |

</details>
<details open>
<summary><h3>[Inv]S-Box circuits added around 2020</h3></summary>

#### S-box circuits added around 2020

The following historical table was retrieved from an old version of the [NIST Circuit Complexity list of circuits](https://csrc.nist.gov/Projects/circuit-complexity/list-of-circuits).

| File | Direction | #AND<br>(*A*) | #Gates<br>(*G*) | Gate<br>depth | XX<br>(*x*+*x'*) | #XOR<br>(*x*) | #XNOR<br>(*x'*) |
|---|---|---:|---:|---:|---:|---:|---:|
| [slp](./old-2020/aes-sbox-fwd-g115-a32-d28-ad6.slp) | Forward | 32 | 115 | 28 | 83 | 79 | 4 |
| [slp](./old-2020/aes-sbox-fwd-g113-a32-d27-ad6.slp) | Forward | 32 | 113 | 27 | 81 | 77 | 4 |
| [slp](./old-2020/aes-sbox-fwd-g128-a34-d16-ad4.slp) | Forward | 34 | 128 | 16 | 94 | 90 | 4 |
| [slp](./old-2020/aes-sbox-rev-g121-a34-d21-ad4.slp) | Inverse | 34 | 121 | 21 | 87 | 83 | 4 |
| [slp](./old-2020/aes-sbox-rev-g127-a34-d16-ad4.slp) | Inverse | 34 | 127 | 16 | 93 | 83 | 10 |

</details>
<details open>
<summary><h3>AES MixColumns and InvMixColumns</h3></summary>

Example circuits for AES MixColumns (32-bit to 32-bit linear functions), for both Forward and Inverse directions.

Note: For each direction, dark-blue, bold, underlined values indicate the lowest depth and lowest #XOR within the displayed selection.

| File | Direction | Depth | #XOR |
|---|---|---:|---:|
| [circ.txt](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | Forward | [<ins><strong>3</strong></ins>](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 97 |
| [circ.txt](./aes-mixcols/aes-mixcols-xor90-depth4.circ.txt) | Forward | 4 | 90 |
| [circ.txt](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | Forward | 5 | [<ins><strong>88</strong></ins>](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) |
| [circ.txt](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | Inverse | [<ins><strong>5</strong></ins>](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 146 |
| [circ.txt](./aes-invmixcols/aes-invmixcols-xor135-depth6.circ.txt) | Inverse | 6 | 135 |
| [circ.txt](./aes-invmixcols/aes-invmixcols-xor124-depth7.circ.txt) | Inverse | 7 | 124 |
| [circ.txt](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | Inverse | 8 | [<ins><strong>114</strong></ins>](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) |




</details>
</details>
<details open>
<summary><h2>Nomenclature for large AES circuits</h2></summary>

- **Key-expansion**
  - `keyexp`: implements `KeyExpansion()` from FIPS 197. Note that both `encipher` and `decipher` integrate a key expansion within them.
- **From plaintext to ciphertext**
  - `cipher`: implements `Cipher()` from FIPS 197; one of its inputs is an expanded key.
  - `encipher`: implements the `AES-128()`, `AES-192()`, and `AES-256()` from FIPS 197. One of the inputs is a *non*-expanded key.
- **From ciphertext to plaintext**
  - `invcipher`: implements `InvCipher()` from FIPS 197. The input includes the original key and an expanded key.
  - `decipher`: One of the inputs is a *non*-expanded key. This integration of `keyexp` and `invcipher` is not defined as a function in FIPS 197, but is here for convenience (`decipher` is to `invcipher` as `encipher` is to `cipher`).

</details>
<details open>
<summary><h2>AES folded circuits</h2></summary>

A circuit is "folded" when some components are not "flattened" to a sequence of basic gates. The circuits below call `aes-sbox`, `aes-invsbox`, `aes-mixcols`, or `aes-invmixcols` as applicable, and use vector operations (`VXOR` and `VXNOR`) to describe XOR-family gates succinctly.

The table counts top-level operations as written in each folded circuit; #XOR and #XNOR give the numbers of scalar XOR and XNOR gates, respectively, represented by the VXOR and VXNOR operations.

| Operation | File | Key<br>size | #sbox | #inv<br>sbox | #mixcols | #inv<br>mixcols | #VXOR | #VXNOR | #XOR | #XNOR |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| KeyExpansion | [aes128-keyexp](./aes-fold1/aes128-keyexp-fold1.circ.txt) | 128 | 40 | — | — | — | 40 | 2 | 1264 | 8 |
| KeyExpansion | [aes192-keyexp](./aes-fold1/aes192-keyexp-fold1.circ.txt) | 192 | 32 | — | — | — | 46 | — | 1464 | 8 |
| KeyExpansion | [aes256-keyexp](./aes-fold1/aes256-keyexp-fold1.circ.txt) | 256 | 52 | — | — | — | 52 | — | 1657 | 7 |
| Cipher | [aes128-cipher](./aes-fold1/aes128-cipher-fold1.circ.txt) | 128 | 160 | — | 36 | — | 11 | — | 1408 | — |
| Cipher | [aes192-cipher](./aes-fold1/aes192-cipher-fold1.circ.txt) | 192 | 192 | — | 44 | — | 14 | — | 1664 | — |
| Cipher | [aes256-cipher](./aes-fold1/aes256-cipher-fold1.circ.txt) | 256 | 224 | — | 52 | — | 15 | — | 1920 | — |
| Encipher | [aes128-encipher](./aes-fold1/aes128-encipher-fold1.circ.txt) | 128 | 200 | — | 36 | — | 51 | 2 | 2672 | 8 |
| Encipher | [aes192-encipher](./aes-fold1/aes192-encipher-fold1.circ.txt) | 192 | 224 | — | 44 | — | 60 | — | 3128 | 8 |
| Encipher | [aes256-encipher](./aes-fold1/aes256-encipher-fold1.circ.txt) | 256 | 276 | — | 52 | — | 67 | — | 3577 | 7 |
| InvCipher | [aes128-invcipher](./aes-fold1/aes128-invcipher-fold1.circ.txt) | 128 | — | 160 | — | 36 | 11 | — | 1408 | — |
| InvCipher | [aes192-invcipher](./aes-fold1/aes192-invcipher-fold1.circ.txt) | 192 | — | 192 | — | 44 | 14 | — | 1664 | — |
| InvCipher | [aes256-invcipher](./aes-fold1/aes256-invcipher-fold1.circ.txt) | 256 | — | 224 | — | 52 | 15 | — | 1920 | — |
| Decipher | [aes128-decipher](./aes-fold1/aes128-decipher-fold1.circ.txt) | 128 | 40 | 160 | — | 36 | 51 | 2 | 2672 | 8 |
| Decipher | [aes192-decipher](./aes-fold1/aes192-decipher-fold1.circ.txt) | 192 | 32 | 192 | — | 44 | 60 | — | 3128 | 8 |
| Decipher | [aes256-decipher](./aes-fold1/aes256-decipher-fold1.circ.txt) | 256 | 52 | 224 | — | 52 | 67 | — | 3577 | 7 |


</details>
<details open>
<summary><h2>AES flat circuits</h2></summary>

Example circuits, flattened to the level of basic Boolean gates (AND, XOR, XNOR).

There are up to 4! = 24 tuple metrics corresponding to the possible orderings of (A, AD, G, GD). For succinctness, the tables consider only four tuple metrics (TM): (1) A-AD-G-GD, (2) AD-GD-G-A, (3) GD-G-AD-A, and (4) G-A-GD-AD.

Note: In selected columns, dark-blue, bold, underlined values indicate the lowest displayed value.


<details open>
<summary><h3>AES-128 KeyExpansion: Flat circuits</h3></summary>

#### KeyExpansion using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| [circ.txt](./aes-keyexp/aes128-keyexp-a1160-ad50-g6840-gd381-xx5680-816.circ.txt) | 128 | [<ins><strong>1160</strong></ins>](./aes-keyexp/aes128-keyexp-a1160-ad50-g6840-gd381-xx5680-816.circ.txt) | 50 | 6840 | 381 | 5680 | 4864 | 816 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | 1 |
| [circ.txt](./aes-keyexp/aes128-keyexp-a1360-ad40-g6400-gd190-xx5040-176.circ.txt) | 128 | 1360 | [<ins><strong>40</strong></ins>](./aes-keyexp/aes128-keyexp-a1360-ad40-g6400-gd190-xx5040-176.circ.txt) | 6400 | [<ins><strong>190</strong></ins>](./aes-keyexp/aes128-keyexp-a1360-ad40-g6400-gd190-xx5040-176.circ.txt) | 5040 | 4864 | 176 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 2,3 |
| [circ.txt](./aes-keyexp/aes128-keyexp-a1280-ad50-g5680-gd270-xx4400-136.circ.txt) | 128 | 1280 | 50 | [<ins><strong>5680</strong></ins>](./aes-keyexp/aes128-keyexp-a1280-ad50-g5680-gd270-xx4400-136.circ.txt) | 270 | 4400 | 4264 | 136 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### KeyExpansion using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| [circ.txt](./aes-keyexp/aes128-keyexp-a1440-ad40-g6800-gd180-xx5360-896.circ.txt) | 128 | [<ins><strong>1440</strong></ins>](./aes-keyexp/aes128-keyexp-a1440-ad40-g6800-gd180-xx5360-896.circ.txt) | 40 | [<ins><strong>6800</strong></ins>](./aes-keyexp/aes128-keyexp-a1440-ad40-g6800-gd180-xx5360-896.circ.txt) | 180 | 5360 | 4464 | 896 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 1,4 |
| [circ.txt](./aes-keyexp/aes128-keyexp-a1880-ad30-g10280-gd190-xx8400-176.circ.txt) | 128 | 1880 | [<ins><strong>30</strong></ins>](./aes-keyexp/aes128-keyexp-a1880-ad30-g10280-gd190-xx8400-176.circ.txt) | 10280 | 190 | 8400 | 8224 | 176 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 2 |
| [circ.txt](./aes-keyexp/aes128-keyexp-a1800-ad40-g7320-gd160-xx5520-1056.circ.txt) | 128 | 1800 | 40 | 7320 | [<ins><strong>160</strong></ins>](./aes-keyexp/aes128-keyexp-a1800-ad40-g7320-gd160-xx5520-1056.circ.txt) | 5520 | 4464 | 1056 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-128 Cipher: Flat circuits</h3></summary>

#### Cipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-cipher/aes128-cipher-a4640-ad50-g26816-gd388-xx22176-3200.circ.txt) | 128 | [<ins><strong>4640</strong></ins>](./aes-cipher/aes128-cipher-a4640-ad50-g26816-gd388-xx22176-3200.circ.txt) | 50 | 26816 | 388 | 22176 | 18976 | 3200 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ.txt](./aes-cipher/aes128-cipher-a5440-ad40-g25380-gd188-xx19940-640.circ.txt) | 128 | 5440 | [<ins><strong>40</strong></ins>](./aes-cipher/aes128-cipher-a5440-ad40-g25380-gd188-xx19940-640.circ.txt) | 25380 | [<ins><strong>188</strong></ins>](./aes-cipher/aes128-cipher-a5440-ad40-g25380-gd188-xx19940-640.circ.txt) | 19940 | 19300 | 640 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ.txt](./aes-cipher/aes128-cipher-a5120-ad50-g22176-gd286-xx17056-480.circ.txt) | 128 | 5120 | 50 | [<ins><strong>22176</strong></ins>](./aes-cipher/aes128-cipher-a5120-ad50-g22176-gd286-xx17056-480.circ.txt) | 286 | 17056 | 16576 | 480 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Cipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-cipher/aes128-cipher-a5760-ad40-g26656-gd196-xx20896-3520.circ.txt) | 128 | [<ins><strong>5760</strong></ins>](./aes-cipher/aes128-cipher-a5760-ad40-g26656-gd196-xx20896-3520.circ.txt) | 40 | [<ins><strong>26656</strong></ins>](./aes-cipher/aes128-cipher-a5760-ad40-g26656-gd196-xx20896-3520.circ.txt) | 196 | 20896 | 17376 | 3520 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1,4 |
| [circ.txt](./aes-cipher/aes128-cipher-a7520-ad30-g40900-gd188-xx33380-640.circ.txt) | 128 | 7520 | [<ins><strong>30</strong></ins>](./aes-cipher/aes128-cipher-a7520-ad30-g40900-gd188-xx33380-640.circ.txt) | 40900 | 188 | 33380 | 32740 | 640 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ.txt](./aes-cipher/aes128-cipher-a7200-ad40-g29060-gd158-xx21860-4160.circ.txt) | 128 | 7200 | 40 | 29060 | [<ins><strong>158</strong></ins>](./aes-cipher/aes128-cipher-a7200-ad40-g29060-gd158-xx21860-4160.circ.txt) | 21860 | 17700 | 4160 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-128 Encipher: Flat circuits</h3></summary>

#### Encipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-encipher/aes128-encipher-a5800-ad50-g33656-gd388-xx27856-4016.circ.txt) | 128 | [<ins><strong>5800</strong></ins>](./aes-encipher/aes128-encipher-a5800-ad50-g33656-gd388-xx27856-4016.circ.txt) | 50 | 33656 | 388 | 27856 | 23840 | 4016 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ.txt](./aes-encipher/aes128-encipher-a6800-ad40-g31780-gd191-xx24980-816.circ.txt) | 128 | 6800 | [<ins><strong>40</strong></ins>](./aes-encipher/aes128-encipher-a6800-ad40-g31780-gd191-xx24980-816.circ.txt) | 31780 | [<ins><strong>191</strong></ins>](./aes-encipher/aes128-encipher-a6800-ad40-g31780-gd191-xx24980-816.circ.txt) | 24980 | 24164 | 816 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ.txt](./aes-encipher/aes128-encipher-a6400-ad50-g27856-gd286-xx21456-616.circ.txt) | 128 | 6400 | 50 | [<ins><strong>27856</strong></ins>](./aes-encipher/aes128-encipher-a6400-ad50-g27856-gd286-xx21456-616.circ.txt) | 286 | 21456 | 20840 | 616 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Encipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-encipher/aes128-encipher-a7200-ad40-g33456-gd196-xx26256-4416.circ.txt) | 128 | [<ins><strong>7200</strong></ins>](./aes-encipher/aes128-encipher-a7200-ad40-g33456-gd196-xx26256-4416.circ.txt) | 40 | [<ins><strong>33456</strong></ins>](./aes-encipher/aes128-encipher-a7200-ad40-g33456-gd196-xx26256-4416.circ.txt) | 196 | 26256 | 21840 | 4416 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1,4 |
| [circ.txt](./aes-encipher/aes128-encipher-a9400-ad30-g51180-gd191-xx41780-816.circ.txt) | 128 | 9400 | [<ins><strong>30</strong></ins>](./aes-encipher/aes128-encipher-a9400-ad30-g51180-gd191-xx41780-816.circ.txt) | 51180 | 191 | 41780 | 40964 | 816 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ.txt](./aes-encipher/aes128-encipher-a9000-ad40-g36380-gd161-xx27380-5216.circ.txt) | 128 | 9000 | 40 | 36380 | [<ins><strong>161</strong></ins>](./aes-encipher/aes128-encipher-a9000-ad40-g36380-gd161-xx27380-5216.circ.txt) | 27380 | 22164 | 5216 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-128 InvCipher: Flat circuits</h3></summary>

#### InvCipher using InvSbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-invcipher/aes128-invcipher-a4640-ad50-g28712-gd403-xx24072-4640.circ.txt) | 128 | [<ins><strong>4640</strong></ins>](./aes-invcipher/aes128-invcipher-a4640-ad50-g28712-gd403-xx24072-4640.circ.txt) | 50 | 28712 | 403 | 24072 | 19432 | 4640 | [Link](./aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ.txt](./aes-invcipher/aes128-invcipher-a5440-ad40-g28104-gd206-xx22664-2880.circ.txt) | 128 | 5440 | [<ins><strong>40</strong></ins>](./aes-invcipher/aes128-invcipher-a5440-ad40-g28104-gd206-xx22664-2880.circ.txt) | 28104 | [<ins><strong>206</strong></ins>](./aes-invcipher/aes128-invcipher-a5440-ad40-g28104-gd206-xx22664-2880.circ.txt) | 22664 | 19784 | 2880 | [Link](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ.txt](./aes-invcipher/aes128-invcipher-a5120-ad50-g23432-gd353-xx18312-1440.circ.txt) | 128 | 5120 | 50 | [<ins><strong>23432</strong></ins>](./aes-invcipher/aes128-invcipher-a5120-ad50-g23432-gd353-xx18312-1440.circ.txt) | 353 | 18312 | 16872 | 1440 | [Link](./aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### InvCipher using InvSbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-invcipher/aes128-invcipher-a5760-ad40-g29032-gd223-xx23272-3840.circ.txt) | 128 | [<ins><strong>5760</strong></ins>](./aes-invcipher/aes128-invcipher-a5760-ad40-g29032-gd223-xx23272-3840.circ.txt) | [<ins><strong>40</strong></ins>](./aes-invcipher/aes128-invcipher-a5760-ad40-g29032-gd223-xx23272-3840.circ.txt) | [<ins><strong>29032</strong></ins>](./aes-invcipher/aes128-invcipher-a5760-ad40-g29032-gd223-xx23272-3840.circ.txt) | 223 | 23272 | 19432 | 3840 | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ.txt](./aes-invcipher/aes128-invcipher-a5760-ad40-g30184-gd196-xx24424-3840.circ.txt) | 128 | [<ins><strong>5760</strong></ins>](./aes-invcipher/aes128-invcipher-a5760-ad40-g30184-gd196-xx24424-3840.circ.txt) | [<ins><strong>40</strong></ins>](./aes-invcipher/aes128-invcipher-a5760-ad40-g30184-gd196-xx24424-3840.circ.txt) | 30184 | [<ins><strong>196</strong></ins>](./aes-invcipher/aes128-invcipher-a5760-ad40-g30184-gd196-xx24424-3840.circ.txt) | 24424 | 20584 | 3840 | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-128 Decipher: Flat circuits</h3></summary>

#### Decipher using [Inv]Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|
| [circ.txt](./aes-decipher/aes128-decipher-a5800-ad100-g35552-gd784-xx29752-5456.circ.txt) | 128 | [<ins><strong>5800</strong></ins>](./aes-decipher/aes128-decipher-a5800-ad100-g35552-gd784-xx29752-5456.circ.txt) | 100 | 35552 | 784 | 29752 | 24296 | 5456 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ.txt](./aes-decipher/aes128-decipher-a6800-ad80-g34504-gd396-xx27704-3056.circ.txt) | 128 | 6800 | [<ins><strong>80</strong></ins>](./aes-decipher/aes128-decipher-a6800-ad80-g34504-gd396-xx27704-3056.circ.txt) | 34504 | [<ins><strong>396</strong></ins>](./aes-decipher/aes128-decipher-a6800-ad80-g34504-gd396-xx27704-3056.circ.txt) | 27704 | 24648 | 3056 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ.txt](./aes-decipher/aes128-decipher-a6400-ad100-g29112-gd623-xx22712-1576.circ.txt) | 128 | 6400 | 100 | [<ins><strong>29112</strong></ins>](./aes-decipher/aes128-decipher-a6400-ad100-g29112-gd623-xx22712-1576.circ.txt) | 623 | 22712 | 21136 | 1576 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Decipher using [Inv]Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|
| [circ.txt](./aes-decipher/aes128-decipher-a7200-ad80-g35832-gd403-xx28632-4736.circ.txt) | 128 | [<ins><strong>7200</strong></ins>](./aes-decipher/aes128-decipher-a7200-ad80-g35832-gd403-xx28632-4736.circ.txt) | 80 | [<ins><strong>35832</strong></ins>](./aes-decipher/aes128-decipher-a7200-ad80-g35832-gd403-xx28632-4736.circ.txt) | 403 | 28632 | 23896 | 4736 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ.txt](./aes-decipher/aes128-decipher-a7640-ad70-g40464-gd386-xx32824-4016.circ.txt) | 128 | 7640 | [<ins><strong>70</strong></ins>](./aes-decipher/aes128-decipher-a7640-ad70-g40464-gd386-xx32824-4016.circ.txt) | 40464 | 386 | 32824 | 28808 | 4016 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ.txt](./aes-decipher/aes128-decipher-a7560-ad80-g37504-gd356-xx29944-4896.circ.txt) | 128 | 7560 | 80 | 37504 | [<ins><strong>356</strong></ins>](./aes-decipher/aes128-decipher-a7560-ad80-g37504-gd356-xx29944-4896.circ.txt) | 29944 | 25048 | 4896 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-192 KeyExpansion: Flat circuits</h3></summary>

#### KeyExpansion using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| [circ.txt](./aes-keyexp/aes192-keyexp-a928-ad40-g5920-gd319-xx4992-648.circ.txt) | 192 | [<ins><strong>928</strong></ins>](./aes-keyexp/aes192-keyexp-a928-ad40-g5920-gd319-xx4992-648.circ.txt) | 40 | 5920 | 319 | 4992 | 4344 | 648 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | 1 |
| [circ.txt](./aes-keyexp/aes192-keyexp-a1088-ad32-g5568-gd166-xx4480-136.circ.txt) | 192 | 1088 | [<ins><strong>32</strong></ins>](./aes-keyexp/aes192-keyexp-a1088-ad32-g5568-gd166-xx4480-136.circ.txt) | 5568 | [<ins><strong>166</strong></ins>](./aes-keyexp/aes192-keyexp-a1088-ad32-g5568-gd166-xx4480-136.circ.txt) | 4480 | 4344 | 136 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 2,3 |
| [circ.txt](./aes-keyexp/aes192-keyexp-a1024-ad40-g4992-gd230-xx3968-104.circ.txt) | 192 | 1024 | 40 | [<ins><strong>4992</strong></ins>](./aes-keyexp/aes192-keyexp-a1024-ad40-g4992-gd230-xx3968-104.circ.txt) | 230 | 3968 | 3864 | 104 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### KeyExpansion using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| [circ.txt](./aes-keyexp/aes192-keyexp-a1152-ad32-g5888-gd158-xx4736-712.circ.txt) | 192 | [<ins><strong>1152</strong></ins>](./aes-keyexp/aes192-keyexp-a1152-ad32-g5888-gd158-xx4736-712.circ.txt) | 32 | [<ins><strong>5888</strong></ins>](./aes-keyexp/aes192-keyexp-a1152-ad32-g5888-gd158-xx4736-712.circ.txt) | 158 | 4736 | 4024 | 712 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 1,4 |
| [circ.txt](./aes-keyexp/aes192-keyexp-a1504-ad24-g8672-gd166-xx7168-136.circ.txt) | 192 | 1504 | [<ins><strong>24</strong></ins>](./aes-keyexp/aes192-keyexp-a1504-ad24-g8672-gd166-xx7168-136.circ.txt) | 8672 | 166 | 7168 | 7032 | 136 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 2 |
| [circ.txt](./aes-keyexp/aes192-keyexp-a1440-ad32-g6304-gd142-xx4864-840.circ.txt) | 192 | 1440 | 32 | 6304 | [<ins><strong>142</strong></ins>](./aes-keyexp/aes192-keyexp-a1440-ad32-g6304-gd142-xx4864-840.circ.txt) | 4864 | 4024 | 840 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-192 Cipher: Flat circuits</h3></summary>

#### Cipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-cipher/aes192-cipher-a5568-ad60-g32224-gd466-xx26656-3840.circ.txt) | 192 | [<ins><strong>5568</strong></ins>](./aes-cipher/aes192-cipher-a5568-ad60-g32224-gd466-xx26656-3840.circ.txt) | 60 | 32224 | 466 | 26656 | 22816 | 3840 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ.txt](./aes-cipher/aes192-cipher-a6528-ad48-g30508-gd226-xx23980-768.circ.txt) | 192 | 6528 | [<ins><strong>48</strong></ins>](./aes-cipher/aes192-cipher-a6528-ad48-g30508-gd226-xx23980-768.circ.txt) | 30508 | [<ins><strong>226</strong></ins>](./aes-cipher/aes192-cipher-a6528-ad48-g30508-gd226-xx23980-768.circ.txt) | 23980 | 23212 | 768 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ.txt](./aes-cipher/aes192-cipher-a6144-ad60-g26656-gd344-xx20512-576.circ.txt) | 192 | 6144 | 60 | [<ins><strong>26656</strong></ins>](./aes-cipher/aes192-cipher-a6144-ad60-g26656-gd344-xx20512-576.circ.txt) | 344 | 20512 | 19936 | 576 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Cipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-cipher/aes192-cipher-a6912-ad48-g32032-gd236-xx25120-4224.circ.txt) | 192 | [<ins><strong>6912</strong></ins>](./aes-cipher/aes192-cipher-a6912-ad48-g32032-gd236-xx25120-4224.circ.txt) | 48 | [<ins><strong>32032</strong></ins>](./aes-cipher/aes192-cipher-a6912-ad48-g32032-gd236-xx25120-4224.circ.txt) | 236 | 25120 | 20896 | 4224 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1,4 |
| [circ.txt](./aes-cipher/aes192-cipher-a9024-ad36-g49132-gd226-xx40108-768.circ.txt) | 192 | 9024 | [<ins><strong>36</strong></ins>](./aes-cipher/aes192-cipher-a9024-ad36-g49132-gd226-xx40108-768.circ.txt) | 49132 | 226 | 40108 | 39340 | 768 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ.txt](./aes-cipher/aes192-cipher-a8640-ad48-g34924-gd190-xx26284-4992.circ.txt) | 192 | 8640 | 48 | 34924 | [<ins><strong>190</strong></ins>](./aes-cipher/aes192-cipher-a8640-ad48-g34924-gd190-xx26284-4992.circ.txt) | 26284 | 21292 | 4992 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-192 Encipher: Flat circuits</h3></summary>

#### Encipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-encipher/aes192-encipher-a6496-ad60-g38144-gd466-xx31648-4488.circ.txt) | 192 | [<ins><strong>6496</strong></ins>](./aes-encipher/aes192-encipher-a6496-ad60-g38144-gd466-xx31648-4488.circ.txt) | 60 | 38144 | 466 | 31648 | 27160 | 4488 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ.txt](./aes-encipher/aes192-encipher-a7616-ad48-g36076-gd226-xx28460-904.circ.txt) | 192 | 7616 | [<ins><strong>48</strong></ins>](./aes-encipher/aes192-encipher-a7616-ad48-g36076-gd226-xx28460-904.circ.txt) | 36076 | [<ins><strong>226</strong></ins>](./aes-encipher/aes192-encipher-a7616-ad48-g36076-gd226-xx28460-904.circ.txt) | 28460 | 27556 | 904 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ.txt](./aes-encipher/aes192-encipher-a7168-ad60-g31648-gd344-xx24480-680.circ.txt) | 192 | 7168 | 60 | [<ins><strong>31648</strong></ins>](./aes-encipher/aes192-encipher-a7168-ad60-g31648-gd344-xx24480-680.circ.txt) | 344 | 24480 | 23800 | 680 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Encipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-encipher/aes192-encipher-a8064-ad48-g37920-gd236-xx29856-4936.circ.txt) | 192 | [<ins><strong>8064</strong></ins>](./aes-encipher/aes192-encipher-a8064-ad48-g37920-gd236-xx29856-4936.circ.txt) | 48 | [<ins><strong>37920</strong></ins>](./aes-encipher/aes192-encipher-a8064-ad48-g37920-gd236-xx29856-4936.circ.txt) | 236 | 29856 | 24920 | 4936 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1,4 |
| [circ.txt](./aes-encipher/aes192-encipher-a10528-ad36-g57804-gd226-xx47276-904.circ.txt) | 192 | 10528 | [<ins><strong>36</strong></ins>](./aes-encipher/aes192-encipher-a10528-ad36-g57804-gd226-xx47276-904.circ.txt) | 57804 | 226 | 47276 | 46372 | 904 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ.txt](./aes-encipher/aes192-encipher-a10080-ad48-g41228-gd190-xx31148-5832.circ.txt) | 192 | 10080 | 48 | 41228 | [<ins><strong>190</strong></ins>](./aes-encipher/aes192-encipher-a10080-ad48-g41228-gd190-xx31148-5832.circ.txt) | 31148 | 25316 | 5832 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-192 InvCipher: Flat circuits</h3></summary>

#### InvCipher using InvSbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-invcipher/aes192-invcipher-a5568-ad60-g34520-gd485-xx28952-5568.circ.txt) | 192 | [<ins><strong>5568</strong></ins>](./aes-invcipher/aes192-invcipher-a5568-ad60-g34520-gd485-xx28952-5568.circ.txt) | 60 | 34520 | 485 | 28952 | 23384 | 5568 | [Link](./aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ.txt](./aes-invcipher/aes192-invcipher-a6528-ad48-g33816-gd248-xx27288-3456.circ.txt) | 192 | 6528 | [<ins><strong>48</strong></ins>](./aes-invcipher/aes192-invcipher-a6528-ad48-g33816-gd248-xx27288-3456.circ.txt) | 33816 | [<ins><strong>248</strong></ins>](./aes-invcipher/aes192-invcipher-a6528-ad48-g33816-gd248-xx27288-3456.circ.txt) | 27288 | 23832 | 3456 | [Link](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ.txt](./aes-invcipher/aes192-invcipher-a6144-ad60-g28184-gd425-xx22040-1728.circ.txt) | 192 | 6144 | 60 | [<ins><strong>28184</strong></ins>](./aes-invcipher/aes192-invcipher-a6144-ad60-g28184-gd425-xx22040-1728.circ.txt) | 425 | 22040 | 20312 | 1728 | [Link](./aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### InvCipher using InvSbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-invcipher/aes192-invcipher-a6912-ad48-g34904-gd269-xx27992-4608.circ.txt) | 192 | [<ins><strong>6912</strong></ins>](./aes-invcipher/aes192-invcipher-a6912-ad48-g34904-gd269-xx27992-4608.circ.txt) | [<ins><strong>48</strong></ins>](./aes-invcipher/aes192-invcipher-a6912-ad48-g34904-gd269-xx27992-4608.circ.txt) | [<ins><strong>34904</strong></ins>](./aes-invcipher/aes192-invcipher-a6912-ad48-g34904-gd269-xx27992-4608.circ.txt) | 269 | 27992 | 23384 | 4608 | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ.txt](./aes-invcipher/aes192-invcipher-a6912-ad48-g36312-gd236-xx29400-4608.circ.txt) | 192 | [<ins><strong>6912</strong></ins>](./aes-invcipher/aes192-invcipher-a6912-ad48-g36312-gd236-xx29400-4608.circ.txt) | [<ins><strong>48</strong></ins>](./aes-invcipher/aes192-invcipher-a6912-ad48-g36312-gd236-xx29400-4608.circ.txt) | 36312 | [<ins><strong>236</strong></ins>](./aes-invcipher/aes192-invcipher-a6912-ad48-g36312-gd236-xx29400-4608.circ.txt) | 29400 | 24792 | 4608 | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-192 Decipher: Flat circuits</h3></summary>

#### Decipher using [Inv]Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|
| [circ.txt](./aes-decipher/aes192-decipher-a6496-ad100-g40440-gd804-xx33944-6216.circ.txt) | 192 | [<ins><strong>6496</strong></ins>](./aes-decipher/aes192-decipher-a6496-ad100-g40440-gd804-xx33944-6216.circ.txt) | 100 | 40440 | 804 | 33944 | 27728 | 6216 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ.txt](./aes-decipher/aes192-decipher-a7616-ad80-g39384-gd414-xx31768-3592.circ.txt) | 192 | 7616 | [<ins><strong>80</strong></ins>](./aes-decipher/aes192-decipher-a7616-ad80-g39384-gd414-xx31768-3592.circ.txt) | 39384 | [<ins><strong>414</strong></ins>](./aes-decipher/aes192-decipher-a7616-ad80-g39384-gd414-xx31768-3592.circ.txt) | 31768 | 28176 | 3592 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ.txt](./aes-decipher/aes192-decipher-a7168-ad100-g33176-gd655-xx26008-1832.circ.txt) | 192 | 7168 | 100 | [<ins><strong>33176</strong></ins>](./aes-decipher/aes192-decipher-a7168-ad100-g33176-gd655-xx26008-1832.circ.txt) | 655 | 26008 | 24176 | 1832 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Decipher using [Inv]Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|
| [circ.txt](./aes-decipher/aes192-decipher-a8064-ad80-g40792-gd427-xx32728-5320.circ.txt) | 192 | [<ins><strong>8064</strong></ins>](./aes-decipher/aes192-decipher-a8064-ad80-g40792-gd427-xx32728-5320.circ.txt) | 80 | [<ins><strong>40792</strong></ins>](./aes-decipher/aes192-decipher-a8064-ad80-g40792-gd427-xx32728-5320.circ.txt) | 427 | 32728 | 27408 | 5320 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ.txt](./aes-decipher/aes192-decipher-a8416-ad72-g44984-gd402-xx36568-4744.circ.txt) | 192 | 8416 | [<ins><strong>72</strong></ins>](./aes-decipher/aes192-decipher-a8416-ad72-g44984-gd402-xx36568-4744.circ.txt) | 44984 | 402 | 36568 | 31824 | 4744 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ.txt](./aes-decipher/aes192-decipher-a8352-ad80-g42616-gd378-xx34264-5448.circ.txt) | 192 | 8352 | 80 | 42616 | [<ins><strong>378</strong></ins>](./aes-decipher/aes192-decipher-a8352-ad80-g42616-gd378-xx34264-5448.circ.txt) | 34264 | 28816 | 5448 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-256 KeyExpansion: Flat circuits</h3></summary>

#### KeyExpansion using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| [circ.txt](./aes-keyexp/aes256-keyexp-a1508-ad65-g8892-gd495-xx7384-1047.circ.txt) | 256 | [<ins><strong>1508</strong></ins>](./aes-keyexp/aes256-keyexp-a1508-ad65-g8892-gd495-xx7384-1047.circ.txt) | 65 | 8892 | 495 | 7384 | 6337 | 1047 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | 1 |
| [circ.txt](./aes-keyexp/aes256-keyexp-a1768-ad52-g8320-gd247-xx6552-215.circ.txt) | 256 | 1768 | [<ins><strong>52</strong></ins>](./aes-keyexp/aes256-keyexp-a1768-ad52-g8320-gd247-xx6552-215.circ.txt) | 8320 | [<ins><strong>247</strong></ins>](./aes-keyexp/aes256-keyexp-a1768-ad52-g8320-gd247-xx6552-215.circ.txt) | 6552 | 6337 | 215 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 2,3 |
| [circ.txt](./aes-keyexp/aes256-keyexp-a1664-ad65-g7384-gd351-xx5720-163.circ.txt) | 256 | 1664 | 65 | [<ins><strong>7384</strong></ins>](./aes-keyexp/aes256-keyexp-a1664-ad65-g7384-gd351-xx5720-163.circ.txt) | 351 | 5720 | 5557 | 163 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### KeyExpansion using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| [circ.txt](./aes-keyexp/aes256-keyexp-a1872-ad52-g8840-gd234-xx6968-1151.circ.txt) | 256 | [<ins><strong>1872</strong></ins>](./aes-keyexp/aes256-keyexp-a1872-ad52-g8840-gd234-xx6968-1151.circ.txt) | 52 | [<ins><strong>8840</strong></ins>](./aes-keyexp/aes256-keyexp-a1872-ad52-g8840-gd234-xx6968-1151.circ.txt) | 234 | 6968 | 5817 | 1151 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 1,4 |
| [circ.txt](./aes-keyexp/aes256-keyexp-a2444-ad39-g13364-gd247-xx10920-215.circ.txt) | 256 | 2444 | [<ins><strong>39</strong></ins>](./aes-keyexp/aes256-keyexp-a2444-ad39-g13364-gd247-xx10920-215.circ.txt) | 13364 | 247 | 10920 | 10705 | 215 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 2 |
| [circ.txt](./aes-keyexp/aes256-keyexp-a2340-ad52-g9516-gd208-xx7176-1359.circ.txt) | 256 | 2340 | 52 | 9516 | [<ins><strong>208</strong></ins>](./aes-keyexp/aes256-keyexp-a2340-ad52-g9516-gd208-xx7176-1359.circ.txt) | 7176 | 5817 | 1359 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-256 Cipher: Flat circuits</h3></summary>

#### Cipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-cipher/aes256-cipher-a6496-ad70-g37632-gd544-xx31136-4480.circ.txt) | 256 | [<ins><strong>6496</strong></ins>](./aes-cipher/aes256-cipher-a6496-ad70-g37632-gd544-xx31136-4480.circ.txt) | 70 | 37632 | 544 | 31136 | 26656 | 4480 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ.txt](./aes-cipher/aes256-cipher-a7616-ad56-g35636-gd264-xx28020-896.circ.txt) | 256 | 7616 | [<ins><strong>56</strong></ins>](./aes-cipher/aes256-cipher-a7616-ad56-g35636-gd264-xx28020-896.circ.txt) | 35636 | [<ins><strong>264</strong></ins>](./aes-cipher/aes256-cipher-a7616-ad56-g35636-gd264-xx28020-896.circ.txt) | 28020 | 27124 | 896 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ.txt](./aes-cipher/aes256-cipher-a7168-ad70-g31136-gd402-xx23968-672.circ.txt) | 256 | 7168 | 70 | [<ins><strong>31136</strong></ins>](./aes-cipher/aes256-cipher-a7168-ad70-g31136-gd402-xx23968-672.circ.txt) | 402 | 23968 | 23296 | 672 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Cipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-cipher/aes256-cipher-a8064-ad56-g37408-gd276-xx29344-4928.circ.txt) | 256 | [<ins><strong>8064</strong></ins>](./aes-cipher/aes256-cipher-a8064-ad56-g37408-gd276-xx29344-4928.circ.txt) | 56 | [<ins><strong>37408</strong></ins>](./aes-cipher/aes256-cipher-a8064-ad56-g37408-gd276-xx29344-4928.circ.txt) | 276 | 29344 | 24416 | 4928 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1,4 |
| [circ.txt](./aes-cipher/aes256-cipher-a10528-ad42-g57364-gd264-xx46836-896.circ.txt) | 256 | 10528 | [<ins><strong>42</strong></ins>](./aes-cipher/aes256-cipher-a10528-ad42-g57364-gd264-xx46836-896.circ.txt) | 57364 | 264 | 46836 | 45940 | 896 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ.txt](./aes-cipher/aes256-cipher-a10080-ad56-g40788-gd222-xx30708-5824.circ.txt) | 256 | 10080 | 56 | 40788 | [<ins><strong>222</strong></ins>](./aes-cipher/aes256-cipher-a10080-ad56-g40788-gd222-xx30708-5824.circ.txt) | 30708 | 24884 | 5824 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-256 Encipher: Flat circuits</h3></summary>

#### Encipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-encipher/aes256-encipher-a8004-ad70-g46524-gd544-xx38520-5527.circ.txt) | 256 | [<ins><strong>8004</strong></ins>](./aes-encipher/aes256-encipher-a8004-ad70-g46524-gd544-xx38520-5527.circ.txt) | 70 | 46524 | 544 | 38520 | 32993 | 5527 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ.txt](./aes-encipher/aes256-encipher-a9384-ad56-g43956-gd264-xx34572-1111.circ.txt) | 256 | 9384 | [<ins><strong>56</strong></ins>](./aes-encipher/aes256-encipher-a9384-ad56-g43956-gd264-xx34572-1111.circ.txt) | 43956 | [<ins><strong>264</strong></ins>](./aes-encipher/aes256-encipher-a9384-ad56-g43956-gd264-xx34572-1111.circ.txt) | 34572 | 33461 | 1111 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ.txt](./aes-encipher/aes256-encipher-a8832-ad70-g38520-gd402-xx29688-835.circ.txt) | 256 | 8832 | 70 | [<ins><strong>38520</strong></ins>](./aes-encipher/aes256-encipher-a8832-ad70-g38520-gd402-xx29688-835.circ.txt) | 402 | 29688 | 28853 | 835 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Encipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-encipher/aes256-encipher-a9936-ad56-g46248-gd276-xx36312-6079.circ.txt) | 256 | [<ins><strong>9936</strong></ins>](./aes-encipher/aes256-encipher-a9936-ad56-g46248-gd276-xx36312-6079.circ.txt) | 56 | [<ins><strong>46248</strong></ins>](./aes-encipher/aes256-encipher-a9936-ad56-g46248-gd276-xx36312-6079.circ.txt) | 276 | 36312 | 30233 | 6079 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1,4 |
| [circ.txt](./aes-encipher/aes256-encipher-a12972-ad42-g70728-gd264-xx57756-1111.circ.txt) | 256 | 12972 | [<ins><strong>42</strong></ins>](./aes-encipher/aes256-encipher-a12972-ad42-g70728-gd264-xx57756-1111.circ.txt) | 70728 | 264 | 57756 | 56645 | 1111 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ.txt](./aes-encipher/aes256-encipher-a12420-ad56-g50304-gd222-xx37884-7183.circ.txt) | 256 | 12420 | 56 | 50304 | [<ins><strong>222</strong></ins>](./aes-encipher/aes256-encipher-a12420-ad56-g50304-gd222-xx37884-7183.circ.txt) | 37884 | 30701 | 7183 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [Link](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-256 InvCipher: Flat circuits</h3></summary>

#### InvCipher using InvSbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-invcipher/aes256-invcipher-a6496-ad70-g40328-gd567-xx33832-6496.circ.txt) | 256 | [<ins><strong>6496</strong></ins>](./aes-invcipher/aes256-invcipher-a6496-ad70-g40328-gd567-xx33832-6496.circ.txt) | 70 | 40328 | 567 | 33832 | 27336 | 6496 | [Link](./aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ.txt](./aes-invcipher/aes256-invcipher-a7616-ad56-g39528-gd290-xx31912-4032.circ.txt) | 256 | 7616 | [<ins><strong>56</strong></ins>](./aes-invcipher/aes256-invcipher-a7616-ad56-g39528-gd290-xx31912-4032.circ.txt) | 39528 | [<ins><strong>290</strong></ins>](./aes-invcipher/aes256-invcipher-a7616-ad56-g39528-gd290-xx31912-4032.circ.txt) | 31912 | 27880 | 4032 | [Link](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ.txt](./aes-invcipher/aes256-invcipher-a7168-ad70-g32936-gd497-xx25768-2016.circ.txt) | 256 | 7168 | 70 | [<ins><strong>32936</strong></ins>](./aes-invcipher/aes256-invcipher-a7168-ad70-g32936-gd497-xx25768-2016.circ.txt) | 497 | 25768 | 23752 | 2016 | [Link](./aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### InvCipher using InvSbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| [circ.txt](./aes-invcipher/aes256-invcipher-a8064-ad56-g40776-gd315-xx32712-5376.circ.txt) | 256 | [<ins><strong>8064</strong></ins>](./aes-invcipher/aes256-invcipher-a8064-ad56-g40776-gd315-xx32712-5376.circ.txt) | [<ins><strong>56</strong></ins>](./aes-invcipher/aes256-invcipher-a8064-ad56-g40776-gd315-xx32712-5376.circ.txt) | [<ins><strong>40776</strong></ins>](./aes-invcipher/aes256-invcipher-a8064-ad56-g40776-gd315-xx32712-5376.circ.txt) | 315 | 32712 | 27336 | 5376 | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ.txt](./aes-invcipher/aes256-invcipher-a8064-ad56-g42440-gd276-xx34376-5376.circ.txt) | 256 | [<ins><strong>8064</strong></ins>](./aes-invcipher/aes256-invcipher-a8064-ad56-g42440-gd276-xx34376-5376.circ.txt) | [<ins><strong>56</strong></ins>](./aes-invcipher/aes256-invcipher-a8064-ad56-g42440-gd276-xx34376-5376.circ.txt) | 42440 | [<ins><strong>276</strong></ins>](./aes-invcipher/aes256-invcipher-a8064-ad56-g42440-gd276-xx34376-5376.circ.txt) | 34376 | 29000 | 5376 | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-256 Decipher: Flat circuits</h3></summary>

#### Decipher using [Inv]Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|
| [circ.txt](./aes-decipher/aes256-decipher-a8004-ad135-g49220-gd1062-xx41216-7543.circ.txt) | 256 | [<ins><strong>8004</strong></ins>](./aes-decipher/aes256-decipher-a8004-ad135-g49220-gd1062-xx41216-7543.circ.txt) | 135 | 49220 | 1062 | 41216 | 33673 | 7543 | [Link](./aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ.txt](./aes-decipher/aes256-decipher-a9384-ad108-g47848-gd537-xx38464-4247.circ.txt) | 256 | 9384 | [<ins><strong>108</strong></ins>](./aes-decipher/aes256-decipher-a9384-ad108-g47848-gd537-xx38464-4247.circ.txt) | 47848 | [<ins><strong>537</strong></ins>](./aes-decipher/aes256-decipher-a9384-ad108-g47848-gd537-xx38464-4247.circ.txt) | 38464 | 34217 | 4247 | [Link](./aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ.txt](./aes-decipher/aes256-decipher-a8832-ad135-g40320-gd848-xx31488-2179.circ.txt) | 256 | 8832 | 135 | [<ins><strong>40320</strong></ins>](./aes-decipher/aes256-decipher-a8832-ad135-g40320-gd848-xx31488-2179.circ.txt) | 848 | 31488 | 29309 | 2179 | [Link](./aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Decipher using [Inv]Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|
| [circ.txt](./aes-decipher/aes256-decipher-a9936-ad108-g49616-gd549-xx39680-6527.circ.txt) | 256 | [<ins><strong>9936</strong></ins>](./aes-decipher/aes256-decipher-a9936-ad108-g49616-gd549-xx39680-6527.circ.txt) | 108 | [<ins><strong>49616</strong></ins>](./aes-decipher/aes256-decipher-a9936-ad108-g49616-gd549-xx39680-6527.circ.txt) | 549 | 39680 | 33153 | 6527 | [Link](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ.txt](./aes-decipher/aes256-decipher-a10508-ad95-g55804-gd523-xx45296-5591.circ.txt) | 256 | 10508 | [<ins><strong>95</strong></ins>](./aes-decipher/aes256-decipher-a10508-ad95-g55804-gd523-xx45296-5591.circ.txt) | 55804 | 523 | 45296 | 39705 | 5591 | [Link](./aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ.txt](./aes-decipher/aes256-decipher-a10404-ad108-g51956-gd484-xx41552-6735.circ.txt) | 256 | 10404 | 108 | 51956 | [<ins><strong>484</strong></ins>](./aes-decipher/aes256-decipher-a10404-ad108-g51956-gd484-xx41552-6735.circ.txt) | 41552 | 34817 | 6735 | [Link](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [Link](./aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [Link](./aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.
</details>
</details>
