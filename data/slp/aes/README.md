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

Example circuits for the AES S-box (8-bit to 8-bit function).

#### S-box with #AND ≤ 34

| File | #AND<br>(*A*) | AND<br>depth | #Gates<br>(*G*) | Gate<br>depth | XX<br>(*x*+*x'*) | #XOR<br>(*x*) | #XNOR<br>(*x'*) |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](./aes-sbox/aes-sbox-a28-ad4-g131-gd32-xx103-47.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 131 | 32 | 103 | 56 | 47 |
| [circ](./aes-sbox/aes-sbox-a28-ad4-g154-gd15-xx126-28.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 154 | 15 | 126 | 98 | 28 |
| [circ](./aes-sbox/aes-sbox-a28-ad4-g177-gd14-xx149-24.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 177 | <ins><strong>14</strong></ins> | 149 | 125 | 24 |
| [circ](./aes-sbox/aes-sbox-a28-ad5-g124-gd27-xx96-26.circ.txt) | <ins><strong>28</strong></ins> | 5 | 124 | 27 | 96 | 70 | 26 |
| [circ](./aes-sbox/aes-sbox-a28-ad5-g150-gd15-xx122-23.circ.txt) | <ins><strong>28</strong></ins> | 5 | 150 | 15 | 122 | 99 | 23 |
| [circ](./aes-sbox/aes-sbox-a32-ad5-g110-gd22-xx78-3.circ.txt) | 32 | 5 | <ins><strong>110</strong></ins> | 22 | 78 | 75 | 3 |
| [circ](./aes-sbox/aes-sbox-a34-ad4-g110-gd22-xx76-3.circ.txt) | 34 | <ins><strong>4</strong></ins> | <ins><strong>110</strong></ins> | 22 | <ins><strong>76</strong></ins> | 73 | 3 |
| [circ](./aes-sbox/aes-sbox-a34-ad4-g125-gd15-xx91-4.circ.txt) | 34 | <ins><strong>4</strong></ins> | 125 | 15 | 91 | 87 | 4 |

Note: The A28 sbox circuits listed in the table above are either (AD4/\{G154/GD15, G177/GD14\}, AD5\{G124/GD27,G150/GD15\}) contributions communicated by Milad Nasr (@ Anthropic) on 2026-09-24, or were derived therefrom via linear optimization.

#### S-box with #AND > 34

| File | #AND<br>(*A*) | AND<br>depth | #Gates<br>(*G*) | Gate<br>depth | XX<br>(*x*+*x'*) | #XOR<br>(*x*) | #XNOR<br>(*x'*) |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](./aes-sbox/aes-sbox-a35-ad3-g139-gd27-xx104-38.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 139 | 27 | 104 | 66 | 38 |
| [circ](./aes-sbox/aes-sbox-a35-ad3-g157-gd15-xx122-31.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 157 | 15 | 122 | 91 | 31 |
| [circ](./aes-sbox/aes-sbox-a35-ad3-g183-gd13-xx148-29.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 183 | 13 | 148 | 119 | 29 |
| [circ](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 36 | 4 | <ins><strong>138</strong></ins> | 14 | <ins><strong>102</strong></ins> | 80 | 22 |
| [circ](./aes-sbox/aes-sbox-a38-ad4-g148-gd13-xx110-24.circ.txt) | 38 | 4 | 148 | 13 | 110 | 86 | 24 |
| [circ](./aes-sbox/aes-sbox-a39-ad3-g233-gd12-xx194-30.circ.txt) | 39 | <ins><strong>3</strong></ins> | 233 | 12 | 194 | 164 | 30 |
| [circ](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 45 | 4 | 151 | 12 | 106 | 80 | 26 |
| [circ](./aes-sbox/aes-sbox-a49-ad4-g164-gd11-xx115-32.circ.txt) | 49 | 4 | 164 | <ins><strong>11</strong></ins> | 115 | 83 | 32 |
| [circ](./aes-sbox/aes-sbox-a51-ad4-g157-gd11-xx106-28.circ.txt) | 51 | 4 | 157 | <ins><strong>11</strong></ins> | 106 | 78 | 28 |

Note: The A35/AD3 sbox circuits listed in the table above are either (G157/GD15, G183/GD13) contributions communicated by Milad Nasr (@ Anthropic) on 2026-09-24, or were derived therefrom via linear optimization.

</details>
<details open>
<summary><h3>AES Inverse S-box</h3></summary>

Example circuits for the AES Inverse S-box (8-bit to 8-bit function).

#### Inverse S-box with #AND ≤ 34

| File | #AND<br>(*A*) | AND<br>depth | #Gates<br>(*G*) | Gate<br>depth | XX<br>(*x*+*x'*) | #XOR<br>(*x*) | #XNOR<br>(*x'*) |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](./aes-invsbox/aes-invsbox-a28-ad4-g141-gd31-xx113-28.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 141 | 31 | 113 | 85 | 28 |
| [circ](./aes-invsbox/aes-invsbox-a28-ad4-g174-gd15-xx146-33.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 174 | 15 | 146 | 113 | 33 |
| [circ](./aes-invsbox/aes-invsbox-a28-ad4-g214-gd14-xx186-32.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 214 | <ins><strong>14</strong></ins> | 186 | 154 | 32 |
| [circ](./aes-invsbox/aes-invsbox-a28-ad5-g126-gd27-xx98-26.circ.txt) | <ins><strong>28</strong></ins> | 5 | 126 | 27 | 98 | 72 | 26 |
| [circ](./aes-invsbox/aes-invsbox-a28-ad5-g152-gd15-xx124-18.circ.txt) | <ins><strong>28</strong></ins> | 5 | 152 | 15 | 124 | 106 | 18 |
| [circ](./aes-invsbox/aes-invsbox-a32-ad5-g112-gd26-xx80-9.circ.txt) | 32 | 5 | <ins><strong>112</strong></ins> | 26 | <ins><strong>80</strong></ins> | 71 | 9 |
| [circ](./aes-invsbox/aes-invsbox-a34-ad4-g114-gd26-xx80-11.circ.txt) | 34 | <ins><strong>4</strong></ins> | 114 | 26 | <ins><strong>80</strong></ins> | 69 | 11 |
| [circ](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | 34 | <ins><strong>4</strong></ins> | 134 | 15 | 100 | 82 | 18 |
| [circ](./aes-invsbox/aes-invsbox-a34-ad4-g194-gd14-xx160-20.circ.txt) | 34 | <ins><strong>4</strong></ins> | 194 | <ins><strong>14</strong></ins> | 160 | 140 | 20 |

Note: The A28 invsbox circuits listed in the table above are either (AD5/\{G126/GD27,G152/GD15\}) contributions communicated on 2026-09-24 by Milad Nasr (@ Anthropic), or were derived later on via linear optimization of contributed A28 [sbox circuits](./aes-sbox/).

#### Inverse S-box with #AND > 34

| File | #AND<br>(*A*) | AND<br>depth | #Gates<br>(*G*) | Gate<br>depth | XX<br>(*x*+*x'*) | #XOR<br>(*x*) | #XNOR<br>(*x'*) |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](./aes-invsbox/aes-invsbox-a35-ad3-g158-gd28-xx123-34.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 158 | 28 | 123 | 89 | 34 |
| [circ](./aes-invsbox/aes-invsbox-a35-ad3-g185-gd15-xx150-26.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 185 | 15 | 150 | 124 | 26 |
| [circ](./aes-invsbox/aes-invsbox-a35-ad3-g203-gd14-xx168-32.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 203 | 14 | 168 | 136 | 32 |
| [circ](./aes-invsbox/aes-invsbox-a36-ad4-g140-gd14-xx104-24.circ.txt) | 36 | 4 | <ins><strong>140</strong></ins> | 14 | <ins><strong>104</strong></ins> | 80 | 24 |
| [circ](./aes-invsbox/aes-invsbox-a37-ad4-g172-gd13-xx135-24.circ.txt) | 37 | 4 | 172 | <ins><strong>13</strong></ins> | 135 | 111 | 24 |
| [circ](./aes-invsbox/aes-invsbox-a39-ad4-g171-gd13-xx132-22.circ.txt) | 39 | 4 | 171 | <ins><strong>13</strong></ins> | 132 | 110 | 22 |

Note: The A35/AD3 invsbox circuits listed in the table above were derived by linear transformation and optimization of A35/AD3 [sbox circuits](./aes-sbox/) (G139/GD31, G157/GD15, G183/GD13) externally contributed by Milad Nasr (@ Anthropic) on 2026-09-24.

</details>
<details open>
<summary><h3>[Inv]S-Box circuits added around 2020</h3></summary>

The following historical table was retrieved/adapted from an old version of the [NIST Circuit Complexity list of circuits](https://csrc.nist.gov/Projects/circuit-complexity/list-of-circuits).

| File | Direction | #AND<br>(*A*) | AND<br>depth | #Gates<br>(*G*) | Gate<br>depth | XX<br>(*x*+*x'*) | #XOR<br>(*x*) | #XNOR<br>(*x'*) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| [slp](./old-2020/aes-sbox-fwd-g115-a32-d28-ad6.slp) | Forward | <ins><strong>32</strong></ins> | 6 | 115 | 28 | 83 | 79 | 4 |
| [slp](./old-2020/aes-sbox-fwd-g113-a32-d27-ad6.slp) | Forward | <ins><strong>32</strong></ins> | 6 | <ins><strong>113</strong></ins> | 27 | 81 | 77 | 4 |
| [slp](./old-2020/aes-sbox-fwd-g128-a34-d16-ad4.slp) | Forward | 34 | <ins><strong>4</strong></ins> | 128 | <ins><strong>16</strong></ins> | 94 | 90 | 4 |
| [slp](./old-2020/aes-sbox-rev-g121-a34-d21-ad4.slp) | Inverse | 34 | 4 | <ins><strong>121</strong></ins> | 21 | 87 | 83 | 4 |
| [slp](./old-2020/aes-sbox-rev-g127-a34-d16-ad4.slp) | Inverse | 34 | 4 | 127 | <ins><strong>16</strong></ins> | 93 | 83 | 10 |

</details>
<details open>
<summary><h3>AES MixColumns</h3></summary>

Example circuits for AES MixColumns (32-bit to 32-bit linear function).

| File | Depth | #XOR |
|---|---:|---:|
| [circ](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | <ins><strong>3</strong></ins> | 97 |
| [circ](./aes-mixcols/aes-mixcols-xor90-depth4.circ.txt) | 4 | 90 |
| [circ](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 5 | <ins><strong>88</strong></ins> |

</details>
<details open>
<summary><h3>AES InvMixColumns</h3></summary>

Example circuits for AES Inverse MixColumns (32-bit to 32-bit linear function).

| File | Depth | #XOR |
|---|---:|---:|
| [circ](./aes-invmixcols/aes-invmixcols-xor127-depth5.circ.txt) | <ins><strong>5</strong></ins> | 127 |
| [circ](./aes-invmixcols/aes-invmixcols-xor107-depth6.circ.txt) | 6 | 107 |
| [circ](./aes-invmixcols/aes-invmixcols-xor102-depth7.circ.txt) | 7 | 102 |
| [circ](./aes-invmixcols/aes-invmixcols-xor96-depth8.circ.txt) | 8 | 96 |
| [circ](./aes-invmixcols/aes-invmixcols-xor95-depth9.circ.txt) | 9 | 95 |
| [circ](./aes-invmixcols/aes-invmixcols-xor93-depth12.circ.txt) | 12 | 93 |
| [circ](./aes-invmixcols/aes-invmixcols-xor92-depth13.circ.txt) | 13 | 92 |
| [circ](./aes-invmixcols/aes-invmixcols-xor91-depth15.circ.txt) | 15 | <ins><strong>91</strong></ins> |


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

`VXOR` and `VXNOR` are succinct notation for parallel execution of `XOR` and `XNOR` gates, respectively.

| Operation | File | Key<br>size | #sbox | #inv<br>sbox | #mixcols | #inv<br>mixcols | #VXOR<br>(#XOR) | #VXNOR<br>(#XNOR) | #XNOR |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| KeyExpansion | [aes128-keyexp](./aes-fold1/aes128-keyexp-fold1.circ.txt) | 128 | 40 | — | — | — | 40 (1264) | 2 (8) | 8 |
| KeyExpansion | [aes192-keyexp](./aes-fold1/aes192-keyexp-fold1.circ.txt) | 192 | 32 | — | — | — | 46 (1464) | — | 8 |
| KeyExpansion | [aes256-keyexp](./aes-fold1/aes256-keyexp-fold1.circ.txt) | 256 | 52 | — | — | — | 52 (1657) | — | 7 |
| Cipher | [aes128-cipher](./aes-fold1/aes128-cipher-fold1.circ.txt) | 128 | 160 | — | 36 | — | 11 (1408) | — | — |
| Cipher | [aes192-cipher](./aes-fold1/aes192-cipher-fold1.circ.txt) | 192 | 192 | — | 44 | — | 14 (1664) | — | — |
| Cipher | [aes256-cipher](./aes-fold1/aes256-cipher-fold1.circ.txt) | 256 | 224 | — | 52 | — | 15 (1920) | — | — |
| Encipher | [aes128-encipher](./aes-fold1/aes128-encipher-fold1.circ.txt) | 128 | 200 | — | 36 | — | 51 (2672) | 2 (8) | 8 |
| Encipher | [aes192-encipher](./aes-fold1/aes192-encipher-fold1.circ.txt) | 192 | 224 | — | 44 | — | 60 (3128) | — | 8 |
| Encipher | [aes256-encipher](./aes-fold1/aes256-encipher-fold1.circ.txt) | 256 | 276 | — | 52 | — | 67 (3577) | — | 7 |
| InvCipher | [aes128-invcipher](./aes-fold1/aes128-invcipher-fold1.circ.txt) | 128 | — | 160 | — | 36 | 11 (1408) | — | — |
| InvCipher | [aes192-invcipher](./aes-fold1/aes192-invcipher-fold1.circ.txt) | 192 | — | 192 | — | 44 | 14 (1664) | — | — |
| InvCipher | [aes256-invcipher](./aes-fold1/aes256-invcipher-fold1.circ.txt) | 256 | — | 224 | — | 52 | 15 (1920) | — | — |
| Decipher | [aes128-decipher](./aes-fold1/aes128-decipher-fold1.circ.txt) | 128 | 40 | 160 | — | 36 | 51 (2672) | 2 (8) | 8 |
| Decipher | [aes192-decipher](./aes-fold1/aes192-decipher-fold1.circ.txt) | 192 | 32 | 192 | — | 44 | 60 (3128) | — | 8 |
| Decipher | [aes256-decipher](./aes-fold1/aes256-decipher-fold1.circ.txt) | 256 | 52 | 224 | — | 52 | 67 (3577) | — | 7 |


</details>
<details open>
<summary><h2>AES flat circuits</h2></summary>

This section contains "flat circuits", built by flattening the intermediate components (\[inv\]sbox and \[inv\]mixcols, VX\[N\]OR) of folded circuits, so that the entire circuit uses only basic Boolean gates (AND, XOR, XNOR).

The exemplified compilations may be using components that are no longer optimal, since new optimized components may be posted more frequently than flat circuits.

There are up to 4! = 24 tuple metrics corresponding to the possible orderings of (A, AD, G, GD). For succinctness, the tables consider only four tuple metrics (TM): (1) A-AD-G-GD, (2) AD-GD-G-A, (3) GD-G-AD-A, and (4) G-A-GD-AD.

Note: In selected columns, values <ins><strong>underlined in bold</ins></strong> indicate the lowest displayed value.


<details open>
<summary><h3>AES-128 KeyExpansion: Flat circuits</h3></summary>

#### KeyExpansion using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| [circ](./aes-keyexp/aes128-keyexp-a1160-ad50-g6840-gd381-xx5680-816.circ.txt) | 128 | <ins><strong>1160</strong></ins> | 50 | 6840 | 381 | 5680 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | 1 |
| [circ](./aes-keyexp/aes128-keyexp-a1360-ad40-g6400-gd190-xx5040-176.circ.txt) | 128 | 1360 | <ins><strong>40</strong></ins> | 6400 | <ins><strong>190</strong></ins> | 5040 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 2,3 |
| [circ](./aes-keyexp/aes128-keyexp-a1280-ad50-g5680-gd270-xx4400-136.circ.txt) | 128 | 1280 | 50 | <ins><strong>5680</strong></ins> | 270 | 4400 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### KeyExpansion using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| [circ](./aes-keyexp/aes128-keyexp-a1440-ad40-g6800-gd180-xx5360-896.circ.txt) | 128 | <ins><strong>1440</strong></ins> | 40 | <ins><strong>6800</strong></ins> | 180 | 5360 | [A36/AD4/G138/GD14](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 1,4 |
| [circ](./aes-keyexp/aes128-keyexp-a1880-ad30-g10280-gd190-xx8400-176.circ.txt) | 128 | 1880 | <ins><strong>30</strong></ins> | 10280 | 190 | 8400 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 2 |
| [circ](./aes-keyexp/aes128-keyexp-a1800-ad40-g7320-gd160-xx5520-1056.circ.txt) | 128 | 1800 | 40 | 7320 | <ins><strong>160</strong></ins> | 5520 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-128 Cipher: Flat circuits</h3></summary>

#### Cipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-cipher/aes128-cipher-a4640-ad50-g26816-gd388-xx22176-3200.circ.txt) | 128 | <ins><strong>4640</strong></ins> | 50 | 26816 | 388 | 22176 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](./aes-cipher/aes128-cipher-a5440-ad40-g25380-gd188-xx19940-640.circ.txt) | 128 | 5440 | <ins><strong>40</strong></ins> | 25380 | <ins><strong>188</strong></ins> | 19940 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](./aes-cipher/aes128-cipher-a5120-ad50-g22176-gd286-xx17056-480.circ.txt) | 128 | 5120 | 50 | <ins><strong>22176</strong></ins> | 286 | 17056 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Cipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-cipher/aes128-cipher-a7520-ad30-g40900-gd188-xx33380-640.circ.txt) | 128 | 7520 | <ins><strong>30</strong></ins> | 40900 | 188 | 33380 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](./aes-cipher/aes128-cipher-a7200-ad40-g29060-gd158-xx21860-4160.circ.txt) | 128 | <ins><strong>7200</strong></ins> | 40 | <ins><strong>29060</strong></ins> | <ins><strong>158</strong></ins> | 21860 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-128 Encipher: Flat circuits</h3></summary>

#### Encipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-encipher/aes128-encipher-a5800-ad50-g33656-gd388-xx27856-4016.circ.txt) | 128 | <ins><strong>5800</strong></ins> | 50 | 33656 | 388 | 27856 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](./aes-encipher/aes128-encipher-a6800-ad40-g31780-gd191-xx24980-816.circ.txt) | 128 | 6800 | <ins><strong>40</strong></ins> | 31780 | <ins><strong>191</strong></ins> | 24980 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](./aes-encipher/aes128-encipher-a6400-ad50-g27856-gd286-xx21456-616.circ.txt) | 128 | 6400 | 50 | <ins><strong>27856</strong></ins> | 286 | 21456 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Encipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-encipher/aes128-encipher-a9400-ad30-g51180-gd191-xx41780-816.circ.txt) | 128 | 9400 | <ins><strong>30</strong></ins> | 51180 | 191 | 41780 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](./aes-encipher/aes128-encipher-a9000-ad40-g36380-gd161-xx27380-5216.circ.txt) | 128 | <ins><strong>9000</strong></ins> | 40 | <ins><strong>36380</strong></ins> | <ins><strong>161</strong></ins> | 27380 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-128 InvCipher: Flat circuits</h3></summary>

#### InvCipher using InvSbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-invcipher/aes128-invcipher-a4640-ad50-g28712-gd403-xx24072-4640.circ.txt) | 128 | <ins><strong>4640</strong></ins> | 50 | 28712 | 403 | 24072 | [A29/AD5/G145/GD32](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A29-AD5-ADP-9-1-2-2-15/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](./aes-invcipher/aes128-invcipher-a5440-ad40-g28104-gd206-xx22664-2880.circ.txt) | 128 | 5440 | <ins><strong>40</strong></ins> | 28104 | <ins><strong>206</strong></ins> | 22664 | [A34/AD4/G134/GD15](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](./aes-invcipher/aes128-invcipher-a5120-ad50-g23432-gd353-xx18312-1440.circ.txt) | 128 | 5120 | 50 | <ins><strong>23432</strong></ins> | 353 | 18312 | [A32/AD5/G112/GD27](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A32-AD5-ADP-9-1-2-8-12/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### InvCipher using InvSbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-invcipher/aes128-invcipher-a5760-ad40-g30184-gd196-xx24424-3840.circ.txt) | 128 | 5760 | 40 | 30184 | 196 | 24424 | [A36/AD4/G147/GD14](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A-greater-than-34/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 1,2,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-128 Decipher: Flat circuits</h3></summary>

#### Decipher using [Inv]Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: |
| [circ](./aes-decipher/aes128-decipher-a5800-ad100-g35552-gd784-xx29752-5456.circ.txt) | 128 | <ins><strong>5800</strong></ins> | 100 | 35552 | 784 | 29752 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [A29/AD5/G145/GD32](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A29-AD5-ADP-9-1-2-2-15/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](./aes-decipher/aes128-decipher-a6800-ad80-g34504-gd396-xx27704-3056.circ.txt) | 128 | 6800 | <ins><strong>80</strong></ins> | 34504 | <ins><strong>396</strong></ins> | 27704 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [A34/AD4/G134/GD15](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](./aes-decipher/aes128-decipher-a6400-ad100-g29112-gd623-xx22712-1576.circ.txt) | 128 | 6400 | 100 | <ins><strong>29112</strong></ins> | 623 | 22712 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [A32/AD5/G112/GD27](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A32-AD5-ADP-9-1-2-8-12/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Decipher using [Inv]Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: |
| [circ](./aes-decipher/aes128-decipher-a7640-ad70-g40464-gd386-xx32824-4016.circ.txt) | 128 | 7640 | <ins><strong>70</strong></ins> | 40464 | 386 | 32824 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [A36/AD4/G147/GD14](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A-greater-than-34/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ](./aes-decipher/aes128-decipher-a7560-ad80-g37504-gd356-xx29944-4896.circ.txt) | 128 | <ins><strong>7560</strong></ins> | 80 | <ins><strong>37504</strong></ins> | <ins><strong>356</strong></ins> | 29944 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [A36/AD4/G147/GD14](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A-greater-than-34/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 1,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-192 KeyExpansion: Flat circuits</h3></summary>

#### KeyExpansion using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| [circ](./aes-keyexp/aes192-keyexp-a928-ad40-g5920-gd319-xx4992-648.circ.txt) | 192 | <ins><strong>928</strong></ins> | 40 | 5920 | 319 | 4992 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | 1 |
| [circ](./aes-keyexp/aes192-keyexp-a1088-ad32-g5568-gd166-xx4480-136.circ.txt) | 192 | 1088 | <ins><strong>32</strong></ins> | 5568 | <ins><strong>166</strong></ins> | 4480 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 2,3 |
| [circ](./aes-keyexp/aes192-keyexp-a1024-ad40-g4992-gd230-xx3968-104.circ.txt) | 192 | 1024 | 40 | <ins><strong>4992</strong></ins> | 230 | 3968 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### KeyExpansion using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| [circ](./aes-keyexp/aes192-keyexp-a1152-ad32-g5888-gd158-xx4736-712.circ.txt) | 192 | <ins><strong>1152</strong></ins> | 32 | <ins><strong>5888</strong></ins> | 158 | 4736 | [A36/AD4/G138/GD14](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 1,4 |
| [circ](./aes-keyexp/aes192-keyexp-a1504-ad24-g8672-gd166-xx7168-136.circ.txt) | 192 | 1504 | <ins><strong>24</strong></ins> | 8672 | 166 | 7168 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 2 |
| [circ](./aes-keyexp/aes192-keyexp-a1440-ad32-g6304-gd142-xx4864-840.circ.txt) | 192 | 1440 | 32 | 6304 | <ins><strong>142</strong></ins> | 4864 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-192 Cipher: Flat circuits</h3></summary>

#### Cipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-cipher/aes192-cipher-a5568-ad60-g32224-gd466-xx26656-3840.circ.txt) | 192 | <ins><strong>5568</strong></ins> | 60 | 32224 | 466 | 26656 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](./aes-cipher/aes192-cipher-a6528-ad48-g30508-gd226-xx23980-768.circ.txt) | 192 | 6528 | <ins><strong>48</strong></ins> | 30508 | <ins><strong>226</strong></ins> | 23980 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](./aes-cipher/aes192-cipher-a6144-ad60-g26656-gd344-xx20512-576.circ.txt) | 192 | 6144 | 60 | <ins><strong>26656</strong></ins> | 344 | 20512 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Cipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-cipher/aes192-cipher-a9024-ad36-g49132-gd226-xx40108-768.circ.txt) | 192 | 9024 | <ins><strong>36</strong></ins> | 49132 | 226 | 40108 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](./aes-cipher/aes192-cipher-a8640-ad48-g34924-gd190-xx26284-4992.circ.txt) | 192 | <ins><strong>8640</strong></ins> | 48 | <ins><strong>34924</strong></ins> | <ins><strong>190</strong></ins> | 26284 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-192 Encipher: Flat circuits</h3></summary>

#### Encipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-encipher/aes192-encipher-a6496-ad60-g38144-gd466-xx31648-4488.circ.txt) | 192 | <ins><strong>6496</strong></ins> | 60 | 38144 | 466 | 31648 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](./aes-encipher/aes192-encipher-a7616-ad48-g36076-gd226-xx28460-904.circ.txt) | 192 | 7616 | <ins><strong>48</strong></ins> | 36076 | <ins><strong>226</strong></ins> | 28460 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](./aes-encipher/aes192-encipher-a7168-ad60-g31648-gd344-xx24480-680.circ.txt) | 192 | 7168 | 60 | <ins><strong>31648</strong></ins> | 344 | 24480 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Encipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-encipher/aes192-encipher-a10528-ad36-g57804-gd226-xx47276-904.circ.txt) | 192 | 10528 | <ins><strong>36</strong></ins> | 57804 | 226 | 47276 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](./aes-encipher/aes192-encipher-a10080-ad48-g41228-gd190-xx31148-5832.circ.txt) | 192 | <ins><strong>10080</strong></ins> | 48 | <ins><strong>41228</strong></ins> | <ins><strong>190</strong></ins> | 31148 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-192 InvCipher: Flat circuits</h3></summary>

#### InvCipher using InvSbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-invcipher/aes192-invcipher-a5568-ad60-g34520-gd485-xx28952-5568.circ.txt) | 192 | <ins><strong>5568</strong></ins> | 60 | 34520 | 485 | 28952 | [A29/AD5/G145/GD32](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A29-AD5-ADP-9-1-2-2-15/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](./aes-invcipher/aes192-invcipher-a6528-ad48-g33816-gd248-xx27288-3456.circ.txt) | 192 | 6528 | <ins><strong>48</strong></ins> | 33816 | <ins><strong>248</strong></ins> | 27288 | [A34/AD4/G134/GD15](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](./aes-invcipher/aes192-invcipher-a6144-ad60-g28184-gd425-xx22040-1728.circ.txt) | 192 | 6144 | 60 | <ins><strong>28184</strong></ins> | 425 | 22040 | [A32/AD5/G112/GD27](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A32-AD5-ADP-9-1-2-8-12/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### InvCipher using InvSbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-invcipher/aes192-invcipher-a6912-ad48-g36312-gd236-xx29400-4608.circ.txt) | 192 | 6912 | 48 | 36312 | 236 | 29400 | [A36/AD4/G147/GD14](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A-greater-than-34/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 1,2,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-192 Decipher: Flat circuits</h3></summary>

#### Decipher using [Inv]Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: |
| [circ](./aes-decipher/aes192-decipher-a6496-ad100-g40440-gd804-xx33944-6216.circ.txt) | 192 | <ins><strong>6496</strong></ins> | 100 | 40440 | 804 | 33944 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [A29/AD5/G145/GD32](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A29-AD5-ADP-9-1-2-2-15/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](./aes-decipher/aes192-decipher-a7616-ad80-g39384-gd414-xx31768-3592.circ.txt) | 192 | 7616 | <ins><strong>80</strong></ins> | 39384 | <ins><strong>414</strong></ins> | 31768 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [A34/AD4/G134/GD15](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](./aes-decipher/aes192-decipher-a7168-ad100-g33176-gd655-xx26008-1832.circ.txt) | 192 | 7168 | 100 | <ins><strong>33176</strong></ins> | 655 | 26008 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [A32/AD5/G112/GD27](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A32-AD5-ADP-9-1-2-8-12/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Decipher using [Inv]Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: |
| [circ](./aes-decipher/aes192-decipher-a8416-ad72-g44984-gd402-xx36568-4744.circ.txt) | 192 | 8416 | <ins><strong>72</strong></ins> | 44984 | 402 | 36568 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [A36/AD4/G147/GD14](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A-greater-than-34/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ](./aes-decipher/aes192-decipher-a8352-ad80-g42616-gd378-xx34264-5448.circ.txt) | 192 | <ins><strong>8352</strong></ins> | 80 | <ins><strong>42616</strong></ins> | <ins><strong>378</strong></ins> | 34264 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [A36/AD4/G147/GD14](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A-greater-than-34/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 1,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-256 KeyExpansion: Flat circuits</h3></summary>

#### KeyExpansion using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| [circ](./aes-keyexp/aes256-keyexp-a1508-ad65-g8892-gd495-xx7384-1047.circ.txt) | 256 | <ins><strong>1508</strong></ins> | 65 | 8892 | 495 | 7384 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | 1 |
| [circ](./aes-keyexp/aes256-keyexp-a1768-ad52-g8320-gd247-xx6552-215.circ.txt) | 256 | 1768 | <ins><strong>52</strong></ins> | 8320 | <ins><strong>247</strong></ins> | 6552 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 2,3 |
| [circ](./aes-keyexp/aes256-keyexp-a1664-ad65-g7384-gd351-xx5720-163.circ.txt) | 256 | 1664 | 65 | <ins><strong>7384</strong></ins> | 351 | 5720 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### KeyExpansion using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| [circ](./aes-keyexp/aes256-keyexp-a1872-ad52-g8840-gd234-xx6968-1151.circ.txt) | 256 | <ins><strong>1872</strong></ins> | 52 | <ins><strong>8840</strong></ins> | 234 | 6968 | [A36/AD4/G138/GD14](./aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 1,4 |
| [circ](./aes-keyexp/aes256-keyexp-a2444-ad39-g13364-gd247-xx10920-215.circ.txt) | 256 | 2444 | <ins><strong>39</strong></ins> | 13364 | 247 | 10920 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 2 |
| [circ](./aes-keyexp/aes256-keyexp-a2340-ad52-g9516-gd208-xx7176-1359.circ.txt) | 256 | 2340 | 52 | 9516 | <ins><strong>208</strong></ins> | 7176 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 3 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-256 Cipher: Flat circuits</h3></summary>

#### Cipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-cipher/aes256-cipher-a6496-ad70-g37632-gd544-xx31136-4480.circ.txt) | 256 | <ins><strong>6496</strong></ins> | 70 | 37632 | 544 | 31136 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](./aes-cipher/aes256-cipher-a7616-ad56-g35636-gd264-xx28020-896.circ.txt) | 256 | 7616 | <ins><strong>56</strong></ins> | 35636 | <ins><strong>264</strong></ins> | 28020 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](./aes-cipher/aes256-cipher-a7168-ad70-g31136-gd402-xx23968-672.circ.txt) | 256 | 7168 | 70 | <ins><strong>31136</strong></ins> | 402 | 23968 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Cipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-cipher/aes256-cipher-a10528-ad42-g57364-gd264-xx46836-896.circ.txt) | 256 | 10528 | <ins><strong>42</strong></ins> | 57364 | 264 | 46836 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](./aes-cipher/aes256-cipher-a10080-ad56-g40788-gd222-xx30708-5824.circ.txt) | 256 | <ins><strong>10080</strong></ins> | 56 | <ins><strong>40788</strong></ins> | <ins><strong>222</strong></ins> | 30708 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-256 Encipher: Flat circuits</h3></summary>

#### Encipher using Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-encipher/aes256-encipher-a8004-ad70-g46524-gd544-xx38520-5527.circ.txt) | 256 | <ins><strong>8004</strong></ins> | 70 | 46524 | 544 | 38520 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](./aes-encipher/aes256-encipher-a9384-ad56-g43956-gd264-xx34572-1111.circ.txt) | 256 | 9384 | <ins><strong>56</strong></ins> | 43956 | <ins><strong>264</strong></ins> | 34572 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](./aes-encipher/aes256-encipher-a8832-ad70-g38520-gd402-xx29688-835.circ.txt) | 256 | 8832 | 70 | <ins><strong>38520</strong></ins> | 402 | 29688 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](./aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Encipher using Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-encipher/aes256-encipher-a12972-ad42-g70728-gd264-xx57756-1111.circ.txt) | 256 | 12972 | <ins><strong>42</strong></ins> | 70728 | 264 | 57756 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](./aes-encipher/aes256-encipher-a12420-ad56-g50304-gd222-xx37884-7183.circ.txt) | 256 | <ins><strong>12420</strong></ins> | 56 | <ins><strong>50304</strong></ins> | <ins><strong>222</strong></ins> | 37884 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](./aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-256 InvCipher: Flat circuits</h3></summary>

#### InvCipher using InvSbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-invcipher/aes256-invcipher-a6496-ad70-g40328-gd567-xx33832-6496.circ.txt) | 256 | <ins><strong>6496</strong></ins> | 70 | 40328 | 567 | 33832 | [A29/AD5/G145/GD32](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A29-AD5-ADP-9-1-2-2-15/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](./aes-invcipher/aes256-invcipher-a7616-ad56-g39528-gd290-xx31912-4032.circ.txt) | 256 | 7616 | <ins><strong>56</strong></ins> | 39528 | <ins><strong>290</strong></ins> | 31912 | [A34/AD4/G134/GD15](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](./aes-invcipher/aes256-invcipher-a7168-ad70-g32936-gd497-xx25768-2016.circ.txt) | 256 | 7168 | 70 | <ins><strong>32936</strong></ins> | 497 | 25768 | [A32/AD5/G112/GD27](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A32-AD5-ADP-9-1-2-8-12/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### InvCipher using InvSbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| [circ](./aes-invcipher/aes256-invcipher-a8064-ad56-g42440-gd276-xx34376-5376.circ.txt) | 256 | 8064 | 56 | 42440 | 276 | 34376 | [A36/AD4/G147/GD14](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A-greater-than-34/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 1,2,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

</details>
<details open>
<summary><h3>AES-256 Decipher: Flat circuits</h3></summary>

#### Decipher using [Inv]Sbox with A<=34

| File | \|k\| | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: |
| [circ](./aes-decipher/aes256-decipher-a8004-ad135-g49220-gd1062-xx41216-7543.circ.txt) | 256 | <ins><strong>8004</strong></ins> | 135 | 49220 | 1062 | 41216 | [A29/AD5/G139/GD35](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A29-AD5-ADP-9-1-2-2-15/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [A29/AD5/G145/GD32](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A29-AD5-ADP-9-1-2-2-15/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](./aes-decipher/aes256-decipher-a9384-ad108-g47848-gd537-xx38464-4247.circ.txt) | 256 | 9384 | <ins><strong>108</strong></ins> | 47848 | <ins><strong>537</strong></ins> | 38464 | [A34/AD4/G128/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A34-AD4-ADP-9-3-4-18/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [A34/AD4/G134/GD15](./aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](./aes-decipher/aes256-decipher-a8832-ad135-g40320-gd848-xx31488-2179.circ.txt) | 256 | 8832 | 135 | <ins><strong>40320</strong></ins> | 848 | 31488 | [A32/AD5/G110/GD23](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A32-AD5-ADP-9-1-2-8-12/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [A32/AD5/G112/GD27](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A32-AD5-ADP-9-1-2-8-12/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [G114/GD8](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.

#### Decipher using [Inv]Sbox with A>34

| File | \|k\| | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: |
| [circ](./aes-decipher/aes256-decipher-a10508-ad95-g55804-gd523-xx45296-5591.circ.txt) | 256 | 10508 | <ins><strong>95</strong></ins> | 55804 | 523 | 45296 | [A47/AD3/G225/GD15](../02-circs-discovered/aes-sbox-with-AND-XOR/sbox-fwd-A-greater-than-34/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [A36/AD4/G147/GD14](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A-greater-than-34/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ](./aes-decipher/aes256-decipher-a10404-ad108-g51956-gd484-xx41552-6735.circ.txt) | 256 | <ins><strong>10404</strong></ins> | 108 | <ins><strong>51956</strong></ins> | <ins><strong>484</strong></ins> | 41552 | [A45/AD4/G151/GD12](./aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [A36/AD4/G147/GD14](../02-circs-discovered/aes-invsbox-with-AND-XOR/invsbox-A-greater-than-34/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../02-circs-discovered/aes-invmixcols-with-AND-XOR/aes-invmixcols-xor146-depth5.circ.txt) | 1,3,4 |

\|k\| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.
</details>
</details>
