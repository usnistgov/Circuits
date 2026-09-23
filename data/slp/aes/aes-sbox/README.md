# AES S-box based on {AND,XOR,XNOR}

These circuits implement the AES S-box over the basis `{AND, XOR, XNOR, NOT}`.

## S-box with #AND ≤ 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | <ins><strong>29</strong></ins> | 5 | 139 | 35 | 110 | 90 | 20 |
| [circ](aes-sbox-a29-ad5-g184-gd16-xx155-26.circ.txt) | <ins><strong>29</strong></ins> | 5 | 184 | 16 | 155 | 129 | 26 |
| [circ](aes-sbox-a30-ad4-g147-gd36-xx117-13.circ.txt) | 30 | <ins><strong>4</strong></ins> | 147 | 36 | 117 | 104 | 13 |
| [circ](aes-sbox-a30-ad4-g205-gd15-xx175-54.circ.txt) | 30 | <ins><strong>4</strong></ins> | 205 | <ins><strong>15</strong></ins> | 175 | 121 | 54 |
| [circ](aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 32 | 5 | <ins><strong>110</strong></ins> | 23 | 78 | 75 | 3 |
| [circ](aes-sbox-a34-ad4-g110-gd22-xx76-3.circ.txt) | 34 | <ins><strong>4</strong></ins> | <ins><strong>110</strong></ins> | 22 | 76 | 73 | 3 |
| [circ](aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 34 | <ins><strong>4</strong></ins> | 128 | <ins><strong>15</strong></ins> | 94 | 90 | 4 |

## S-box with #AND > 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | <ins><strong>36</strong></ins> | 4 | <ins><strong>138</strong></ins> | 14 | 102 | 80 | 22 |
| [circ](aes-sbox-a39-ad4-g148-gd13-xx109-24.circ.txt) | 39 | 4 | 148 | 13 | 109 | 85 | 24 |
| [circ](aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 45 | 4 | 151 | <ins><strong>12</strong></ins> | 106 | 80 | 26 |
| [circ](aes-sbox-a46-ad3-g250-gd15-xx204-4.circ.txt) | 46 | <ins><strong>3</strong></ins> | 250 | 15 | 204 | 200 | 4 |
| [circ](aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 47 | <ins><strong>3</strong></ins> | 225 | 15 | 178 | 174 | 4 |

Notes: The table presents selected circuits from two disjoint spaces (#AND ≤ 34 and #AND > 34). The circuits with #AND ≤ 31 were obtained from transformations applied to the [AES S-Box circuit](https://github.com/umizame/S-box_29-AND/blob/master/circuits/aes-sbox-fwd-g228-a29-d35-ad6.slp) from [@umizame](https://github.com/umizame/S-box_29-AND), with 29 AND, 195 XOR, 4 NOT, depth 35, and AND-depth 6.

`XX` is the sum of XOR and XNOR gates.
