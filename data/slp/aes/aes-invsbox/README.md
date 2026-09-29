# Inverse AES S-box based on {AND,XOR,XNOR}

These circuits implement the Inverse AES S-box over the basis `{AND, XOR, XNOR, NOT}`.

## Inverse S-box with #AND ≤ 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a28-ad5-g126-gd27-xx98-26.circ.txt) | <ins><strong>28</strong></ins> | 5 | 126 | 27 | 98 | 72 | 26 |
| [circ](aes-invsbox-a28-ad5-g152-gd15-xx124-18.circ.txt) | <ins><strong>28</strong></ins> | 5 | 152 | <ins><strong>15</strong></ins> | 124 | 106 | 18 |
| [circ](aes-invsbox-a30-ad4-g150-gd30-xx120-33.circ.txt) | 30 | <ins><strong>4</strong></ins> | 150 | 30 | 120 | 87 | 33 |
| [circ](aes-invsbox-a30-ad4-g203-gd15-xx173-52.circ.txt) | 30 | <ins><strong>4</strong></ins> | 203 | <ins><strong>15</strong></ins> | 173 | 121 | 52 |
| [circ](aes-invsbox-a32-ad5-g112-gd26-xx80-9.circ.txt) | 32 | 5 | <ins><strong>112</strong></ins> | 26 | 80 | 71 | 9 |
| [circ](aes-invsbox-a34-ad4-g114-gd26-xx80-11.circ.txt) | 34 | <ins><strong>4</strong></ins> | 114 | 26 | 80 | 69 | 11 |
| [circ](aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | 34 | <ins><strong>4</strong></ins> | 134 | <ins><strong>15</strong></ins> | 100 | 82 | 18 |

## Inverse S-box with #AND > 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a36-ad4-g140-gd14-xx104-24.circ.txt) | <ins><strong>36</strong></ins> | 4 | <ins><strong>140</strong></ins> | <ins><strong>14</strong></ins> | 104 | 80 | 24 |
| [circ](aes-invsbox-a47-ad3-g201-gd25-xx154-3.circ.txt) | 47 | <ins><strong>3</strong></ins> | 201 | 25 | 154 | 151 | 3 |

Notes: The table presents selected circuits from two disjoint spaces (#AND ≤ 34 and #AND > 34). The A28 circuits were externally contributed by Milad Nasr (@ Anthropic) on 2026-09-24. The other circuits with #AND ≤ 31 were obtained from transformations applied to the [AES S-Box circuit](https://github.com/umizame/S-box_29-AND/blob/master/circuits/aes-sbox-fwd-g228-a29-d35-ad6.slp) from [@umizame](https://github.com/umizame/S-box_29-AND), with 29 AND, 195 XOR, 4 NOT, depth 35, and AND-depth 6. `XX` is the sum of XOR and XNOR gates.
