# Inverse AES S-box based on {AND,XOR,XNOR}

These circuits implement the Inverse AES S-box over the basis `{AND, XOR, XNOR, NOT}`.

## Inverse S-box with #AND ≤ 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | <ins><strong>29</strong></ins> | 5 | 145 | 32 | 116 | 87 | 29 |
| [circ](aes-invsbox-a29-ad5-g186-gd16-xx157-43.circ.txt) | <ins><strong>29</strong></ins> | 5 | 186 | 16 | 157 | 114 | 43 |
| [circ](aes-invsbox-a30-ad4-g150-gd30-xx120-33.circ.txt) | 30 | <ins><strong>4</strong></ins> | 150 | 30 | 120 | 87 | 33 |
| [circ](aes-invsbox-a30-ad4-g203-gd15-xx173-52.circ.txt) | 30 | <ins><strong>4</strong></ins> | 203 | <ins><strong>15</strong></ins> | 173 | 121 | 52 |
| [circ](aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | 32 | 5 | <ins><strong>112</strong></ins> | 27 | 80 | 71 | 9 |
| [circ](aes-invsbox-a34-ad4-g114-gd27-xx80-11.circ.txt) | 34 | <ins><strong>4</strong></ins> | 114 | 27 | 80 | 69 | 11 |
| [circ](aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | 34 | <ins><strong>4</strong></ins> | 134 | <ins><strong>15</strong></ins> | 100 | 82 | 18 |

## Inverse S-box with #AND > 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | 36 | 4 | 147 | 14 | 111 | 87 | 24 |

Notes: The table presents selected circuits from two disjoint spaces (#AND ≤ 34 and #AND > 34). The circuits with #AND ≤ 31 were obtained from transformations applied to the [AES S-Box circuit](https://github.com/umizame/S-box_29-AND/blob/master/circuits/aes-sbox-fwd-g228-a29-d35-ad6.slp) from [@umizame](https://github.com/umizame/S-box_29-AND), with 29 AND, 195 XOR, 4 NOT, depth 35, and AND-depth 6.

`XX` is the sum of XOR and XNOR gates.
