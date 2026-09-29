# AES S-box based on {AND,XOR,XNOR}

These circuits implement the AES S-box over the basis `{AND, XOR, XNOR, NOT}`.

## S-box with #AND ≤ 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a28-ad4-g131-gd36-xx103-28.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 131 | 36 | 103 | 75 | 28 |
| [circ](aes-sbox-a28-ad4-g154-gd15-xx126-28.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 154 | 15 | 126 | 98 | 28 |
| [circ](aes-sbox-a28-ad4-g177-gd14-xx149-24.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 177 | <ins><strong>14</strong></ins> | 149 | 125 | 24 |
| [circ](aes-sbox-a28-ad5-g124-gd27-xx96-26.circ.txt) | <ins><strong>28</strong></ins> | 5 | 124 | 27 | 96 | 70 | 26 |
| [circ](aes-sbox-a28-ad5-g150-gd15-xx122-23.circ.txt) | <ins><strong>28</strong></ins> | 5 | 150 | 15 | 122 | 99 | 23 |
| [circ](aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 32 | 5 | <ins><strong>110</strong></ins> | 23 | 78 | 75 | 3 |
| [circ](aes-sbox-a34-ad4-g110-gd22-xx76-3.circ.txt) | 34 | <ins><strong>4</strong></ins> | <ins><strong>110</strong></ins> | 22 | 76 | 73 | 3 |
| [circ](aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 34 | <ins><strong>4</strong></ins> | 128 | 15 | 94 | 90 | 4 |

## S-box with #AND > 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a35-ad3-g139-gd31-xx104-25.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 139 | 31 | 104 | 79 | 25 |
| [circ](aes-sbox-a35-ad3-g157-gd15-xx122-31.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 157 | 15 | 122 | 91 | 31 |
| [circ](aes-sbox-a35-ad3-g183-gd13-xx148-29.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 183 | 13 | 148 | 119 | 29 |
| [circ](aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 36 | 4 | <ins><strong>138</strong></ins> | 14 | 102 | 80 | 22 |
| [circ](aes-sbox-a39-ad4-g148-gd13-xx109-24.circ.txt) | 39 | 4 | 148 | 13 | 109 | 85 | 24 |
| [circ](aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 45 | 4 | 151 | 12 | 106 | 80 | 26 |
| [circ](aes-sbox-a53-ad4-g157-gd11-xx104-31.circ.txt) | 53 | 4 | 157 | <ins><strong>11</strong></ins> | 104 | 73 | 31 |

Notes: The table presents selected circuits from two disjoint spaces (#AND ≤ 34 and #AND > 34). The A28 and A35/AD3 circuits were externally contributed by Milad Nasr (@ Anthropic) on 2026-09-24. `XX` is the sum of XOR and XNOR gates.
