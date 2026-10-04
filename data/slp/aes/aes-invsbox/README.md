# AES Inverse S-box based on {AND,XOR,XNOR}

The tables include selected circuits from two disjoint spaces: #AND ≤ 34 and #AND > 34.

XX denotes the sum of the numbers of XOR and XNOR gates.

## Inverse S-box with #AND ≤ 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a28-ad4-g141-gd31-xx113-28.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 141 | 31 | 113 | 85 | 28 |
| [circ](aes-invsbox-a28-ad4-g174-gd15-xx146-33.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 174 | 15 | 146 | 113 | 33 |
| [circ](aes-invsbox-a28-ad5-g126-gd27-xx98-26.circ.txt) | <ins><strong>28</strong></ins> | 5 | 126 | 27 | 98 | 72 | 26 |
| [circ](aes-invsbox-a28-ad5-g152-gd15-xx124-18.circ.txt) | <ins><strong>28</strong></ins> | 5 | 152 | 15 | 124 | 106 | 18 |
| [circ](aes-invsbox-a32-ad5-g112-gd26-xx80-9.circ.txt) | 32 | 5 | <ins><strong>112</strong></ins> | 26 | 80 | 71 | 9 |
| [circ](aes-invsbox-a34-ad4-g114-gd26-xx80-11.circ.txt) | 34 | <ins><strong>4</strong></ins> | 114 | 26 | 80 | 69 | 11 |
| [circ](aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | 34 | <ins><strong>4</strong></ins> | 134 | 15 | 100 | 82 | 18 |
| [circ](aes-invsbox-a34-ad4-g194-gd14-xx160-20.circ.txt) | 34 | <ins><strong>4</strong></ins> | 194 | <ins><strong>14</strong></ins> | 160 | 140 | 20 |

Note: The A28 invsbox circuits listed in the table above are either (AD5/\{G126/GD27,G152/GD15\}) contributions communicated on 2026-09-24 by Milad Nasr (@ Anthropic), or were derived later on via linear optimization of contributed A28 [sbox circuits](../aes-sbox/).

## Inverse S-box with #AND > 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a35-ad3-g158-gd28-xx123-34.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 158 | 28 | 123 | 89 | 34 |
| [circ](aes-invsbox-a35-ad3-g185-gd15-xx150-26.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 185 | 15 | 150 | 124 | 26 |
| [circ](aes-invsbox-a36-ad4-g140-gd14-xx104-24.circ.txt) | 36 | 4 | <ins><strong>140</strong></ins> | 14 | 104 | 80 | 24 |
| [circ](aes-invsbox-a37-ad4-g172-gd13-xx135-24.circ.txt) | 37 | 4 | 172 | <ins><strong>13</strong></ins> | 135 | 111 | 24 |

Note: The A35/AD3 invsbox circuits listed in the table above were derived by linear transformation and optimization of A35/AD3 [sbox circuits](../aes-sbox/) (G139/GD31, G157/GD15, G183/GD13) externally contributed by Milad Nasr (@ Anthropic) on 2026-09-24.
