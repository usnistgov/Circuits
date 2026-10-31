# AES Inverse S-box based on {AND,XOR,XNOR}

The tables include selected circuits from three disjoint spaces: #AND = 28, 32 ≤ #AND ≤ 34, and #AND > 34.

XX denotes the sum of the numbers of XOR and XNOR gates.

## Inverse S-box with #AND = 28

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a28-ad4-g140-gd30-xx112-28.circ.txt) | 28 | <ins><strong>4</strong></ins> | 140 | 30 | 112 | 84 | 28 |
| [circ](aes-invsbox-a28-ad4-g147-gd19-xx119-26.circ.txt) | 28 | <ins><strong>4</strong></ins> | 147 | 19 | 119 | 93 | 26 |
| [circ](aes-invsbox-a28-ad4-g167-gd15-xx139-31.circ.txt) | 28 | <ins><strong>4</strong></ins> | 167 | <ins><strong>15</strong></ins> | 139 | 108 | 31 |
| [circ](aes-invsbox-a28-ad5-g126-gd27-xx98-26.circ.txt) | 28 | 5 | <ins><strong>126</strong></ins> | 27 | 98 | 72 | 26 |
| [circ](aes-invsbox-a28-ad5-g131-gd20-xx103-29.circ.txt) | 28 | 5 | 131 | 20 | 103 | 74 | 29 |
| [circ](aes-invsbox-a28-ad5-g152-gd15-xx124-18.circ.txt) | 28 | 5 | 152 | <ins><strong>15</strong></ins> | 124 | 106 | 18 |

Note: The A28 invsbox circuits listed in the table above are either (AD5/\{G126/GD27,G152/GD15\}) contributions communicated on 2026-09-24 by Milad Nasr (@ Anthropic), or were derived later on via linear optimization of contributed A28 [sbox circuits](../aes-sbox/).

## Inverse S-box with 32 ≤ #AND ≤ 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a32-ad5-g106-gd25-xx74-9.circ.txt) | <ins><strong>32</strong></ins> | 5 | 106 | 25 | 74 | 65 | 9 |
| [circ](aes-invsbox-a32-ad5-g112-gd20-xx80-7.circ.txt) | <ins><strong>32</strong></ins> | 5 | 112 | 20 | 80 | 73 | 7 |
| [circ](aes-invsbox-a32-ad5-g128-gd16-xx96-23.circ.txt) | <ins><strong>32</strong></ins> | 5 | 128 | 16 | 96 | 73 | 23 |
| [circ](aes-invsbox-a32-ad7-g105-gd27-xx73-14.circ.txt) | <ins><strong>32</strong></ins> | 7 | <ins><strong>105</strong></ins> | 27 | 73 | 59 | 14 |
| [circ](aes-invsbox-a33-ad4-g107-gd21-xx74-10.circ.txt) | 33 | <ins><strong>4</strong></ins> | 107 | 21 | 74 | 64 | 10 |
| [circ](aes-invsbox-a33-ad4-g120-gd15-xx87-8.circ.txt) | 33 | <ins><strong>4</strong></ins> | 120 | 15 | 87 | 79 | 8 |
| [circ](aes-invsbox-a34-ad4-g139-gd14-xx105-23.circ.txt) | 34 | <ins><strong>4</strong></ins> | 139 | <ins><strong>14</strong></ins> | 105 | 82 | 23 |

## Inverse S-box with #AND > 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a35-ad3-g149-gd26-xx114-32.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 149 | 26 | 114 | 82 | 32 |
| [circ](aes-invsbox-a35-ad3-g155-gd20-xx120-29.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 155 | 20 | 120 | 91 | 29 |
| [circ](aes-invsbox-a35-ad3-g178-gd15-xx143-30.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 178 | 15 | 143 | 113 | 30 |
| [circ](aes-invsbox-a35-ad3-g190-gd14-xx155-27.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 190 | 14 | 155 | 128 | 27 |
| [circ](aes-invsbox-a36-ad4-g138-gd14-xx102-24.circ.txt) | 36 | 4 | <ins><strong>138</strong></ins> | 14 | 102 | 78 | 24 |
| [circ](aes-invsbox-a37-ad4-g157-gd13-xx120-19.circ.txt) | 37 | 4 | 157 | <ins><strong>13</strong></ins> | 120 | 101 | 19 |

Note: The A35/AD3 invsbox circuits listed in the table above were derived by linear transformation and optimization of A35/AD3 [sbox circuits](../aes-sbox/) (G139/GD31, G157/GD15, G183/GD13) externally contributed by Milad Nasr (@ Anthropic) on 2026-09-24.
