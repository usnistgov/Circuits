# AES S-box based on {AND,XOR,XNOR}

The tables include selected circuits from two disjoint spaces: #AND ≤ 34 and #AND > 34.

XX denotes the sum of the numbers of XOR and XNOR gates.

## S-box with #AND ≤ 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a28-ad4-g131-gd32-xx103-47.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 131 | 32 | 103 | 56 | 47 |
| [circ](aes-sbox-a28-ad4-g154-gd15-xx126-28.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 154 | 15 | 126 | 98 | 28 |
| [circ](aes-sbox-a28-ad4-g177-gd14-xx149-24.circ.txt) | <ins><strong>28</strong></ins> | <ins><strong>4</strong></ins> | 177 | <ins><strong>14</strong></ins> | 149 | 125 | 24 |
| [circ](aes-sbox-a28-ad5-g124-gd27-xx96-26.circ.txt) | <ins><strong>28</strong></ins> | 5 | 124 | 27 | 96 | 70 | 26 |
| [circ](aes-sbox-a28-ad5-g150-gd15-xx122-23.circ.txt) | <ins><strong>28</strong></ins> | 5 | 150 | 15 | 122 | 99 | 23 |
| [circ](aes-sbox-a32-ad5-g110-gd22-xx78-3.circ.txt) | 32 | 5 | <ins><strong>110</strong></ins> | 22 | 78 | 75 | 3 |
| [circ](aes-sbox-a34-ad4-g110-gd22-xx76-3.circ.txt) | 34 | <ins><strong>4</strong></ins> | <ins><strong>110</strong></ins> | 22 | 76 | 73 | 3 |
| [circ](aes-sbox-a34-ad4-g125-gd15-xx91-4.circ.txt) | 34 | <ins><strong>4</strong></ins> | 125 | 15 | 91 | 87 | 4 |

Note: The A28 sbox circuits listed in the table above are either (AD4/\{G154/GD15, G177/GD14\}, AD5\{G124/GD27,G150/GD15\}) contributions communicated by Milad Nasr (@ Anthropic) on 2026-09-24, or were derived therefrom via linear optimization.


## S-box with #AND > 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a35-ad3-g139-gd27-xx104-38.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 139 | 27 | 104 | 66 | 38 |
| [circ](aes-sbox-a35-ad3-g157-gd15-xx122-31.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 157 | 15 | 122 | 91 | 31 |
| [circ](aes-sbox-a35-ad3-g183-gd13-xx148-29.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 183 | 13 | 148 | 119 | 29 |
| [circ](aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 36 | 4 | <ins><strong>138</strong></ins> | 14 | 102 | 80 | 22 |
| [circ](aes-sbox-a39-ad4-g148-gd13-xx109-24.circ.txt) | 39 | 4 | 148 | 13 | 109 | 85 | 24 |
| [circ](aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 45 | 4 | 151 | 12 | 106 | 80 | 26 |
| [circ](aes-sbox-a53-ad4-g157-gd11-xx104-31.circ.txt) | 53 | 4 | 157 | <ins><strong>11</strong></ins> | 104 | 73 | 31 |

Note: The A35/AD3 sbox circuits listed in the table above are either (G157/GD15, G183/GD13) contributions communicated by Milad Nasr (@ Anthropic) on 2026-09-24, or were derived therefrom via linear optimization.
