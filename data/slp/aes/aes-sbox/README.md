# AES S-box based on {AND,XOR,XNOR}

The tables include selected circuits from three disjoint spaces: #AND = 28, 32 ≤ #AND ≤ 34, and #AND > 34.

XX denotes the sum of the numbers of XOR and XNOR gates.

## S-box with #AND = 28

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a28-ad4-g131-gd31-xx103-43.circ.txt) | 28 | <ins><strong>4</strong></ins> | 131 | 31 | 103 | 60 | 43 |
| [circ](aes-sbox-a28-ad4-g134-gd20-xx106-41.circ.txt) | 28 | <ins><strong>4</strong></ins> | 134 | 20 | 106 | 65 | 41 |
| [circ](aes-sbox-a28-ad4-g141-gd16-xx113-27.circ.txt) | 28 | <ins><strong>4</strong></ins> | 141 | 16 | 113 | 86 | 27 |
| [circ](aes-sbox-a28-ad4-g153-gd15-xx125-29.circ.txt) | 28 | <ins><strong>4</strong></ins> | 153 | 15 | 125 | 96 | 29 |
| [circ](aes-sbox-a28-ad4-g177-gd14-xx149-24.circ.txt) | 28 | <ins><strong>4</strong></ins> | 177 | <ins><strong>14</strong></ins> | 149 | 125 | 24 |
| [circ](aes-sbox-a28-ad5-g124-gd27-xx96-26.circ.txt) | 28 | 5 | <ins><strong>124</strong></ins> | 27 | 96 | 70 | 26 |
| [circ](aes-sbox-a28-ad5-g128-gd20-xx100-27.circ.txt) | 28 | 5 | 128 | 20 | 100 | 73 | 27 |
| [circ](aes-sbox-a28-ad5-g150-gd15-xx122-23.circ.txt) | 28 | 5 | 150 | 15 | 122 | 99 | 23 |

Note: The A28 sbox circuits listed in the table above are either (AD4/\{G154/GD15, G177/GD14\}, AD5\{G124/GD27,G150/GD15\}) contributions communicated by Milad Nasr (@ Anthropic) on 2026-09-24, or were derived therefrom via linear optimization.


## S-box with 32 ≤ #AND ≤ 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a32-ad4-g137-gd18-xx105-7.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 137 | 18 | 105 | 98 | 7 |
| [circ](aes-sbox-a32-ad5-g104-gd30-xx72-5.circ.txt) | <ins><strong>32</strong></ins> | 5 | 104 | 30 | 72 | 67 | 5 |
| [circ](aes-sbox-a32-ad5-g105-gd26-xx73-5.circ.txt) | <ins><strong>32</strong></ins> | 5 | 105 | 26 | 73 | 68 | 5 |
| [circ](aes-sbox-a32-ad5-g109-gd20-xx77-4.circ.txt) | <ins><strong>32</strong></ins> | 5 | 109 | 20 | 77 | 73 | 4 |
| [circ](aes-sbox-a32-ad5-g128-gd15-xx96-21.circ.txt) | <ins><strong>32</strong></ins> | 5 | 128 | <ins><strong>15</strong></ins> | 96 | 75 | 21 |
| [circ](aes-sbox-a32-ad7-g103-gd32-xx71-10.circ.txt) | <ins><strong>32</strong></ins> | 7 | <ins><strong>103</strong></ins> | 32 | 71 | 61 | 10 |
| [circ](aes-sbox-a33-ad4-g105-gd25-xx72-6.circ.txt) | 33 | <ins><strong>4</strong></ins> | 105 | 25 | 72 | 66 | 6 |
| [circ](aes-sbox-a33-ad4-g108-gd20-xx75-5.circ.txt) | 33 | <ins><strong>4</strong></ins> | 108 | 20 | 75 | 70 | 5 |
| [circ](aes-sbox-a33-ad4-g116-gd15-xx83-7.circ.txt) | 33 | <ins><strong>4</strong></ins> | 116 | <ins><strong>15</strong></ins> | 83 | 76 | 7 |

## S-box with #AND > 34

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a35-ad3-g139-gd27-xx104-38.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 139 | 27 | 104 | 66 | 38 |
| [circ](aes-sbox-a35-ad3-g157-gd15-xx122-31.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 157 | 15 | 122 | 91 | 31 |
| [circ](aes-sbox-a35-ad3-g183-gd13-xx148-29.circ.txt) | <ins><strong>35</strong></ins> | <ins><strong>3</strong></ins> | 183 | 13 | 148 | 119 | 29 |
| [circ](aes-sbox-a36-ad4-g137-gd14-xx101-22.circ.txt) | 36 | 4 | <ins><strong>137</strong></ins> | 14 | 101 | 79 | 22 |
| [circ](aes-sbox-a38-ad4-g148-gd13-xx110-24.circ.txt) | 38 | 4 | 148 | 13 | 110 | 86 | 24 |
| [circ](aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 45 | 4 | 151 | 12 | 106 | 80 | 26 |
| [circ](aes-sbox-a49-ad4-g164-gd11-xx115-32.circ.txt) | 49 | 4 | 164 | <ins><strong>11</strong></ins> | 115 | 83 | 32 |
| [circ](aes-sbox-a51-ad4-g157-gd11-xx106-28.circ.txt) | 51 | 4 | 157 | <ins><strong>11</strong></ins> | 106 | 78 | 28 |

Note: The A35/AD3 sbox circuits listed in the table above are either (G157/GD15, G183/GD13) contributions communicated by Milad Nasr (@ Anthropic) on 2026-09-24, or were derived therefrom via linear optimization.
