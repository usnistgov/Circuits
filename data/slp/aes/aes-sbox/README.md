# AES S-box based on {AND,XOR,XNOR}

The AES S-box is an 8-bit to 8-bit function. The following tables list [selected circuits](../show.md) into four categories: A28, A≤36/AD≥4, A≤36/AD3, GD≤12.

Notation: XX = #XOR + #XNOR; XXN = XX + #NOT.

## S-box with #AND = 28

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a28-ad4-g132-gd22-xx104-40.circ.txt) | 28 | <ins><strong>4</strong></ins> | 132 | 22 | 104 | 64 | 40 |
| [circ](aes-sbox-a28-ad4-g134-gd20-xx106-41.circ.txt) | 28 | <ins><strong>4</strong></ins> | 134 | 20 | 106 | 65 | 41 |
| [circ](aes-sbox-a28-ad4-g141-gd16-xx113-27.circ.txt) | 28 | <ins><strong>4</strong></ins> | 141 | 16 | 113 | 86 | 27 |
| [circ](aes-sbox-a28-ad4-g153-gd15-xx125-29.circ.txt) | 28 | <ins><strong>4</strong></ins> | 153 | 15 | 125 | 96 | 29 |
| [circ](aes-sbox-a28-ad4-g176-gd14-xx148-25.circ.txt) | 28 | <ins><strong>4</strong></ins> | 176 | <ins><strong>14</strong></ins> | 148 | 123 | 25 |
| [circ](aes-sbox-a28-ad5-g128-gd20-xx100-27.circ.txt) | 28 | 5 | <ins><strong>128</strong></ins> | 20 | 100 | 73 | 27 |
| [circ](aes-sbox-a28-ad5-g150-gd15-xx122-23.circ.txt) | 28 | 5 | 150 | 15 | 122 | 99 | 23 |

Note: The A28 sbox circuits listed in the table above are either contributions (e.g., AD4/{G141/GD16}, AD5/{G128/GD20,G150/GD15}) communicated by Milad Nasr (@ Anthropic) on 2026-09-24, or were derived therefrom via linear optimization.

## S-box with #AND≤36/AD≥4

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a32-ad4-g117-gd19-xx85-6.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 117 | 19 | 85 | 79 | 6 |
| [circ](aes-sbox-a32-ad4-g127-gd15-xx95-7.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 127 | 15 | 95 | 88 | 7 |
| [circ](aes-sbox-a32-ad4-g132-gd14-xx100-26.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 132 | 14 | 100 | 74 | 26 |
| [circ](aes-sbox-a32-ad4-g149-gd13-xx117-41.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 149 | <ins><strong>13</strong></ins> | 117 | 76 | 41 |
| [circ](aes-sbox-a32-ad5-g108-gd20-xx76-3.circ.txt) | <ins><strong>32</strong></ins> | 5 | 108 | 20 | 76 | 73 | 3 |
| [circ](aes-sbox-a32-ad5-g119-gd15-xx87-33.circ.txt) | <ins><strong>32</strong></ins> | 5 | 119 | 15 | 87 | 54 | 33 |
| [circ](aes-sbox-a32-ad7-g103-gd32-xx71-10.circ.txt) | <ins><strong>32</strong></ins> | 7 | <ins><strong>103</strong></ins> | 32 | 71 | 61 | 10 |
| [circ](aes-sbox-a33-ad4-g107-gd20-xx74-5.circ.txt) | 33 | <ins><strong>4</strong></ins> | 107 | 20 | 74 | 69 | 5 |
| [circ](aes-sbox-a33-ad4-g115-gd15-xx82-7.circ.txt) | 33 | <ins><strong>4</strong></ins> | 115 | 15 | 82 | 75 | 7 |
| [circ](aes-sbox-a33-ad4-g119-gd14-xx86-25.circ.txt) | 33 | <ins><strong>4</strong></ins> | 119 | 14 | 86 | 61 | 25 |
| [circ](aes-sbox-a33-ad4-g124-gd13-xx91-46.circ.txt) | 33 | <ins><strong>4</strong></ins> | 124 | <ins><strong>13</strong></ins> | 91 | 45 | 46 |
| [circ](aes-sbox-a35-ad4-g123-gd13-xx88-32.circ.txt) | 35 | <ins><strong>4</strong></ins> | 123 | <ins><strong>13</strong></ins> | 88 | 56 | 32 |

## S-box with #AND≤36/AD3

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a35-ad3-g139-gd27-xx104-38.circ.txt) | 35 | 3 | <ins><strong>139</strong></ins> | 27 | 104 | 66 | 38 |
| [circ](aes-sbox-a35-ad3-g141-gd25-xx106-37.circ.txt) | 35 | 3 | 141 | 25 | 106 | 69 | 37 |
| [circ](aes-sbox-a35-ad3-g147-gd20-xx112-36.circ.txt) | 35 | 3 | 147 | 20 | 112 | 76 | 36 |
| [circ](aes-sbox-a35-ad3-g157-gd15-xx122-31.circ.txt) | 35 | 3 | 157 | 15 | 122 | 91 | 31 |
| [circ](aes-sbox-a35-ad3-g172-gd14-xx137-30.circ.txt) | 35 | 3 | 172 | 14 | 137 | 107 | 30 |
| [circ](aes-sbox-a35-ad3-g183-gd13-xx148-29.circ.txt) | 35 | 3 | 183 | <ins><strong>13</strong></ins> | 148 | 119 | 29 |

Note: The A35/AD3 sbox circuits listed in the table above are either (G157/GD15, G183/GD13) contributions communicated by Milad Nasr (@ Anthropic) on 2026-09-24, or were derived therefrom via linear optimization.

## S-box with GD≤12

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XXN | #XOR | #XNOR | #NOT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-a35-ad4-g136-gd12-xx101-33.circ.txt) | 35 | 4 | 136 | 12 | 101 | 68 | 33 | 0 |
| [circ](aes-sbox-a36-ad4-g132-gd12-xx96-32.circ.txt) | 36 | 4 | 132 | 12 | 96 | 64 | 32 | 0 |
| [circ](aes-sbox-a37-ad4-g128-gd12-xx91-29.circ.txt) | 37 | 4 | <ins><strong>128</strong></ins> | 12 | 91 | 62 | 29 | 0 |
| [circ](aes-sbox-a48-ad4-g144-gd11-xx96-34.circ.txt) | 48 | 4 | 144 | 11 | 96 | 62 | 34 | 0 |
| [circ](aes-sbox-a59-ad4-g159-gd10-xxn100-30-3.circ.txt) | 59 | 4 | 159 | <ins><strong>10</strong></ins> | 100 | 67 | 30 | 3 |

