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
| [circ](aes-sbox-a32-ad5-g107-gd22-xx75-4.circ.txt) | <ins><strong>32</strong></ins> | 5 | 107 | 22 | 75 | 71 | 4 |
| [circ](aes-sbox-a32-ad5-g109-gd20-xx77-4.circ.txt) | <ins><strong>32</strong></ins> | 5 | 109 | 20 | 77 | 73 | 4 |
| [circ](aes-sbox-a32-ad5-g122-gd16-xx90-15.circ.txt) | <ins><strong>32</strong></ins> | 5 | 122 | 16 | 90 | 75 | 15 |
| [circ](aes-sbox-a32-ad7-g103-gd32-xx71-10.circ.txt) | <ins><strong>32</strong></ins> | 7 | <ins><strong>103</strong></ins> | 32 | 71 | 61 | 10 |
| [circ](aes-sbox-a33-ad4-g106-gd22-xx73-5.circ.txt) | 33 | <ins><strong>4</strong></ins> | 106 | 22 | 73 | 68 | 5 |
| [circ](aes-sbox-a33-ad4-g108-gd20-xx75-5.circ.txt) | 33 | <ins><strong>4</strong></ins> | 108 | 20 | 75 | 70 | 5 |
| [circ](aes-sbox-a33-ad4-g116-gd15-xx83-7.circ.txt) | 33 | <ins><strong>4</strong></ins> | 116 | 15 | 83 | 76 | 7 |
| [circ](aes-sbox-a33-ad4-g127-gd14-xx94-23.circ.txt) | 33 | <ins><strong>4</strong></ins> | 127 | 14 | 94 | 71 | 23 |
| [circ](aes-sbox-a34-ad4-g140-gd13-xx106-27.circ.txt) | 34 | <ins><strong>4</strong></ins> | 140 | <ins><strong>13</strong></ins> | 106 | 79 | 27 |
| [circ](aes-sbox-a35-ad4-g133-gd13-xx98-21.circ.txt) | 35 | <ins><strong>4</strong></ins> | 133 | <ins><strong>13</strong></ins> | 98 | 77 | 21 |

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
| [circ](aes-sbox-a45-ad4-g141-gd12-xx96-25.circ.txt) | 45 | 4 | <ins><strong>141</strong></ins> | 12 | 96 | 71 | 25 | 0 |
| [circ](aes-sbox-a49-ad4-g148-gd11-xx99-23.circ.txt) | 49 | 4 | 148 | 11 | 99 | 76 | 23 | 0 |
| [circ](aes-sbox-a59-ad4-g159-gd10-xxn100-30-3.circ.txt) | 59 | 4 | 159 | <ins><strong>10</strong></ins> | 100 | 67 | 30 | 3 |

