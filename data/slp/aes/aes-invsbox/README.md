# AES Inverse S-box based on {AND,XOR,XNOR}

The AES S-box is an 8-bit to 8-bit function. The following tables list [selected circuits](../show.md) into four categories: A28, A≤36/AD≥4, A≤36/AD3, GD≤12.

Notation: XX = #XOR + #XNOR; XXN = XX + #NOT.

## Inverse S-box with #AND = 28

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a28-ad4-g140-gd19-xx112-30.circ.txt) | 28 | <ins><strong>4</strong></ins> | 140 | 19 | 112 | 82 | 30 |
| [circ](aes-invsbox-a28-ad4-g150-gd16-xx122-27.circ.txt) | 28 | <ins><strong>4</strong></ins> | 150 | 16 | 122 | 95 | 27 |
| [circ](aes-invsbox-a28-ad4-g156-gd15-xx128-28.circ.txt) | 28 | <ins><strong>4</strong></ins> | 156 | <ins><strong>15</strong></ins> | 128 | 100 | 28 |
| [circ](aes-invsbox-a28-ad5-g129-gd22-xx101-30.circ.txt) | 28 | 5 | <ins><strong>129</strong></ins> | 22 | 101 | 71 | 30 |
| [circ](aes-invsbox-a28-ad5-g131-gd20-xx103-29.circ.txt) | 28 | 5 | 131 | 20 | 103 | 74 | 29 |
| [circ](aes-invsbox-a28-ad5-g147-gd16-xx119-25.circ.txt) | 28 | 5 | 147 | 16 | 119 | 94 | 25 |
| [circ](aes-invsbox-a28-ad5-g152-gd15-xx124-18.circ.txt) | 28 | 5 | 152 | <ins><strong>15</strong></ins> | 124 | 106 | 18 |

Note: The A28 invsbox circuits listed in the table above are either contributions (e.g., AD5/{G152/GD15}) communicated by Milad Nasr (@ Anthropic) on 2026-09-24, or were derived therefrom via linear transformation and optimization of contributed A28 [sbox circuits](../aes-sbox/) or invsbox circuits.

## Inverse S-box with #AND≤36/AD≥4

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a32-ad4-g117-gd21-xx85-10.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 117 | 21 | 85 | 75 | 10 |
| [circ](aes-invsbox-a32-ad4-g119-gd20-xx87-10.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 119 | 20 | 87 | 77 | 10 |
| [circ](aes-invsbox-a32-ad4-g128-gd16-xx96-11.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 128 | 16 | 96 | 85 | 11 |
| [circ](aes-invsbox-a32-ad4-g135-gd15-xx103-11.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 135 | 15 | 103 | 92 | 11 |
| [circ](aes-invsbox-a32-ad5-g107-gd23-xx75-6.circ.txt) | <ins><strong>32</strong></ins> | 5 | 107 | 23 | 75 | 69 | 6 |
| [circ](aes-invsbox-a32-ad5-g110-gd20-xx78-7.circ.txt) | <ins><strong>32</strong></ins> | 5 | 110 | 20 | 78 | 71 | 7 |
| [circ](aes-invsbox-a32-ad5-g125-gd15-xx93-21.circ.txt) | <ins><strong>32</strong></ins> | 5 | 125 | 15 | 93 | 72 | 21 |
| [circ](aes-invsbox-a32-ad7-g105-gd27-xx73-14.circ.txt) | <ins><strong>32</strong></ins> | 7 | <ins><strong>105</strong></ins> | 27 | 73 | 59 | 14 |
| [circ](aes-invsbox-a33-ad4-g107-gd21-xx74-10.circ.txt) | 33 | <ins><strong>4</strong></ins> | 107 | 21 | 74 | 64 | 10 |
| [circ](aes-invsbox-a33-ad4-g108-gd20-xx75-9.circ.txt) | 33 | <ins><strong>4</strong></ins> | 108 | 20 | 75 | 66 | 9 |
| [circ](aes-invsbox-a33-ad4-g118-gd15-xx85-9.circ.txt) | 33 | <ins><strong>4</strong></ins> | 118 | 15 | 85 | 76 | 9 |
| [circ](aes-invsbox-a33-ad4-g123-gd14-xx90-16.circ.txt) | 33 | <ins><strong>4</strong></ins> | 123 | 14 | 90 | 74 | 16 |
| [circ](aes-invsbox-a34-ad4-g136-gd13-xx102-23.circ.txt) | 34 | <ins><strong>4</strong></ins> | 136 | <ins><strong>13</strong></ins> | 102 | 79 | 23 |
| [circ](aes-invsbox-a35-ad4-g133-gd13-xx98-18.circ.txt) | 35 | <ins><strong>4</strong></ins> | 133 | <ins><strong>13</strong></ins> | 98 | 80 | 18 |

## Inverse S-box with #AND≤36/AD3

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a35-ad3-g151-gd22-xx116-30.circ.txt) | 35 | 3 | <ins><strong>151</strong></ins> | 22 | 116 | 86 | 30 |
| [circ](aes-invsbox-a35-ad3-g154-gd20-xx119-30.circ.txt) | 35 | 3 | 154 | 20 | 119 | 89 | 30 |
| [circ](aes-invsbox-a35-ad3-g166-gd16-xx131-30.circ.txt) | 35 | 3 | 166 | 16 | 131 | 101 | 30 |
| [circ](aes-invsbox-a35-ad3-g172-gd15-xx137-30.circ.txt) | 35 | 3 | 172 | <ins><strong>15</strong></ins> | 137 | 107 | 30 |

Note: The A35/AD3 invsbox circuits listed in the table above were derived by linear transformation and optimization of A35/AD3 [sbox circuits](../aes-sbox/) (G139/GD31, G157/GD15, G183/GD13) externally contributed by Milad Nasr (@ Anthropic) on 2026-09-24.

## Inverse S-box with GD≤12

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XXN | #XOR | #XNOR | #NOT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a45-ad4-g141-gd12-xx96-25.circ.txt) | 45 | 4 | <ins><strong>141</strong></ins> | 12 | 96 | 71 | 25 | 0 |
| [circ](aes-invsbox-a55-ad4-g148-gd10-xxn93-26-2.circ.txt) | 55 | 4 | 148 | <ins><strong>10</strong></ins> | 93 | 65 | 26 | 2 |

