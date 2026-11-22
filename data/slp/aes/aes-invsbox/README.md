# AES Inverse S-box based on {AND,XOR,XNOR}

The AES S-box is an 8-bit to 8-bit function. The following tables list [selected circuits](../show.md) into four categories: A28, A≤36/AD≥4, A≤36/AD3, GD≤12.

Notation: XX = #XOR + #XNOR; XXN = XX + #NOT.

## Inverse S-box with #AND = 28

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a28-ad4-g135-gd20-xx107-35.circ.txt) | 28 | <ins><strong>4</strong></ins> | 135 | 20 | 107 | 72 | 35 |
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
| [circ](aes-invsbox-a32-ad4-g124-gd15-xx92-37.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 124 | 15 | 92 | 55 | 37 |
| [circ](aes-invsbox-a32-ad4-g131-gd14-xx99-36.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 131 | 14 | 99 | 63 | 36 |
| [circ](aes-invsbox-a32-ad5-g107-gd23-xx75-6.circ.txt) | <ins><strong>32</strong></ins> | 5 | 107 | 23 | 75 | 69 | 6 |
| [circ](aes-invsbox-a32-ad5-g110-gd20-xx78-7.circ.txt) | <ins><strong>32</strong></ins> | 5 | 110 | 20 | 78 | 71 | 7 |
| [circ](aes-invsbox-a32-ad5-g121-gd15-xx89-21.circ.txt) | <ins><strong>32</strong></ins> | 5 | 121 | 15 | 89 | 68 | 21 |
| [circ](aes-invsbox-a32-ad7-g105-gd27-xx73-14.circ.txt) | <ins><strong>32</strong></ins> | 7 | <ins><strong>105</strong></ins> | 27 | 73 | 59 | 14 |
| [circ](aes-invsbox-a33-ad4-g107-gd21-xx74-10.circ.txt) | 33 | <ins><strong>4</strong></ins> | 107 | 21 | 74 | 64 | 10 |
| [circ](aes-invsbox-a33-ad4-g108-gd20-xx75-9.circ.txt) | 33 | <ins><strong>4</strong></ins> | 108 | 20 | 75 | 66 | 9 |
| [circ](aes-invsbox-a33-ad4-g115-gd15-xx82-37.circ.txt) | 33 | <ins><strong>4</strong></ins> | 115 | 15 | 82 | 45 | 37 |
| [circ](aes-invsbox-a33-ad4-g117-gd14-xx84-24.circ.txt) | 33 | <ins><strong>4</strong></ins> | 117 | 14 | 84 | 60 | 24 |
| [circ](aes-invsbox-a33-ad4-g125-gd13-xx92-44.circ.txt) | 33 | <ins><strong>4</strong></ins> | 125 | <ins><strong>13</strong></ins> | 92 | 48 | 44 |

## Inverse S-box with #AND≤36/AD3

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a35-ad3-g151-gd22-xx116-30.circ.txt) | 35 | 3 | <ins><strong>151</strong></ins> | 22 | 116 | 86 | 30 |
| [circ](aes-invsbox-a35-ad3-g154-gd20-xx119-30.circ.txt) | 35 | 3 | 154 | 20 | 119 | 89 | 30 |
| [circ](aes-invsbox-a35-ad3-g157-gd15-xx122-27.circ.txt) | 35 | 3 | 157 | <ins><strong>15</strong></ins> | 122 | 95 | 27 |

Note: The A35/AD3 invsbox circuits listed in the table above were derived by linear transformation and optimization of A35/AD3 [sbox circuits](../aes-sbox/) (G139/GD31, G157/GD15, G183/GD13) externally contributed by Milad Nasr (@ Anthropic) on 2026-09-24.

## Inverse S-box with GD≤12

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XXN | #XOR | #XNOR | #NOT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a35-ad4-g140-gd12-xx105-40.circ.txt) | 35 | 4 | 140 | 12 | 105 | 65 | 40 | 0 |
| [circ](aes-invsbox-a36-ad4-g135-gd12-xx99-34.circ.txt) | 36 | 4 | 135 | 12 | 99 | 65 | 34 | 0 |
| [circ](aes-invsbox-a37-ad4-g133-gd12-xx96-27.circ.txt) | 37 | 4 | <ins><strong>133</strong></ins> | 12 | 96 | 69 | 27 | 0 |
| [circ](aes-invsbox-a55-ad4-g148-gd10-xxn93-26-2.circ.txt) | 55 | 4 | 148 | <ins><strong>10</strong></ins> | 93 | 65 | 26 | 2 |

