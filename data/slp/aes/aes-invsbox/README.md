# AES Inverse S-box based on {AND,XOR,XNOR}

The AES S-box is an 8-bit to 8-bit function. The following tables list [selected circuits](../show.md) into four categories: A28, A≤36/AD≥4, A≤35/AD3, GD≤12.

Notation: XX = #XOR + #XNOR; XXN = XX + #NOT.

## Inverse S-box with A=28

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a28-ad4-g135-gd20-xx107-35.circ.txt) | 28 | <ins><strong>4</strong></ins> | 135 | 20 | 107 | 72 | 35 |
| [circ](aes-invsbox-a28-ad4-g156-gd15-xx128-28.circ.txt) | 28 | <ins><strong>4</strong></ins> | 156 | <ins><strong>15</strong></ins> | 128 | 100 | 28 |
| [circ](aes-invsbox-a28-ad5-g131-gd20-xx103-29.circ.txt) | 28 | 5 | <ins><strong>131</strong></ins> | 20 | 103 | 74 | 29 |
| [circ](aes-invsbox-a28-ad5-g152-gd15-xx124-18.circ.txt) | 28 | 5 | 152 | <ins><strong>15</strong></ins> | 124 | 106 | 18 |

Note: The A28 invsbox circuits listed in the table above are either contributions (e.g., AD5/{G152/GD15}) communicated by Milad Nasr (@ Anthropic) on 2026-09-24, or were derived therefrom via further optimization of contributed A28 [sbox circuits](../aes-sbox/) or invsbox circuits.

## Inverse S-box with A≤36/AD≥4

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a32-ad4-g117-gd21-xx85-10.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 117 | 21 | 85 | 75 | 10 |
| [circ](aes-invsbox-a32-ad4-g119-gd20-xx87-10.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 119 | 20 | 87 | 77 | 10 |
| [circ](aes-invsbox-a32-ad4-g123-gd15-xx91-21.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 123 | 15 | 91 | 70 | 21 |
| [circ](aes-invsbox-a32-ad4-g125-gd14-xx93-21.circ.txt) | <ins><strong>32</strong></ins> | <ins><strong>4</strong></ins> | 125 | 14 | 93 | 72 | 21 |
| [circ](aes-invsbox-a32-ad5-g110-gd20-xx78-7.circ.txt) | <ins><strong>32</strong></ins> | 5 | 110 | 20 | 78 | 71 | 7 |
| [circ](aes-invsbox-a32-ad5-g119-gd15-xx87-27.circ.txt) | <ins><strong>32</strong></ins> | 5 | 119 | 15 | 87 | 60 | 27 |
| [circ](aes-invsbox-a32-ad7-g105-gd27-xx73-14.circ.txt) | <ins><strong>32</strong></ins> | 7 | <ins><strong>105</strong></ins> | 27 | 73 | 59 | 14 |
| [circ](aes-invsbox-a33-ad4-g108-gd20-xx75-9.circ.txt) | 33 | <ins><strong>4</strong></ins> | 108 | 20 | 75 | 66 | 9 |
| [circ](aes-invsbox-a33-ad4-g115-gd15-xx82-37.circ.txt) | 33 | <ins><strong>4</strong></ins> | 115 | 15 | 82 | 45 | 37 |
| [circ](aes-invsbox-a33-ad4-g117-gd14-xx84-24.circ.txt) | 33 | <ins><strong>4</strong></ins> | 117 | 14 | 84 | 60 | 24 |
| [circ](aes-invsbox-a33-ad4-g123-gd13-xx90-34.circ.txt) | 33 | <ins><strong>4</strong></ins> | 123 | <ins><strong>13</strong></ins> | 90 | 56 | 34 |
| [circ](aes-invsbox-a34-ad4-g119-gd13-xx85-17.circ.txt) | 34 | <ins><strong>4</strong></ins> | 119 | <ins><strong>13</strong></ins> | 85 | 68 | 17 |

## Inverse S-box with A≤35/AD3

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a35-ad3-g151-gd22-xx116-30.circ.txt) | 35 | 3 | <ins><strong>151</strong></ins> | 22 | 116 | 86 | 30 |
| [circ](aes-invsbox-a35-ad3-g154-gd20-xx119-30.circ.txt) | 35 | 3 | 154 | 20 | 119 | 89 | 30 |
| [circ](aes-invsbox-a35-ad3-g157-gd15-xx122-27.circ.txt) | 35 | 3 | 157 | <ins><strong>15</strong></ins> | 122 | 95 | 27 |

Note: The A35/AD3 invsbox circuits listed in the table above were derived by optimization of A35/AD3 [sbox circuits](../aes-sbox/) externally contributed by Milad Nasr (@ Anthropic) on 2026-09-24.

## Inverse S-box with GD≤12

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XXN | #XOR | #XNOR | #NOT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-invsbox-a36-ad4-g133-gd12-xx97-31.circ.txt) | 36 | 4 | 133 | 12 | 97 | 66 | 31 | 0 |
| [circ](aes-invsbox-a39-ad4-g127-gd12-xx88-25.circ.txt) | 39 | 4 | <ins><strong>127</strong></ins> | 12 | 88 | 63 | 25 | 0 |
| [circ](aes-invsbox-a48-ad4-g139-gd11-xx91-22.circ.txt) | 48 | 4 | 139 | 11 | 91 | 69 | 22 | 0 |
| [circ](aes-invsbox-a55-ad4-g148-gd10-xxn93-26-2.circ.txt) | 55 | 4 | 148 | <ins><strong>10</strong></ins> | 93 | 65 | 26 | 2 |

