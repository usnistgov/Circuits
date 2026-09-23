# AES KeyExpansion (flat circuits)

TM = Tuple metric. TM1 = A-AD-G-GD; TM2 = AD-GD-G-A; TM3 = GD-G-AD-A; TM4 = G-A-GD-AD.

A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR.

Within each table, bold underlined values mark the lowest displayed value in selected columns.

## AES-128

### KeyExpansion using Sbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|:---:|
| [circ](aes128-keyexp-a1160-ad50-g6840-gd381-xx5680-816.circ.txt) | <ins><strong>1160</strong></ins> | 50 | 6840 | 381 | 5680 | 4864 | 816 | [circ](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | 1 |
| [circ](aes128-keyexp-a1360-ad40-g6400-gd190-xx5040-176.circ.txt) | 1360 | <ins><strong>40</strong></ins> | 6400 | <ins><strong>190</strong></ins> | 5040 | 4864 | 176 | [circ](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 2,3 |
| [circ](aes128-keyexp-a1280-ad50-g5680-gd270-xx4400-136.circ.txt) | 1280 | 50 | <ins><strong>5680</strong></ins> | 270 | 4400 | 4264 | 136 | [circ](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 4 |

### KeyExpansion using Sbox with A>34

| File | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|:---:|
| [circ](aes128-keyexp-a1440-ad40-g6800-gd180-xx5360-896.circ.txt) | <ins><strong>1440</strong></ins> | 40 | <ins><strong>6800</strong></ins> | 180 | 5360 | 4464 | 896 | [circ](../aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 1,4 |
| [circ](aes128-keyexp-a1880-ad30-g10280-gd190-xx8400-176.circ.txt) | 1880 | <ins><strong>30</strong></ins> | 10280 | 190 | 8400 | 8224 | 176 | [circ](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 2 |
| [circ](aes128-keyexp-a1800-ad40-g7320-gd160-xx5520-1056.circ.txt) | 1800 | 40 | 7320 | <ins><strong>160</strong></ins> | 5520 | 4464 | 1056 | [circ](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 3 |

## AES-192

### KeyExpansion using Sbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|:---:|
| [circ](aes192-keyexp-a928-ad40-g5920-gd319-xx4992-648.circ.txt) | <ins><strong>928</strong></ins> | 40 | 5920 | 319 | 4992 | 4344 | 648 | [circ](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | 1 |
| [circ](aes192-keyexp-a1088-ad32-g5568-gd166-xx4480-136.circ.txt) | 1088 | <ins><strong>32</strong></ins> | 5568 | <ins><strong>166</strong></ins> | 4480 | 4344 | 136 | [circ](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 2,3 |
| [circ](aes192-keyexp-a1024-ad40-g4992-gd230-xx3968-104.circ.txt) | 1024 | 40 | <ins><strong>4992</strong></ins> | 230 | 3968 | 3864 | 104 | [circ](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 4 |

### KeyExpansion using Sbox with A>34

| File | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|:---:|
| [circ](aes192-keyexp-a1152-ad32-g5888-gd158-xx4736-712.circ.txt) | <ins><strong>1152</strong></ins> | 32 | <ins><strong>5888</strong></ins> | 158 | 4736 | 4024 | 712 | [circ](../aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 1,4 |
| [circ](aes192-keyexp-a1504-ad24-g8672-gd166-xx7168-136.circ.txt) | 1504 | <ins><strong>24</strong></ins> | 8672 | 166 | 7168 | 7032 | 136 | [circ](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 2 |
| [circ](aes192-keyexp-a1440-ad32-g6304-gd142-xx4864-840.circ.txt) | 1440 | 32 | 6304 | <ins><strong>142</strong></ins> | 4864 | 4024 | 840 | [circ](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 3 |

## AES-256

### KeyExpansion using Sbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|:---:|
| [circ](aes256-keyexp-a1508-ad65-g8892-gd495-xx7384-1047.circ.txt) | <ins><strong>1508</strong></ins> | 65 | 8892 | 495 | 7384 | 6337 | 1047 | [circ](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | 1 |
| [circ](aes256-keyexp-a1768-ad52-g8320-gd247-xx6552-215.circ.txt) | 1768 | <ins><strong>52</strong></ins> | 8320 | <ins><strong>247</strong></ins> | 6552 | 6337 | 215 | [circ](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | 2,3 |
| [circ](aes256-keyexp-a1664-ad65-g7384-gd351-xx5720-163.circ.txt) | 1664 | 65 | <ins><strong>7384</strong></ins> | 351 | 5720 | 5557 | 163 | [circ](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | 4 |

### KeyExpansion using Sbox with A>34

| File | A | AD | G | GD | XX | X | X' | Sbox | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|:---:|
| [circ](aes256-keyexp-a1872-ad52-g8840-gd234-xx6968-1151.circ.txt) | <ins><strong>1872</strong></ins> | 52 | <ins><strong>8840</strong></ins> | 234 | 6968 | 5817 | 1151 | [circ](../aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | 1,4 |
| [circ](aes256-keyexp-a2444-ad39-g13364-gd247-xx10920-215.circ.txt) | 2444 | <ins><strong>39</strong></ins> | 13364 | 247 | 10920 | 10705 | 215 | [circ](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | 2 |
| [circ](aes256-keyexp-a2340-ad52-g9516-gd208-xx7176-1359.circ.txt) | 2340 | 52 | 9516 | <ins><strong>208</strong></ins> | 7176 | 5817 | 1359 | [circ](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | 3 |
