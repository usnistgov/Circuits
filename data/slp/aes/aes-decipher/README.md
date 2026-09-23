# AES Decipher (flat circuits)

TM = Tuple metric. TM1 = A-AD-G-GD; TM2 = AD-GD-G-A; TM3 = GD-G-AD-A; TM4 = G-A-GD-AD.

A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR.

Within each table, bold underlined values mark the lowest displayed value in selected columns.

## AES-128

### Decipher using [Inv]Sbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes128-decipher-a5800-ad100-g35552-gd784-xx29752-5456.circ.txt) | <ins><strong>5800</strong></ins> | 100 | 35552 | 784 | 29752 | 24296 | 5456 | [circ](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](aes128-decipher-a6800-ad80-g34504-gd396-xx27704-3056.circ.txt) | 6800 | <ins><strong>80</strong></ins> | 34504 | <ins><strong>396</strong></ins> | 27704 | 24648 | 3056 | [circ](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](aes128-decipher-a6400-ad100-g29112-gd623-xx22712-1576.circ.txt) | 6400 | 100 | <ins><strong>29112</strong></ins> | 623 | 22712 | 21136 | 1576 | [circ](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

### Decipher using [Inv]Sbox with A>34

| File | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes128-decipher-a7200-ad80-g35832-gd403-xx28632-4736.circ.txt) | <ins><strong>7200</strong></ins> | 80 | <ins><strong>35832</strong></ins> | 403 | 28632 | 23896 | 4736 | [circ](../aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ](aes128-decipher-a7640-ad70-g40464-gd386-xx32824-4016.circ.txt) | 7640 | <ins><strong>70</strong></ins> | 40464 | 386 | 32824 | 28808 | 4016 | [circ](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ](aes128-decipher-a7560-ad80-g37504-gd356-xx29944-4896.circ.txt) | 7560 | 80 | 37504 | <ins><strong>356</strong></ins> | 29944 | 25048 | 4896 | [circ](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 3 |

## AES-192

### Decipher using [Inv]Sbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes192-decipher-a6496-ad100-g40440-gd804-xx33944-6216.circ.txt) | <ins><strong>6496</strong></ins> | 100 | 40440 | 804 | 33944 | 27728 | 6216 | [circ](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](aes192-decipher-a7616-ad80-g39384-gd414-xx31768-3592.circ.txt) | 7616 | <ins><strong>80</strong></ins> | 39384 | <ins><strong>414</strong></ins> | 31768 | 28176 | 3592 | [circ](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](aes192-decipher-a7168-ad100-g33176-gd655-xx26008-1832.circ.txt) | 7168 | 100 | <ins><strong>33176</strong></ins> | 655 | 26008 | 24176 | 1832 | [circ](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

### Decipher using [Inv]Sbox with A>34

| File | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes192-decipher-a8064-ad80-g40792-gd427-xx32728-5320.circ.txt) | <ins><strong>8064</strong></ins> | 80 | <ins><strong>40792</strong></ins> | 427 | 32728 | 27408 | 5320 | [circ](../aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ](aes192-decipher-a8416-ad72-g44984-gd402-xx36568-4744.circ.txt) | 8416 | <ins><strong>72</strong></ins> | 44984 | 402 | 36568 | 31824 | 4744 | [circ](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ](aes192-decipher-a8352-ad80-g42616-gd378-xx34264-5448.circ.txt) | 8352 | 80 | 42616 | <ins><strong>378</strong></ins> | 34264 | 28816 | 5448 | [circ](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 3 |

## AES-256

### Decipher using [Inv]Sbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes256-decipher-a8004-ad135-g49220-gd1062-xx41216-7543.circ.txt) | <ins><strong>8004</strong></ins> | 135 | 49220 | 1062 | 41216 | 33673 | 7543 | [circ](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](aes256-decipher-a9384-ad108-g47848-gd537-xx38464-4247.circ.txt) | 9384 | <ins><strong>108</strong></ins> | 47848 | <ins><strong>537</strong></ins> | 38464 | 34217 | 4247 | [circ](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](aes256-decipher-a8832-ad135-g40320-gd848-xx31488-2179.circ.txt) | 8832 | 135 | <ins><strong>40320</strong></ins> | 848 | 31488 | 29309 | 2179 | [circ](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

### Decipher using [Inv]Sbox with A>34

| File | A | AD | G | GD | XX | X | X' | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes256-decipher-a9936-ad108-g49616-gd549-xx39680-6527.circ.txt) | <ins><strong>9936</strong></ins> | 108 | <ins><strong>49616</strong></ins> | 549 | 39680 | 33153 | 6527 | [circ](../aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ](aes256-decipher-a10508-ad95-g55804-gd523-xx45296-5591.circ.txt) | 10508 | <ins><strong>95</strong></ins> | 55804 | 523 | 45296 | 39705 | 5591 | [circ](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ](aes256-decipher-a10404-ad108-g51956-gd484-xx41552-6735.circ.txt) | 10404 | 108 | 51956 | <ins><strong>484</strong></ins> | 41552 | 34817 | 6735 | [circ](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 3 |
