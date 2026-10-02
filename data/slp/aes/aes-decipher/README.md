# AES Decipher (flat circuits)

TM = Tuple metric. TM1 = A-AD-G-GD; TM2 = AD-GD-G-A; TM3 = GD-G-AD-A; TM4 = G-A-GD-AD.

A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR.

Within each table, bold underlined values mark the lowest displayed value in selected columns.

## AES-128

### Decipher using [Inv]Sbox with A<=34

| File | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes128-decipher-a5800-ad100-g35552-gd784-xx29752-5456.circ.txt) | <ins><strong>5800</strong></ins> | 100 | 35552 | 784 | 29752 | [A29/AD5/G139/GD35](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [A29/AD5/G145/GD32](../aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [G114/GD8](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](aes128-decipher-a6800-ad80-g34504-gd396-xx27704-3056.circ.txt) | 6800 | <ins><strong>80</strong></ins> | 34504 | <ins><strong>396</strong></ins> | 27704 | [A34/AD4/G128/GD15](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [A34/AD4/G134/GD15](../aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [G146/GD5](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](aes128-decipher-a6400-ad100-g29112-gd623-xx22712-1576.circ.txt) | 6400 | 100 | <ins><strong>29112</strong></ins> | 623 | 22712 | [A32/AD5/G110/GD23](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [A32/AD5/G112/GD27](../aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [G114/GD8](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

### Decipher using [Inv]Sbox with A>34

| File | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes128-decipher-a7640-ad70-g40464-gd386-xx32824-4016.circ.txt) | 7640 | <ins><strong>70</strong></ins> | 40464 | 386 | 32824 | [A47/AD3/G225/GD15](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [A36/AD4/G147/GD14](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ](aes128-decipher-a7560-ad80-g37504-gd356-xx29944-4896.circ.txt) | <ins><strong>7560</strong></ins> | 80 | <ins><strong>37504</strong></ins> | <ins><strong>356</strong></ins> | 29944 | [A45/AD4/G151/GD12](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [A36/AD4/G147/GD14](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 1,3,4 |

## AES-192

### Decipher using [Inv]Sbox with A<=34

| File | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes192-decipher-a6496-ad100-g40440-gd804-xx33944-6216.circ.txt) | <ins><strong>6496</strong></ins> | 100 | 40440 | 804 | 33944 | [A29/AD5/G139/GD35](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [A29/AD5/G145/GD32](../aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [G114/GD8](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](aes192-decipher-a7616-ad80-g39384-gd414-xx31768-3592.circ.txt) | 7616 | <ins><strong>80</strong></ins> | 39384 | <ins><strong>414</strong></ins> | 31768 | [A34/AD4/G128/GD15](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [A34/AD4/G134/GD15](../aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [G146/GD5](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](aes192-decipher-a7168-ad100-g33176-gd655-xx26008-1832.circ.txt) | 7168 | 100 | <ins><strong>33176</strong></ins> | 655 | 26008 | [A32/AD5/G110/GD23](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [A32/AD5/G112/GD27](../aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [G114/GD8](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

### Decipher using [Inv]Sbox with A>34

| File | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes192-decipher-a8416-ad72-g44984-gd402-xx36568-4744.circ.txt) | 8416 | <ins><strong>72</strong></ins> | 44984 | 402 | 36568 | [A47/AD3/G225/GD15](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [A36/AD4/G147/GD14](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ](aes192-decipher-a8352-ad80-g42616-gd378-xx34264-5448.circ.txt) | <ins><strong>8352</strong></ins> | 80 | <ins><strong>42616</strong></ins> | <ins><strong>378</strong></ins> | 34264 | [A45/AD4/G151/GD12](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [A36/AD4/G147/GD14](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 1,3,4 |

## AES-256

### Decipher using [Inv]Sbox with A<=34

| File | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes256-decipher-a8004-ad135-g49220-gd1062-xx41216-7543.circ.txt) | <ins><strong>8004</strong></ins> | 135 | 49220 | 1062 | 41216 | [A29/AD5/G139/GD35](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [A29/AD5/G145/GD32](../aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [G114/GD8](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](aes256-decipher-a9384-ad108-g47848-gd537-xx38464-4247.circ.txt) | 9384 | <ins><strong>108</strong></ins> | 47848 | <ins><strong>537</strong></ins> | 38464 | [A34/AD4/G128/GD15](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [A34/AD4/G134/GD15](../aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [G146/GD5](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](aes256-decipher-a8832-ad135-g40320-gd848-xx31488-2179.circ.txt) | 8832 | 135 | <ins><strong>40320</strong></ins> | 848 | 31488 | [A32/AD5/G110/GD23](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [A32/AD5/G112/GD27](../aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [G114/GD8](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

### Decipher using [Inv]Sbox with A>34

| File | A | AD | G | GD | XX | Sbox | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---|---|---|:---:|
| [circ](aes256-decipher-a10508-ad95-g55804-gd523-xx45296-5591.circ.txt) | 10508 | <ins><strong>95</strong></ins> | 55804 | 523 | 45296 | [A47/AD3/G225/GD15](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [A36/AD4/G147/GD14](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2 |
| [circ](aes256-decipher-a10404-ad108-g51956-gd484-xx41552-6735.circ.txt) | <ins><strong>10404</strong></ins> | 108 | <ins><strong>51956</strong></ins> | <ins><strong>484</strong></ins> | 41552 | [A45/AD4/G151/GD12](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [A36/AD4/G147/GD14](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [G146/GD5](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 1,3,4 |
