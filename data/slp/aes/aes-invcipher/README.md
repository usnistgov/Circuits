# AES InvCipher (flat circuits)

TM = Tuple metric. TM1 = A-AD-G-GD; TM2 = AD-GD-G-A; TM3 = GD-G-AD-A; TM4 = G-A-GD-AD.

A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR.

Within each table, bold underlined values mark the lowest displayed value in selected columns.

## AES-128

### InvCipher using InvSbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes128-invcipher-a4640-ad50-g28712-gd403-xx24072-4640.circ.txt) | <ins><strong>4640</strong></ins> | 50 | 28712 | 403 | 24072 | 19432 | 4640 | [circ](../aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](aes128-invcipher-a5440-ad40-g28104-gd206-xx22664-2880.circ.txt) | 5440 | <ins><strong>40</strong></ins> | 28104 | <ins><strong>206</strong></ins> | 22664 | 19784 | 2880 | [circ](../aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](aes128-invcipher-a5120-ad50-g23432-gd353-xx18312-1440.circ.txt) | 5120 | 50 | <ins><strong>23432</strong></ins> | 353 | 18312 | 16872 | 1440 | [circ](../aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

### InvCipher using InvSbox with A>34

| File | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes128-invcipher-a5760-ad40-g29032-gd223-xx23272-3840.circ.txt) | <ins><strong>5760</strong></ins> | <ins><strong>40</strong></ins> | <ins><strong>29032</strong></ins> | 223 | 23272 | 19432 | 3840 | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ](aes128-invcipher-a5760-ad40-g30184-gd196-xx24424-3840.circ.txt) | <ins><strong>5760</strong></ins> | <ins><strong>40</strong></ins> | 30184 | <ins><strong>196</strong></ins> | 24424 | 20584 | 3840 | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |

## AES-192

### InvCipher using InvSbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes192-invcipher-a5568-ad60-g34520-gd485-xx28952-5568.circ.txt) | <ins><strong>5568</strong></ins> | 60 | 34520 | 485 | 28952 | 23384 | 5568 | [circ](../aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](aes192-invcipher-a6528-ad48-g33816-gd248-xx27288-3456.circ.txt) | 6528 | <ins><strong>48</strong></ins> | 33816 | <ins><strong>248</strong></ins> | 27288 | 23832 | 3456 | [circ](../aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](aes192-invcipher-a6144-ad60-g28184-gd425-xx22040-1728.circ.txt) | 6144 | 60 | <ins><strong>28184</strong></ins> | 425 | 22040 | 20312 | 1728 | [circ](../aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

### InvCipher using InvSbox with A>34

| File | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes192-invcipher-a6912-ad48-g34904-gd269-xx27992-4608.circ.txt) | <ins><strong>6912</strong></ins> | <ins><strong>48</strong></ins> | <ins><strong>34904</strong></ins> | 269 | 27992 | 23384 | 4608 | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ](aes192-invcipher-a6912-ad48-g36312-gd236-xx29400-4608.circ.txt) | <ins><strong>6912</strong></ins> | <ins><strong>48</strong></ins> | 36312 | <ins><strong>236</strong></ins> | 29400 | 24792 | 4608 | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |

## AES-256

### InvCipher using InvSbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes256-invcipher-a6496-ad70-g40328-gd567-xx33832-6496.circ.txt) | <ins><strong>6496</strong></ins> | 70 | 40328 | 567 | 33832 | 27336 | 6496 | [circ](../aes-invsbox/aes-invsbox-a29-ad5-g145-gd32-xx116-29.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1 |
| [circ](aes256-invcipher-a7616-ad56-g39528-gd290-xx31912-4032.circ.txt) | 7616 | <ins><strong>56</strong></ins> | 39528 | <ins><strong>290</strong></ins> | 31912 | 27880 | 4032 | [circ](../aes-invsbox/aes-invsbox-a34-ad4-g134-gd15-xx100-18.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
| [circ](aes256-invcipher-a7168-ad70-g32936-gd497-xx25768-2016.circ.txt) | 7168 | 70 | <ins><strong>32936</strong></ins> | 497 | 25768 | 23752 | 2016 | [circ](../aes-invsbox/aes-invsbox-a32-ad5-g112-gd27-xx80-9.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 4 |

### InvCipher using InvSbox with A>34

| File | A | AD | G | GD | XX | X | X' | Inv<br>Sbox | Inv<br>MixCols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes256-invcipher-a8064-ad56-g40776-gd315-xx32712-5376.circ.txt) | <ins><strong>8064</strong></ins> | <ins><strong>56</strong></ins> | <ins><strong>40776</strong></ins> | 315 | 32712 | 27336 | 5376 | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor114-depth8.circ.txt) | 1,4 |
| [circ](aes256-invcipher-a8064-ad56-g42440-gd276-xx34376-5376.circ.txt) | <ins><strong>8064</strong></ins> | <ins><strong>56</strong></ins> | 42440 | <ins><strong>276</strong></ins> | 34376 | 29000 | 5376 | [circ](../aes-invsbox/aes-invsbox-a36-ad4-g147-gd14-xx111-24.circ.txt) | [circ](../aes-invmixcols/aes-invmixcols-xor146-depth5.circ.txt) | 2,3 |
