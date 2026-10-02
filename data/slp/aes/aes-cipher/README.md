# AES Cipher (flat circuits)

## Notation:
- A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR.
- TM = Tuple metric. TM1 = A-AD-G-GD; TM2 = AD-GD-G-A; TM3 = GD-G-AD-A; TM4 = G-A-GD-AD.
- Within each table, bold underlined values mark the lowest displayed value in selected columns.

## AES-128

### Cipher using Sbox with A<=34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes128-cipher-a4640-ad50-g26816-gd388-xx22176-3200.circ.txt) | <ins><strong>4640</strong></ins> | 50 | 26816 | 388 | 22176 | [A29/AD5/G139/GD35](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](aes128-cipher-a5440-ad40-g25380-gd188-xx19940-640.circ.txt) | 5440 | <ins><strong>40</strong></ins> | 25380 | <ins><strong>188</strong></ins> | 19940 | [A34/AD4/G128/GD15](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](aes128-cipher-a5120-ad50-g22176-gd286-xx17056-480.circ.txt) | 5120 | 50 | <ins><strong>22176</strong></ins> | 286 | 17056 | [A32/AD5/G110/GD23](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

### Cipher using Sbox with A>34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes128-cipher-a7520-ad30-g40900-gd188-xx33380-640.circ.txt) | 7520 | <ins><strong>30</strong></ins> | 40900 | 188 | 33380 | [A47/AD3/G225/GD15](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](aes128-cipher-a7200-ad40-g29060-gd158-xx21860-4160.circ.txt) | <ins><strong>7200</strong></ins> | 40 | <ins><strong>29060</strong></ins> | <ins><strong>158</strong></ins> | 21860 | [A45/AD4/G151/GD12](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |

## AES-192

### Cipher using Sbox with A<=34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes192-cipher-a5568-ad60-g32224-gd466-xx26656-3840.circ.txt) | <ins><strong>5568</strong></ins> | 60 | 32224 | 466 | 26656 | [A29/AD5/G139/GD35](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](aes192-cipher-a6528-ad48-g30508-gd226-xx23980-768.circ.txt) | 6528 | <ins><strong>48</strong></ins> | 30508 | <ins><strong>226</strong></ins> | 23980 | [A34/AD4/G128/GD15](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](aes192-cipher-a6144-ad60-g26656-gd344-xx20512-576.circ.txt) | 6144 | 60 | <ins><strong>26656</strong></ins> | 344 | 20512 | [A32/AD5/G110/GD23](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

### Cipher using Sbox with A>34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes192-cipher-a9024-ad36-g49132-gd226-xx40108-768.circ.txt) | 9024 | <ins><strong>36</strong></ins> | 49132 | 226 | 40108 | [A47/AD3/G225/GD15](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](aes192-cipher-a8640-ad48-g34924-gd190-xx26284-4992.circ.txt) | <ins><strong>8640</strong></ins> | 48 | <ins><strong>34924</strong></ins> | <ins><strong>190</strong></ins> | 26284 | [A45/AD4/G151/GD12](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |

## AES-256

### Cipher using Sbox with A<=34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes256-cipher-a6496-ad70-g37632-gd544-xx31136-4480.circ.txt) | <ins><strong>6496</strong></ins> | 70 | 37632 | 544 | 31136 | [A29/AD5/G139/GD35](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](aes256-cipher-a7616-ad56-g35636-gd264-xx28020-896.circ.txt) | 7616 | <ins><strong>56</strong></ins> | 35636 | <ins><strong>264</strong></ins> | 28020 | [A34/AD4/G128/GD15](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](aes256-cipher-a7168-ad70-g31136-gd402-xx23968-672.circ.txt) | 7168 | 70 | <ins><strong>31136</strong></ins> | 402 | 23968 | [A32/AD5/G110/GD23](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

### Cipher using Sbox with A>34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes256-cipher-a10528-ad42-g57364-gd264-xx46836-896.circ.txt) | 10528 | <ins><strong>42</strong></ins> | 57364 | 264 | 46836 | [A47/AD3/G225/GD15](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](aes256-cipher-a10080-ad56-g40788-gd222-xx30708-5824.circ.txt) | <ins><strong>10080</strong></ins> | 56 | <ins><strong>40788</strong></ins> | <ins><strong>222</strong></ins> | 30708 | [A45/AD4/G151/GD12](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |
