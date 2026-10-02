# AES Encipher (flat circuits)

TM = Tuple metric. TM1 = A-AD-G-GD; TM2 = AD-GD-G-A; TM3 = GD-G-AD-A; TM4 = G-A-GD-AD.

A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = #XOR + #XNOR.

Within each table, bold underlined values mark the lowest displayed value in selected columns.

## AES-128

### Encipher using Sbox with A<=34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes128-encipher-a5800-ad50-g33656-gd388-xx27856-4016.circ.txt) | <ins><strong>5800</strong></ins> | 50 | 33656 | 388 | 27856 | [A29/AD5/G139/GD35](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](aes128-encipher-a6800-ad40-g31780-gd191-xx24980-816.circ.txt) | 6800 | <ins><strong>40</strong></ins> | 31780 | <ins><strong>191</strong></ins> | 24980 | [A34/AD4/G128/GD15](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](aes128-encipher-a6400-ad50-g27856-gd286-xx21456-616.circ.txt) | 6400 | 50 | <ins><strong>27856</strong></ins> | 286 | 21456 | [A32/AD5/G110/GD23](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

### Encipher using Sbox with A>34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes128-encipher-a9400-ad30-g51180-gd191-xx41780-816.circ.txt) | 9400 | <ins><strong>30</strong></ins> | 51180 | 191 | 41780 | [A47/AD3/G225/GD15](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](aes128-encipher-a9000-ad40-g36380-gd161-xx27380-5216.circ.txt) | <ins><strong>9000</strong></ins> | 40 | <ins><strong>36380</strong></ins> | <ins><strong>161</strong></ins> | 27380 | [A45/AD4/G151/GD12](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |

## AES-192

### Encipher using Sbox with A<=34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes192-encipher-a6496-ad60-g38144-gd466-xx31648-4488.circ.txt) | <ins><strong>6496</strong></ins> | 60 | 38144 | 466 | 31648 | [A29/AD5/G139/GD35](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](aes192-encipher-a7616-ad48-g36076-gd226-xx28460-904.circ.txt) | 7616 | <ins><strong>48</strong></ins> | 36076 | <ins><strong>226</strong></ins> | 28460 | [A34/AD4/G128/GD15](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](aes192-encipher-a7168-ad60-g31648-gd344-xx24480-680.circ.txt) | 7168 | 60 | <ins><strong>31648</strong></ins> | 344 | 24480 | [A32/AD5/G110/GD23](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

### Encipher using Sbox with A>34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes192-encipher-a10528-ad36-g57804-gd226-xx47276-904.circ.txt) | 10528 | <ins><strong>36</strong></ins> | 57804 | 226 | 47276 | [A47/AD3/G225/GD15](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](aes192-encipher-a10080-ad48-g41228-gd190-xx31148-5832.circ.txt) | <ins><strong>10080</strong></ins> | 48 | <ins><strong>41228</strong></ins> | <ins><strong>190</strong></ins> | 31148 | [A45/AD4/G151/GD12](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |

## AES-256

### Encipher using Sbox with A<=34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes256-encipher-a8004-ad70-g46524-gd544-xx38520-5527.circ.txt) | <ins><strong>8004</strong></ins> | 70 | 46524 | 544 | 38520 | [A29/AD5/G139/GD35](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](aes256-encipher-a9384-ad56-g43956-gd264-xx34572-1111.circ.txt) | 9384 | <ins><strong>56</strong></ins> | 43956 | <ins><strong>264</strong></ins> | 34572 | [A34/AD4/G128/GD15](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](aes256-encipher-a8832-ad70-g38520-gd402-xx29688-835.circ.txt) | 8832 | 70 | <ins><strong>38520</strong></ins> | 402 | 29688 | [A32/AD5/G110/GD23](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [G88/GD5](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

### Encipher using Sbox with A>34

| File | A | AD | G | GD | XX | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes256-encipher-a12972-ad42-g70728-gd264-xx57756-1111.circ.txt) | 12972 | <ins><strong>42</strong></ins> | 70728 | 264 | 57756 | [A47/AD3/G225/GD15](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](aes256-encipher-a12420-ad56-g50304-gd222-xx37884-7183.circ.txt) | <ins><strong>12420</strong></ins> | 56 | <ins><strong>50304</strong></ins> | <ins><strong>222</strong></ins> | 37884 | [A45/AD4/G151/GD12](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [G97/GD3](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 1,3,4 |
