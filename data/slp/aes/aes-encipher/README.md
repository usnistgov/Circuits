# AES Encipher (flat circuits)

TM = Tuple metric. TM1 = A-AD-G-GD; TM2 = AD-GD-G-A; TM3 = GD-G-AD-A; TM4 = G-A-GD-AD.

A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR.

Within each table, bold underlined values mark the lowest displayed value in selected columns.

## AES-128

### Encipher using Sbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes128-encipher-a5800-ad50-g33656-gd388-xx27856-4016.circ.txt) | <ins><strong>5800</strong></ins> | 50 | 33656 | 388 | 27856 | 23840 | 4016 | [circ](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](aes128-encipher-a6800-ad40-g31780-gd191-xx24980-816.circ.txt) | 6800 | <ins><strong>40</strong></ins> | 31780 | <ins><strong>191</strong></ins> | 24980 | 24164 | 816 | [circ](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](aes128-encipher-a6400-ad50-g27856-gd286-xx21456-616.circ.txt) | 6400 | 50 | <ins><strong>27856</strong></ins> | 286 | 21456 | 20840 | 616 | [circ](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

### Encipher using Sbox with A>34

| File | A | AD | G | GD | XX | X | X' | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes128-encipher-a7200-ad40-g33456-gd196-xx26256-4416.circ.txt) | <ins><strong>7200</strong></ins> | 40 | <ins><strong>33456</strong></ins> | 196 | 26256 | 21840 | 4416 | [circ](../aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1,4 |
| [circ](aes128-encipher-a9400-ad30-g51180-gd191-xx41780-816.circ.txt) | 9400 | <ins><strong>30</strong></ins> | 51180 | 191 | 41780 | 40964 | 816 | [circ](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](aes128-encipher-a9000-ad40-g36380-gd161-xx27380-5216.circ.txt) | 9000 | 40 | 36380 | <ins><strong>161</strong></ins> | 27380 | 22164 | 5216 | [circ](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 3 |

## AES-192

### Encipher using Sbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes192-encipher-a6496-ad60-g38144-gd466-xx31648-4488.circ.txt) | <ins><strong>6496</strong></ins> | 60 | 38144 | 466 | 31648 | 27160 | 4488 | [circ](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](aes192-encipher-a7616-ad48-g36076-gd226-xx28460-904.circ.txt) | 7616 | <ins><strong>48</strong></ins> | 36076 | <ins><strong>226</strong></ins> | 28460 | 27556 | 904 | [circ](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](aes192-encipher-a7168-ad60-g31648-gd344-xx24480-680.circ.txt) | 7168 | 60 | <ins><strong>31648</strong></ins> | 344 | 24480 | 23800 | 680 | [circ](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

### Encipher using Sbox with A>34

| File | A | AD | G | GD | XX | X | X' | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes192-encipher-a8064-ad48-g37920-gd236-xx29856-4936.circ.txt) | <ins><strong>8064</strong></ins> | 48 | <ins><strong>37920</strong></ins> | 236 | 29856 | 24920 | 4936 | [circ](../aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1,4 |
| [circ](aes192-encipher-a10528-ad36-g57804-gd226-xx47276-904.circ.txt) | 10528 | <ins><strong>36</strong></ins> | 57804 | 226 | 47276 | 46372 | 904 | [circ](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](aes192-encipher-a10080-ad48-g41228-gd190-xx31148-5832.circ.txt) | 10080 | 48 | 41228 | <ins><strong>190</strong></ins> | 31148 | 25316 | 5832 | [circ](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 3 |

## AES-256

### Encipher using Sbox with A<=34

| File | A | AD | G | GD | XX | X | X' | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes256-encipher-a8004-ad70-g46524-gd544-xx38520-5527.circ.txt) | <ins><strong>8004</strong></ins> | 70 | 46524 | 544 | 38520 | 32993 | 5527 | [circ](../aes-sbox/aes-sbox-a29-ad5-g139-gd35-xx110-20.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1 |
| [circ](aes256-encipher-a9384-ad56-g43956-gd264-xx34572-1111.circ.txt) | 9384 | <ins><strong>56</strong></ins> | 43956 | <ins><strong>264</strong></ins> | 34572 | 33461 | 1111 | [circ](../aes-sbox/aes-sbox-a34-ad4-g128-gd15-xx94-4.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2,3 |
| [circ](aes256-encipher-a8832-ad70-g38520-gd402-xx29688-835.circ.txt) | 8832 | 70 | <ins><strong>38520</strong></ins> | 402 | 29688 | 28853 | 835 | [circ](../aes-sbox/aes-sbox-a32-ad5-g110-gd23-xx78-3.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 4 |

### Encipher using Sbox with A>34

| File | A | AD | G | GD | XX | X | X' | Sbox | Mix<br>Cols | TM |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|:---:|
| [circ](aes256-encipher-a9936-ad56-g46248-gd276-xx36312-6079.circ.txt) | <ins><strong>9936</strong></ins> | 56 | <ins><strong>46248</strong></ins> | 276 | 36312 | 30233 | 6079 | [circ](../aes-sbox/aes-sbox-a36-ad4-g138-gd14-xx102-22.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor88-depth5.circ.txt) | 1,4 |
| [circ](aes256-encipher-a12972-ad42-g70728-gd264-xx57756-1111.circ.txt) | 12972 | <ins><strong>42</strong></ins> | 70728 | 264 | 57756 | 56645 | 1111 | [circ](../aes-sbox/aes-sbox-a47-ad3-g225-gd15-xx178-4.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 2 |
| [circ](aes256-encipher-a12420-ad56-g50304-gd222-xx37884-7183.circ.txt) | 12420 | 56 | 50304 | <ins><strong>222</strong></ins> | 37884 | 30701 | 7183 | [circ](../aes-sbox/aes-sbox-a45-ad4-g151-gd12-xx106-26.circ.txt) | [circ](../aes-mixcols/aes-mixcols-xor97-depth3.circ.txt) | 3 |
