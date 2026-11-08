# AES S-box Sel (with direction *sel*ection)

The aes-sbox-sel is a 9-bit to 8-bit vectorial Boolean function whose 1st bit (c) selects whether the other 8 input bits are transformed with [aes-sbox](../aes-sbox/) or with [aes-invsbox](../aes-invsbox/).

The following two tables list selected circuits from two categories: A36 and 40≤A≤41

## S-box Sel with #AND = 36

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-sel-a36-ad6-g187-gd27-xx151-42.circ.txt) | 36 | 6 | <ins><strong>187</strong></ins> | 27 | 151 | 109 | 42 |
| [circ](aes-sbox-sel-a36-ad6-g191-gd26-xx155-31.circ.txt) | 36 | 6 | 191 | <ins><strong>26</strong></ins> | 155 | 124 | 31 |

The listed aes-sbox-sel/A36 circuits were obtained from transformations to sbox/A28 circuits.

## S-box Sel with 40 ≤ #AND ≤ 41

| File | #AND<br>(A) | AND<br>depth | #Gates<br>(G) | Gate<br>depth | XX | #XOR | #XNOR |
|---|---:|---:|---:|---:|---:|---:|---:|
| [circ](aes-sbox-sel-a40-ad6-g184-gd23-xx144-10.circ.txt) | <ins><strong>40</strong></ins> | <ins><strong>6</strong></ins> | 184 | 23 | 144 | 134 | 10 |
| [circ](aes-sbox-sel-a40-ad6-g191-gd20-xx151-11.circ.txt) | <ins><strong>40</strong></ins> | <ins><strong>6</strong></ins> | 191 | <ins><strong>20</strong></ins> | 151 | 140 | 11 |
| [circ](aes-sbox-sel-a40-ad7-g153-gd28-xx113-6.circ.txt) | <ins><strong>40</strong></ins> | 7 | 153 | 28 | 113 | 107 | 6 |
| [circ](aes-sbox-sel-a40-ad7-g166-gd25-xx126-23.circ.txt) | <ins><strong>40</strong></ins> | 7 | 166 | 25 | 126 | 103 | 23 |
| [circ](aes-sbox-sel-a41-ad6-g151-gd24-xx110-10.circ.txt) | 41 | <ins><strong>6</strong></ins> | <ins><strong>151</strong></ins> | 24 | 110 | 100 | 10 |
| [circ](aes-sbox-sel-a41-ad6-g165-gd21-xx124-10.circ.txt) | 41 | <ins><strong>6</strong></ins> | 165 | 21 | 124 | 114 | 10 |
| [circ](aes-sbox-sel-a41-ad6-g173-gd20-xx132-10.circ.txt) | 41 | <ins><strong>6</strong></ins> | 173 | <ins><strong>20</strong></ins> | 132 | 122 | 10 |

The listed aes-sbox-sel/A{40,41} circuits were obtained from transformations to sbox/A{32,33} circuits.
