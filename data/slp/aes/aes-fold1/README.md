# AES "folded" circuits

In this page a _folded_ circuit (as opposed to a _flat_ circuit) is a straight-line program (SLP) that uses calls to subcircuits (e.g., `aes-[inv]sbox`, `aes-[inv]mixcols`) as opposed to just basic Boolean gates.

`VXOR` and `VXNOR` are succinct notation for parallel execution of `XOR` and `XNOR` gates, respectively.

The circuit basenames relate to [FIPS 197-upd1](https://csrc.nist.gov/pubs/fips/197/final) functions, but with some differences:

- `keyexp`: Implements `KeyExpansion()`, but with an output that does not include the original key;
- `cipher`: Implements `Cipher()`, whose input includes an expanded key;
- `encipher`: Implements `AES-128/192/256()`, using `KeyExpansion()` in parallel with `Cipher()`;
- `invcipher`: Implements `InvCipher()`, whose input includes an expanded key;
- `decipher`: Not defined in FIPS 197, but composes `KeyExpansion()` with `InvCipher()`.

## Folded AES Circuits

| Circuit | FIPS 197<br>operation | Key<br>size | aes-<br>sbox | aes-<br>invsbox | aes-<br>mixcols | aes-inv<br>mixcols | VXOR<br>(XORs) | VXNOR<br>(XNORs) | XNOR |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| [aes128-keyexp](aes128-keyexp-fold1.circ.txt) | KeyExpansion() | 128 | 40 | — | — | — | 40 (1264) | 2 (8) | 8 |
| [aes192-keyexp](aes192-keyexp-fold1.circ.txt) | KeyExpansion() | 192 | 32 | — | — | — | 46 (1464) | — | 8 |
| [aes256-keyexp](aes256-keyexp-fold1.circ.txt) | KeyExpansion() | 256 | 52 | — | — | — | 52 (1657) | — | 7 |
| [aes128-cipher](aes128-cipher-fold1.circ.txt) | Cipher() | 128 | 160 | — | 36 | — | 11 (1408) | — | — |
| [aes192-cipher](aes192-cipher-fold1.circ.txt) | Cipher() | 192 | 192 | — | 44 | — | 14 (1664) | — | — |
| [aes256-cipher](aes256-cipher-fold1.circ.txt) | Cipher() | 256 | 224 | — | 52 | — | 15 (1920) | — | — |
| [aes128-encipher](aes128-encipher-fold1.circ.txt) | AES-128() | 128 | 200 | — | 36 | — | 51 (2672) | 2 (8) | 8 |
| [aes192-encipher](aes192-encipher-fold1.circ.txt) | AES-192() | 192 | 224 | — | 44 | — | 60 (3128) | — | 8 |
| [aes256-encipher](aes256-encipher-fold1.circ.txt) | AES-256() | 256 | 276 | — | 52 | — | 67 (3577) | — | 7 |
| [aes128-invcipher](aes128-invcipher-fold1.circ.txt) | InvCipher() | 128 | — | 160 | — | 36 | 11 (1408) | — | — |
| [aes192-invcipher](aes192-invcipher-fold1.circ.txt) | InvCipher() | 192 | — | 192 | — | 44 | 14 (1664) | — | — |
| [aes256-invcipher](aes256-invcipher-fold1.circ.txt) | InvCipher() | 256 | — | 224 | — | 52 | 15 (1920) | — | — |
| [aes128-decipher](aes128-decipher-fold1.circ.txt) | — | 128 | 40 | 160 | — | 36 | 51 (2672) | 2 (8) | 8 |
| [aes192-decipher](aes192-decipher-fold1.circ.txt) | — | 192 | 32 | 192 | — | 44 | 60 (3128) | — | 8 |
| [aes256-decipher](aes256-decipher-fold1.circ.txt) | — | 256 | 52 | 224 | — | 52 | 67 (3577) | — | 7 |

Note: In FIPS 197, the input `key` in `KeyExpansion()` also appears in the output; whereas it does not in our `keyexp` circuit.
