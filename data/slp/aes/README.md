<!-- CreatedUTC: 2026-09-02T16:37:55Z -->

# AES Circuits

Examples of Boolean circuits (straight-line programs) related to the *Advanced Encryption Standard* (AES).

<details open>
<summary><h2>AES S-box</h2></summary>

Example circuits for the AES S-box (8-bit to 8-bit function), for both Forward and Inverse directions.

Note: In selected columns, <sub><sup>🟩</sup></sub> <ins><strong>indicates</strong></ins> the lowest displayed value.

<details open>
<summary><h3>S-box Forward with #AND ≤ 34</h3></summary>

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">#AND<br>(<em>A</em>)</th>
      <th align="center">AND<br>depth</th>
      <th align="center">#Gates<br>(<em>G</em>)</th>
      <th align="center">Gate<br>depth</th>
      <th align="center">XX<br>(<em>x</em>+<em>x'</em>)</th>
      <th align="center">#XOR<br>(<em>x</em>)</th>
      <th align="center">#XNOR<br>(<em>x'</em>)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a29-ad5-g139-gd35-xx110-20.circ.txt">circ.txt</a></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29</strong></ins></td>
      <td align="right">5</td>
      <td align="right">139</td>
      <td align="right">35</td>
      <td align="right">110</td>
      <td align="right">90</td>
      <td align="right">20</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a29-ad5-g184-gd16-xx155-26.circ.txt">circ.txt</a></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29</strong></ins></td>
      <td align="right">5</td>
      <td align="right">184</td>
      <td align="right">16</td>
      <td align="right">155</td>
      <td align="right">129</td>
      <td align="right">26</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a29-ad6-g138-gd38-xx109-14.circ.txt">circ.txt</a></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29</strong></ins></td>
      <td align="right">6</td>
      <td align="right">138</td>
      <td align="right">38</td>
      <td align="right">109</td>
      <td align="right">95</td>
      <td align="right">14</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a30-ad4-g147-gd36-xx117-13.circ.txt">circ.txt</a></td>
      <td align="right">30</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>4</strong></ins></td>
      <td align="right">147</td>
      <td align="right">36</td>
      <td align="right">117</td>
      <td align="right">104</td>
      <td align="right">13</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a30-ad4-g205-gd15-xx175-54.circ.txt">circ.txt</a></td>
      <td align="right">30</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>4</strong></ins></td>
      <td align="right">205</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>15</strong></ins></td>
      <td align="right">175</td>
      <td align="right">121</td>
      <td align="right">54</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g110-gd23-xx76-3.circ.txt">circ.txt</a></td>
      <td align="right">34</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>4</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>110</strong></ins></td>
      <td align="right">23</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>76</strong></ins></td>
      <td align="right">73</td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">circ.txt</a></td>
      <td align="right">34</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>4</strong></ins></td>
      <td align="right">128</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>15</strong></ins></td>
      <td align="right">94</td>
      <td align="right">90</td>
      <td align="right">4</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a32-ad5-g110-gd23-xx78-3.circ.txt">circ.txt</a></td>
      <td align="right">32</td>
      <td align="right">5</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>110</strong></ins></td>
      <td align="right">23</td>
      <td align="right">78</td>
      <td align="right">75</td>
      <td align="right">3</td>
    </tr>
  </tbody>
</table>


</details>

<details open>
<summary><h3>S-box Inverse with #AND ≤ 34</h3></summary>

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">#AND<br>(<em>A</em>)</th>
      <th align="center">AND<br>depth</th>
      <th align="center">#Gates<br>(<em>G</em>)</th>
      <th align="center">Gate<br>depth</th>
      <th align="center">XX<br>(<em>x</em>+<em>x'</em>)</th>
      <th align="center">#XOR<br>(<em>x</em>)</th>
      <th align="center">#XNOR<br>(<em>x'</em>)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./sbox/aes-sbox-inv-a29-ad5-g145-gd32-xx116-29.circ.txt">circ.txt</a></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29</strong></ins></td>
      <td align="right">5</td>
      <td align="right">145</td>
      <td align="right">32</td>
      <td align="right">116</td>
      <td align="right">87</td>
      <td align="right">29</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-inv-a29-ad5-g186-gd16-xx157-43.circ.txt">circ.txt</a></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29</strong></ins></td>
      <td align="right">5</td>
      <td align="right">186</td>
      <td align="right">16</td>
      <td align="right">157</td>
      <td align="right">114</td>
      <td align="right">43</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-inv-a29-ad6-g144-gd34-xx115-21.circ.txt">circ.txt</a></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29</strong></ins></td>
      <td align="right">6</td>
      <td align="right">144</td>
      <td align="right">34</td>
      <td align="right">115</td>
      <td align="right">94</td>
      <td align="right">21</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-inv-a30-ad4-g150-gd30-xx120-33.circ.txt">circ.txt</a></td>
      <td align="right">30</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>4</strong></ins></td>
      <td align="right">150</td>
      <td align="right">30</td>
      <td align="right">120</td>
      <td align="right">87</td>
      <td align="right">33</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-inv-a30-ad4-g203-gd15-xx173-52.circ.txt">circ.txt</a></td>
      <td align="right">30</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>4</strong></ins></td>
      <td align="right">203</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>15</strong></ins></td>
      <td align="right">173</td>
      <td align="right">121</td>
      <td align="right">52</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-inv-a34-ad4-g114-gd27-xx80-11.circ.txt">circ.txt</a></td>
      <td align="right">34</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>4</strong></ins></td>
      <td align="right">114</td>
      <td align="right">27</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>80</strong></ins></td>
      <td align="right">69</td>
      <td align="right">11</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-inv-a34-ad4-g134-gd15-xx100-18.circ.txt">circ.txt</a></td>
      <td align="right">34</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>4</strong></ins></td>
      <td align="right">134</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>15</strong></ins></td>
      <td align="right">100</td>
      <td align="right">82</td>
      <td align="right">18</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-inv-a32-ad5-g112-gd27-xx80-9.circ.txt">circ.txt</a></td>
      <td align="right">32</td>
      <td align="right">5</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>112</strong></ins></td>
      <td align="right">27</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>80</strong></ins></td>
      <td align="right">71</td>
      <td align="right">9</td>
    </tr>
  </tbody>
</table>


</details>

<details open>
<summary><h3>S-box Forward with #AND &gt; 34</h3></summary>

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">#AND<br>(<em>A</em>)</th>
      <th align="center">AND<br>depth</th>
      <th align="center">#Gates<br>(<em>G</em>)</th>
      <th align="center">Gate<br>depth</th>
      <th align="center">XX<br>(<em>x</em>+<em>x'</em>)</th>
      <th align="center">#XOR<br>(<em>x</em>)</th>
      <th align="center">#XNOR<br>(<em>x'</em>)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a46-ad3-g250-gd15-xx204-4.circ.txt">circ.txt</a></td>
      <td align="right">46</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>3</strong></ins></td>
      <td align="right">250</td>
      <td align="right">15</td>
      <td align="right">204</td>
      <td align="right">200</td>
      <td align="right">4</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a47-ad3-g225-gd15-xx178-4.circ.txt">circ.txt</a></td>
      <td align="right">47</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>3</strong></ins></td>
      <td align="right">225</td>
      <td align="right">15</td>
      <td align="right">178</td>
      <td align="right">174</td>
      <td align="right">4</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">circ.txt</a></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>36</strong></ins></td>
      <td align="right">4</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>138</strong></ins></td>
      <td align="right">14</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>102</strong></ins></td>
      <td align="right">80</td>
      <td align="right">22</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a39-ad4-g148-gd13-xx109-24.circ.txt">circ.txt</a></td>
      <td align="right">39</td>
      <td align="right">4</td>
      <td align="right">148</td>
      <td align="right">13</td>
      <td align="right">109</td>
      <td align="right">85</td>
      <td align="right">24</td>
    </tr>
    <tr>
      <td><a href="./sbox/aes-sbox-fwd-a45-ad4-g151-gd12-xx106-26.circ.txt">circ.txt</a></td>
      <td align="right">45</td>
      <td align="right">4</td>
      <td align="right">151</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>12</strong></ins></td>
      <td align="right">106</td>
      <td align="right">80</td>
      <td align="right">26</td>
    </tr>
  </tbody>
</table>


</details>

<details open>
<summary><h3>S-box Inverse with #AND &gt; 34</h3></summary>

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">#AND<br>(<em>A</em>)</th>
      <th align="center">AND<br>depth</th>
      <th align="center">#Gates<br>(<em>G</em>)</th>
      <th align="center">Gate<br>depth</th>
      <th align="center">XX<br>(<em>x</em>+<em>x'</em>)</th>
      <th align="center">#XOR<br>(<em>x</em>)</th>
      <th align="center">#XNOR<br>(<em>x'</em>)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">circ.txt</a></td>
      <td align="right">36</td>
      <td align="right">4</td>
      <td align="right">148</td>
      <td align="right">14</td>
      <td align="right">112</td>
      <td align="right">88</td>
      <td align="right">24</td>
    </tr>
  </tbody>
</table>


</details>

<details>
<summary><h3>S-box circuits added around 2020</h3></summary>

The following historical table is reproduced from the <a href="https://csrc.nist.gov/Projects/circuit-complexity/list-of-circuits">NIST Circuit Complexity list of circuits</a>.

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">Direction</th>
      <th align="center">#AND<br>(<em>A</em>)</th>
      <th align="center">#Gates<br>(<em>G</em>)</th>
      <th align="center">Gate<br>depth</th>
      <th align="center">XX<br>(<em>x</em>+<em>x'</em>)</th>
      <th align="center">#XOR<br>(<em>x</em>)</th>
      <th align="center">#XNOR<br>(<em>x'</em>)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./aes-sbox-fwd-g115-a32-d28-ad6.slp">slp</a></td>
      <td>Forward</td>
      <td align="right">32</td>
      <td align="right">115</td>
      <td align="right">28</td>
      <td align="right">83</td>
      <td align="right">79</td>
      <td align="right">4</td>
    </tr>
    <tr>
      <td><a href="./aes-sbox-fwd-g113-a32-d27-ad6.slp">slp</a></td>
      <td>Forward</td>
      <td align="right">32</td>
      <td align="right">113</td>
      <td align="right">27</td>
      <td align="right">81</td>
      <td align="right">77</td>
      <td align="right">4</td>
    </tr>
    <tr>
      <td><a href="./aes-sbox-fwd-g128-a34-d16-ad4.slp">slp</a></td>
      <td>Forward</td>
      <td align="right">34</td>
      <td align="right">128</td>
      <td align="right">16</td>
      <td align="right">94</td>
      <td align="right">90</td>
      <td align="right">4</td>
    </tr>
    <tr>
      <td><a href="./aes-sbox-rev-g121-a34-d21-ad4.slp">slp</a></td>
      <td>Inverse</td>
      <td align="right">34</td>
      <td align="right">121</td>
      <td align="right">21</td>
      <td align="right">87</td>
      <td align="right">83</td>
      <td align="right">4</td>
    </tr>
    <tr>
      <td><a href="./aes-sbox-rev-g127-a34-d16-ad4.slp">slp</a></td>
      <td>Inverse</td>
      <td align="right">34</td>
      <td align="right">127</td>
      <td align="right">16</td>
      <td align="right">93</td>
      <td align="right">83</td>
      <td align="right">10</td>
    </tr>
  </tbody>
</table>


</details>

</details>

<details open>
<summary><h2>AES MixColumns</h2></summary>

Example circuits for AES MixColumns (32-bit to 32-bit linear functions), for both Forward and Inverse directions.

Note: For each direction, <sub><sup>🟩</sup></sub> <ins><strong>indicates</strong></ins>
the lowest depth and lowest #XOR (within the selection of circuits).


<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">Direction</th>
      <th align="center">Depth</th>
      <th align="center">#XOR</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">circ.txt</a></td>
      <th scope="rowgroup" rowspan="3" align="center">Forward</th>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>3</strong></ins></td>
      <td align="right">97</td>
    </tr>
    <tr>
      <td><a href="./mixcols/aes-mixcols-10-xor-90-gd4.circ.txt">circ.txt</a></td>
      <td align="right">4</td>
      <td align="right">90</td>
    </tr>
    <tr>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">circ.txt</a></td>
      <td align="right">5</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>88</strong></ins></td>
    </tr>
  </tbody>
  <tbody>
    <tr>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">circ.txt</a></td>
      <th scope="rowgroup" rowspan="4" align="center">Inverse</th>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>5</strong></ins></td>
      <td align="right">163</td>
    </tr>
    <tr>
      <td><a href="./mixcols/aes-mixcols-inv-xor-149-gd6.circ.txt">circ.txt</a></td>
      <td align="right">6</td>
      <td align="right">149</td>
    </tr>
    <tr>
      <td><a href="./mixcols/aes-mixcols-inv-xor-148-gd7.circ.txt">circ.txt</a></td>
      <td align="right">7</td>
      <td align="right">148</td>
    </tr>
    <tr>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">circ.txt</a></td>
      <td align="right">8</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>142</strong></ins></td>
    </tr>
  </tbody>
</table>

</details>

<details open>
<summary><h2>AES folded circuits</h2></summary>

A circuit is "folded" when some components are not "flattened" to a sequence of basic gates. The circuits below call aes-sbox, aes-sbox-inv, aes-mixcols, or aes-mixcols-inv as applicable, and use vector operations to describe XOR-family gates succinctly.

The table counts top-level operations as written in each folded circuit.

<table>
  <thead>
    <tr>
      <th align="center">File<br>name</th>
      <th align="center">Operation</th>
      <th align="center">Key<br>size</th>
      <th align="center">#sbox</th>
      <th align="center">#sbox-<br>inv</th>
      <th align="center">#mixcols</th>
      <th align="center">#mixcols-<br>inv</th>
      <th align="center">#VXOR</th>
      <th align="center">#VXNOR</th>
      <th align="center">#XNOR</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes128-cip-fold1.circ.txt">aes128-cip</a></td>
      <th scope="rowgroup" rowspan="3" align="center">Cipher</th>
      <td align="right">128</td>
      <td align="right">200</td>
      <td align="right">—</td>
      <td align="right">36</td>
      <td align="right">—</td>
      <td align="right">61</td>
      <td align="right">2</td>
      <td align="right">8</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-cip-fold1.circ.txt">aes192-cip</a></td>
      <td align="right">192</td>
      <td align="right">224</td>
      <td align="right">—</td>
      <td align="right">44</td>
      <td align="right">—</td>
      <td align="right">68</td>
      <td align="right">—</td>
      <td align="right">8</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-cip-fold1.circ.txt">aes256-cip</a></td>
      <td align="right">256</td>
      <td align="right">276</td>
      <td align="right">—</td>
      <td align="right">52</td>
      <td align="right">—</td>
      <td align="right">74</td>
      <td align="right">—</td>
      <td align="right">7</td>
    </tr>
  </tbody>
  <tbody>
    <tr>
      <td><a href="./full/aes128-dec-fold1.circ.txt">aes128-dec</a></td>
      <th scope="rowgroup" rowspan="3" align="center">InvCipher</th>
      <td align="right">128</td>
      <td align="right">40</td>
      <td align="right">160</td>
      <td align="right">—</td>
      <td align="right">36</td>
      <td align="right">61</td>
      <td align="right">2</td>
      <td align="right">8</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-dec-fold1.circ.txt">aes192-dec</a></td>
      <td align="right">192</td>
      <td align="right">32</td>
      <td align="right">192</td>
      <td align="right">—</td>
      <td align="right">44</td>
      <td align="right">68</td>
      <td align="right">—</td>
      <td align="right">8</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-dec-fold1.circ.txt">aes256-dec</a></td>
      <td align="right">256</td>
      <td align="right">52</td>
      <td align="right">224</td>
      <td align="right">—</td>
      <td align="right">52</td>
      <td align="right">74</td>
      <td align="right">—</td>
      <td align="right">7</td>
    </tr>
  </tbody>
</table>

</details>

<details open>
<summary><h2>AES flat circuits</h2></summary>

Example circuits, flattened to the level of basic Boolean gates (AND, XOR, XNOR).

There are up to 4!=24 tuple metrics from possible orderings of (A, AD, G, GD), from which a corresponding circuit can be produced,  on corresponding components. For succinctness of presentation, the tables consider only four tuple-metrics (TM): (1) A-AD-G-GD, (2) AD-GD-G-A, (3) GD-G-AD-A, and (4) G-A-GD-AD (4).

Note: In selected columns, <sub><sup>🟩</sup></sub> <ins><strong>indicates</strong></ins> the lowest displayed value.

<details open>
<summary><h3>AES-128 Cipher (cip): Flat circuits</h3></summary>

<details open>
<summary><h4>Based on SBox filtered to #AND ≤ 34</h4></summary>

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">|k|</th>
      <th align="center">A</th>
      <th align="center">AD</th>
      <th align="center">G</th>
      <th align="center">GD</th>
      <th align="center">XX</th>
      <th align="center">X</th>
      <th align="center">X'</th>
      <th align="center">S-box</th>
      <th align="center">MixCols</th>
      <th align="center">TM</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes128-cip-a5800-ad50-g33656-gd388-xx27856-4016.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>5800</strong></ins></td>
      <td align="right">50</td>
      <td align="right">33656</td>
      <td align="right">388</td>
      <td align="right">27856</td>
      <td align="right">23840</td>
      <td align="right">4016</td>
      <td><a href="./sbox/aes-sbox-fwd-a29-ad5-g139-gd35-xx110-20.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-cip-a6800-ad40-g31780-gd191-xx24980-816.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right">6800</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>40</strong></ins></td>
      <td align="right">31780</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>191</strong></ins></td>
      <td align="right">24980</td>
      <td align="right">24164</td>
      <td align="right">816</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-cip-a6800-ad40-g31780-gd191-xx24980-816.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right">6800</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>40</strong></ins></td>
      <td align="right">31780</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>191</strong></ins></td>
      <td align="right">24980</td>
      <td align="right">24164</td>
      <td align="right">816</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-cip-a6400-ad50-g27856-gd286-xx21456-616.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right">6400</td>
      <td align="right">50</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>27856</strong></ins></td>
      <td align="right">286</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>21456</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>20840</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>616</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a32-ad5-g110-gd23-xx78-3.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

<details open>
<summary><h4>Based on SBox filtered to #AND &gt; 34</h4></summary>

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">|k|</th>
      <th align="center">A</th>
      <th align="center">AD</th>
      <th align="center">G</th>
      <th align="center">GD</th>
      <th align="center">XX</th>
      <th align="center">X</th>
      <th align="center">X'</th>
      <th align="center">S-box</th>
      <th align="center">MixCols</th>
      <th align="center">TM</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes128-cip-a7200-ad40-g33456-gd196-xx26256-4416.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>7200</strong></ins></td>
      <td align="right">40</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>33456</strong></ins></td>
      <td align="right">196</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>26256</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>21840</strong></ins></td>
      <td align="right">4416</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-cip-a9400-ad30-g51180-gd191-xx41780-816.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right">9400</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>30</strong></ins></td>
      <td align="right">51180</td>
      <td align="right">191</td>
      <td align="right">41780</td>
      <td align="right">40964</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>816</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a47-ad3-g225-gd15-xx178-4.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-cip-a9000-ad40-g36380-gd161-xx27380-5216.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right">9000</td>
      <td align="right">40</td>
      <td align="right">36380</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>161</strong></ins></td>
      <td align="right">27380</td>
      <td align="right">22164</td>
      <td align="right">5216</td>
      <td><a href="./sbox/aes-sbox-fwd-a45-ad4-g151-gd12-xx106-26.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-cip-a7200-ad40-g33456-gd196-xx26256-4416.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>7200</strong></ins></td>
      <td align="right">40</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>33456</strong></ins></td>
      <td align="right">196</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>26256</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>21840</strong></ins></td>
      <td align="right">4416</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

</details>

<details open>
<summary><h3>AES-128 Decipher (dec): Flat circuits</h3></summary>

<details open>
<summary><h4>Based on SBox filtered to #AND ≤ 34</h4></summary>

<table>
  <thead>
    <tr>
      <th rowspan="2" align="center">File</th>
      <th rowspan="2" align="center">|k|</th>
      <th rowspan="2" align="center">A</th>
      <th rowspan="2" align="center">AD</th>
      <th rowspan="2" align="center">G</th>
      <th rowspan="2" align="center">GD</th>
      <th rowspan="2" align="center">XX</th>
      <th rowspan="2" align="center">X</th>
      <th rowspan="2" align="center">X'</th>
      <th colspan="2" align="center">S-box</th>
      <th rowspan="2" align="center">MixCols-Inv</th>
      <th rowspan="2" align="center">TM</th>
    </tr>
    <tr><th align="center">Fwd</th><th align="center">Inv</th></tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes128-dec-a5800-ad100-g36560-gd784-xx30760-5456.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>5800</strong></ins></td>
      <td align="right">100</td>
      <td align="right">36560</td>
      <td align="right">784</td>
      <td align="right">30760</td>
      <td align="right">25304</td>
      <td align="right">5456</td>
      <td><a href="./sbox/aes-sbox-fwd-a29-ad5-g139-gd35-xx110-20.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a29-ad5-g145-gd32-xx116-29.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-dec-a6800-ad80-g35116-gd396-xx28316-3056.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right">6800</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>80</strong></ins></td>
      <td align="right">35116</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>396</strong></ins></td>
      <td align="right">28316</td>
      <td align="right">25260</td>
      <td align="right">3056</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a34-ad4-g134-gd15-xx100-18.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-dec-a6800-ad80-g35116-gd396-xx28316-3056.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right">6800</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>80</strong></ins></td>
      <td align="right">35116</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>396</strong></ins></td>
      <td align="right">28316</td>
      <td align="right">25260</td>
      <td align="right">3056</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a34-ad4-g134-gd15-xx100-18.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-dec-a6400-ad100-g30120-gd614-xx23720-1576.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right">6400</td>
      <td align="right">100</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>30120</strong></ins></td>
      <td align="right">614</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>23720</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>22144</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>1576</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a32-ad5-g110-gd23-xx78-3.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a32-ad5-g112-gd27-xx80-9.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

<details open>
<summary><h4>Based on SBox filtered to #AND &gt; 34</h4></summary>

<table>
  <thead>
    <tr>
      <th rowspan="2" align="center">File</th>
      <th rowspan="2" align="center">|k|</th>
      <th rowspan="2" align="center">A</th>
      <th rowspan="2" align="center">AD</th>
      <th rowspan="2" align="center">G</th>
      <th rowspan="2" align="center">GD</th>
      <th rowspan="2" align="center">XX</th>
      <th rowspan="2" align="center">X</th>
      <th rowspan="2" align="center">X'</th>
      <th colspan="2" align="center">S-box</th>
      <th rowspan="2" align="center">MixCols-Inv</th>
      <th rowspan="2" align="center">TM</th>
    </tr>
    <tr><th align="center">Fwd</th><th align="center">Inv</th></tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes128-dec-a7200-ad80-g37000-gd403-xx29800-4736.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>7200</strong></ins></td>
      <td align="right">80</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>37000</strong></ins></td>
      <td align="right">403</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29800</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>25064</strong></ins></td>
      <td align="right">4736</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-dec-a7640-ad70-g41236-gd386-xx33596-4016.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right">7640</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>70</strong></ins></td>
      <td align="right">41236</td>
      <td align="right">386</td>
      <td align="right">33596</td>
      <td align="right">29580</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>4016</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a47-ad3-g225-gd15-xx178-4.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-dec-a7560-ad80-g38276-gd356-xx30716-4896.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right">7560</td>
      <td align="right">80</td>
      <td align="right">38276</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>356</strong></ins></td>
      <td align="right">30716</td>
      <td align="right">25820</td>
      <td align="right">4896</td>
      <td><a href="./sbox/aes-sbox-fwd-a45-ad4-g151-gd12-xx106-26.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes128-dec-a7200-ad80-g37000-gd403-xx29800-4736.circ.txt">circ.txt</a></td>
      <td align="right">128</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>7200</strong></ins></td>
      <td align="right">80</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>37000</strong></ins></td>
      <td align="right">403</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29800</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>25064</strong></ins></td>
      <td align="right">4736</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

</details>

<details open>
<summary><h3>AES-192 Cipher (cip): Flat circuits</h3></summary>

<details open>
<summary><h4>Based on SBox filtered to #AND ≤ 34</h4></summary>

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">|k|</th>
      <th align="center">A</th>
      <th align="center">AD</th>
      <th align="center">G</th>
      <th align="center">GD</th>
      <th align="center">XX</th>
      <th align="center">X</th>
      <th align="center">X'</th>
      <th align="center">S-box</th>
      <th align="center">MixCols</th>
      <th align="center">TM</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes192-cip-a6496-ad60-g38144-gd466-xx31648-4488.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>6496</strong></ins></td>
      <td align="right">60</td>
      <td align="right">38144</td>
      <td align="right">466</td>
      <td align="right">31648</td>
      <td align="right">27160</td>
      <td align="right">4488</td>
      <td><a href="./sbox/aes-sbox-fwd-a29-ad5-g139-gd35-xx110-20.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-cip-a7616-ad48-g36076-gd226-xx28460-904.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right">7616</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>48</strong></ins></td>
      <td align="right">36076</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>226</strong></ins></td>
      <td align="right">28460</td>
      <td align="right">27556</td>
      <td align="right">904</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-cip-a7616-ad48-g36076-gd226-xx28460-904.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right">7616</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>48</strong></ins></td>
      <td align="right">36076</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>226</strong></ins></td>
      <td align="right">28460</td>
      <td align="right">27556</td>
      <td align="right">904</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-cip-a7168-ad60-g31648-gd344-xx24480-680.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right">7168</td>
      <td align="right">60</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>31648</strong></ins></td>
      <td align="right">344</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>24480</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>23800</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>680</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a32-ad5-g110-gd23-xx78-3.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

<details open>
<summary><h4>Based on SBox filtered to #AND &gt; 34</h4></summary>

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">|k|</th>
      <th align="center">A</th>
      <th align="center">AD</th>
      <th align="center">G</th>
      <th align="center">GD</th>
      <th align="center">XX</th>
      <th align="center">X</th>
      <th align="center">X'</th>
      <th align="center">S-box</th>
      <th align="center">MixCols</th>
      <th align="center">TM</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes192-cip-a8064-ad48-g37920-gd236-xx29856-4936.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>8064</strong></ins></td>
      <td align="right">48</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>37920</strong></ins></td>
      <td align="right">236</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29856</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>24920</strong></ins></td>
      <td align="right">4936</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-cip-a10528-ad36-g57804-gd226-xx47276-904.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right">10528</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>36</strong></ins></td>
      <td align="right">57804</td>
      <td align="right">226</td>
      <td align="right">47276</td>
      <td align="right">46372</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>904</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a47-ad3-g225-gd15-xx178-4.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-cip-a10080-ad48-g41228-gd190-xx31148-5832.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right">10080</td>
      <td align="right">48</td>
      <td align="right">41228</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>190</strong></ins></td>
      <td align="right">31148</td>
      <td align="right">25316</td>
      <td align="right">5832</td>
      <td><a href="./sbox/aes-sbox-fwd-a45-ad4-g151-gd12-xx106-26.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-cip-a8064-ad48-g37920-gd236-xx29856-4936.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>8064</strong></ins></td>
      <td align="right">48</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>37920</strong></ins></td>
      <td align="right">236</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29856</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>24920</strong></ins></td>
      <td align="right">4936</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

</details>

<details open>
<summary><h3>AES-192 Decipher (dec): Flat circuits</h3></summary>

<details open>
<summary><h4>Based on SBox filtered to #AND ≤ 34</h4></summary>

<table>
  <thead>
    <tr>
      <th rowspan="2" align="center">File</th>
      <th rowspan="2" align="center">|k|</th>
      <th rowspan="2" align="center">A</th>
      <th rowspan="2" align="center">AD</th>
      <th rowspan="2" align="center">G</th>
      <th rowspan="2" align="center">GD</th>
      <th rowspan="2" align="center">XX</th>
      <th rowspan="2" align="center">X</th>
      <th rowspan="2" align="center">X'</th>
      <th colspan="2" align="center">S-box</th>
      <th rowspan="2" align="center">MixCols-Inv</th>
      <th rowspan="2" align="center">TM</th>
    </tr>
    <tr><th align="center">Fwd</th><th align="center">Inv</th></tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes192-dec-a6496-ad100-g41672-gd804-xx35176-6216.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>6496</strong></ins></td>
      <td align="right">100</td>
      <td align="right">41672</td>
      <td align="right">804</td>
      <td align="right">35176</td>
      <td align="right">28960</td>
      <td align="right">6216</td>
      <td><a href="./sbox/aes-sbox-fwd-a29-ad5-g139-gd35-xx110-20.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a29-ad5-g145-gd32-xx116-29.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-dec-a7616-ad80-g40132-gd414-xx32516-3592.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right">7616</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>80</strong></ins></td>
      <td align="right">40132</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>414</strong></ins></td>
      <td align="right">32516</td>
      <td align="right">28924</td>
      <td align="right">3592</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a34-ad4-g134-gd15-xx100-18.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-dec-a7616-ad80-g40132-gd414-xx32516-3592.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right">7616</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>80</strong></ins></td>
      <td align="right">40132</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>414</strong></ins></td>
      <td align="right">32516</td>
      <td align="right">28924</td>
      <td align="right">3592</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a34-ad4-g134-gd15-xx100-18.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-dec-a7168-ad100-g34408-gd644-xx27240-1832.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right">7168</td>
      <td align="right">100</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>34408</strong></ins></td>
      <td align="right">644</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>27240</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>25408</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>1832</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a32-ad5-g110-gd23-xx78-3.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a32-ad5-g112-gd27-xx80-9.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

<details open>
<summary><h4>Based on SBox filtered to #AND &gt; 34</h4></summary>

<table>
  <thead>
    <tr>
      <th rowspan="2" align="center">File</th>
      <th rowspan="2" align="center">|k|</th>
      <th rowspan="2" align="center">A</th>
      <th rowspan="2" align="center">AD</th>
      <th rowspan="2" align="center">G</th>
      <th rowspan="2" align="center">GD</th>
      <th rowspan="2" align="center">XX</th>
      <th rowspan="2" align="center">X</th>
      <th rowspan="2" align="center">X'</th>
      <th colspan="2" align="center">S-box</th>
      <th rowspan="2" align="center">MixCols-Inv</th>
      <th rowspan="2" align="center">TM</th>
    </tr>
    <tr><th align="center">Fwd</th><th align="center">Inv</th></tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes192-dec-a8064-ad80-g42216-gd427-xx34152-5320.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>8064</strong></ins></td>
      <td align="right">80</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>42216</strong></ins></td>
      <td align="right">427</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>34152</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>28832</strong></ins></td>
      <td align="right">5320</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-dec-a8416-ad72-g45924-gd402-xx37508-4744.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right">8416</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>72</strong></ins></td>
      <td align="right">45924</td>
      <td align="right">402</td>
      <td align="right">37508</td>
      <td align="right">32764</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>4744</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a47-ad3-g225-gd15-xx178-4.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-dec-a8352-ad80-g43556-gd378-xx35204-5448.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right">8352</td>
      <td align="right">80</td>
      <td align="right">43556</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>378</strong></ins></td>
      <td align="right">35204</td>
      <td align="right">29756</td>
      <td align="right">5448</td>
      <td><a href="./sbox/aes-sbox-fwd-a45-ad4-g151-gd12-xx106-26.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes192-dec-a8064-ad80-g42216-gd427-xx34152-5320.circ.txt">circ.txt</a></td>
      <td align="right">192</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>8064</strong></ins></td>
      <td align="right">80</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>42216</strong></ins></td>
      <td align="right">427</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>34152</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>28832</strong></ins></td>
      <td align="right">5320</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

</details>

<details open>
<summary><h3>AES-256 Cipher (cip): Flat circuits</h3></summary>

<details open>
<summary><h4>Based on SBox filtered to #AND ≤ 34</h4></summary>

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">|k|</th>
      <th align="center">A</th>
      <th align="center">AD</th>
      <th align="center">G</th>
      <th align="center">GD</th>
      <th align="center">XX</th>
      <th align="center">X</th>
      <th align="center">X'</th>
      <th align="center">S-box</th>
      <th align="center">MixCols</th>
      <th align="center">TM</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes256-cip-a8004-ad70-g46524-gd544-xx38520-5527.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>8004</strong></ins></td>
      <td align="right">70</td>
      <td align="right">46524</td>
      <td align="right">544</td>
      <td align="right">38520</td>
      <td align="right">32993</td>
      <td align="right">5527</td>
      <td><a href="./sbox/aes-sbox-fwd-a29-ad5-g139-gd35-xx110-20.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-cip-a9384-ad56-g43956-gd264-xx34572-1111.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right">9384</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>56</strong></ins></td>
      <td align="right">43956</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>264</strong></ins></td>
      <td align="right">34572</td>
      <td align="right">33461</td>
      <td align="right">1111</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-cip-a9384-ad56-g43956-gd264-xx34572-1111.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right">9384</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>56</strong></ins></td>
      <td align="right">43956</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>264</strong></ins></td>
      <td align="right">34572</td>
      <td align="right">33461</td>
      <td align="right">1111</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-cip-a8832-ad70-g38520-gd402-xx29688-835.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right">8832</td>
      <td align="right">70</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>38520</strong></ins></td>
      <td align="right">402</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>29688</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>28853</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>835</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a32-ad5-g110-gd23-xx78-3.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

<details open>
<summary><h4>Based on SBox filtered to #AND &gt; 34</h4></summary>

<table>
  <thead>
    <tr>
      <th align="center">File</th>
      <th align="center">|k|</th>
      <th align="center">A</th>
      <th align="center">AD</th>
      <th align="center">G</th>
      <th align="center">GD</th>
      <th align="center">XX</th>
      <th align="center">X</th>
      <th align="center">X'</th>
      <th align="center">S-box</th>
      <th align="center">MixCols</th>
      <th align="center">TM</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes256-cip-a9936-ad56-g46248-gd276-xx36312-6079.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>9936</strong></ins></td>
      <td align="right">56</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>46248</strong></ins></td>
      <td align="right">276</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>36312</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>30233</strong></ins></td>
      <td align="right">6079</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-cip-a12972-ad42-g70728-gd264-xx57756-1111.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right">12972</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>42</strong></ins></td>
      <td align="right">70728</td>
      <td align="right">264</td>
      <td align="right">57756</td>
      <td align="right">56645</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>1111</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a47-ad3-g225-gd15-xx178-4.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-cip-a12420-ad56-g50304-gd222-xx37884-7183.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right">12420</td>
      <td align="right">56</td>
      <td align="right">50304</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>222</strong></ins></td>
      <td align="right">37884</td>
      <td align="right">30701</td>
      <td align="right">7183</td>
      <td><a href="./sbox/aes-sbox-fwd-a45-ad4-g151-gd12-xx106-26.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-xor-97-gd3.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-cip-a9936-ad56-g46248-gd276-xx36312-6079.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>9936</strong></ins></td>
      <td align="right">56</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>46248</strong></ins></td>
      <td align="right">276</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>36312</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>30233</strong></ins></td>
      <td align="right">6079</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-04-xor-88-gd5.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

</details>

<details open>
<summary><h3>AES-256 Decipher (dec): Flat circuits</h3></summary>

<details open>
<summary><h4>Based on SBox filtered to #AND ≤ 34</h4></summary>

<table>
  <thead>
    <tr>
      <th rowspan="2" align="center">File</th>
      <th rowspan="2" align="center">|k|</th>
      <th rowspan="2" align="center">A</th>
      <th rowspan="2" align="center">AD</th>
      <th rowspan="2" align="center">G</th>
      <th rowspan="2" align="center">GD</th>
      <th rowspan="2" align="center">XX</th>
      <th rowspan="2" align="center">X</th>
      <th rowspan="2" align="center">X'</th>
      <th colspan="2" align="center">S-box</th>
      <th rowspan="2" align="center">MixCols-Inv</th>
      <th rowspan="2" align="center">TM</th>
    </tr>
    <tr><th align="center">Fwd</th><th align="center">Inv</th></tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes256-dec-a8004-ad135-g50676-gd1062-xx42672-7543.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>8004</strong></ins></td>
      <td align="right">135</td>
      <td align="right">50676</td>
      <td align="right">1062</td>
      <td align="right">42672</td>
      <td align="right">35129</td>
      <td align="right">7543</td>
      <td><a href="./sbox/aes-sbox-fwd-a29-ad5-g139-gd35-xx110-20.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a29-ad5-g145-gd32-xx116-29.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-dec-a9384-ad108-g48732-gd537-xx39348-4247.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right">9384</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>108</strong></ins></td>
      <td align="right">48732</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>537</strong></ins></td>
      <td align="right">39348</td>
      <td align="right">35101</td>
      <td align="right">4247</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a34-ad4-g134-gd15-xx100-18.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-dec-a9384-ad108-g48732-gd537-xx39348-4247.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right">9384</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>108</strong></ins></td>
      <td align="right">48732</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>537</strong></ins></td>
      <td align="right">39348</td>
      <td align="right">35101</td>
      <td align="right">4247</td>
      <td><a href="./sbox/aes-sbox-fwd-a34-ad4-g128-gd15-xx94-4.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a34-ad4-g134-gd15-xx100-18.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-dec-a8832-ad135-g41776-gd835-xx32944-2179.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right">8832</td>
      <td align="right">135</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>41776</strong></ins></td>
      <td align="right">835</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>32944</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>30765</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>2179</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a32-ad5-g110-gd23-xx78-3.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a32-ad5-g112-gd27-xx80-9.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

<details open>
<summary><h4>Based on SBox filtered to #AND &gt; 34</h4></summary>

<table>
  <thead>
    <tr>
      <th rowspan="2" align="center">File</th>
      <th rowspan="2" align="center">|k|</th>
      <th rowspan="2" align="center">A</th>
      <th rowspan="2" align="center">AD</th>
      <th rowspan="2" align="center">G</th>
      <th rowspan="2" align="center">GD</th>
      <th rowspan="2" align="center">XX</th>
      <th rowspan="2" align="center">X</th>
      <th rowspan="2" align="center">X'</th>
      <th colspan="2" align="center">S-box</th>
      <th rowspan="2" align="center">MixCols-Inv</th>
      <th rowspan="2" align="center">TM</th>
    </tr>
    <tr><th align="center">Fwd</th><th align="center">Inv</th></tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="./full/aes256-dec-a9936-ad108-g51296-gd549-xx41360-6527.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>9936</strong></ins></td>
      <td align="right">108</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>51296</strong></ins></td>
      <td align="right">549</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>41360</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>34833</strong></ins></td>
      <td align="right">6527</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">1</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-dec-a10508-ad95-g56912-gd523-xx46404-5591.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right">10508</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>95</strong></ins></td>
      <td align="right">56912</td>
      <td align="right">523</td>
      <td align="right">46404</td>
      <td align="right">40813</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>5591</strong></ins></td>
      <td><a href="./sbox/aes-sbox-fwd-a47-ad3-g225-gd15-xx178-4.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">2</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-dec-a10404-ad108-g53064-gd484-xx42660-6735.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right">10404</td>
      <td align="right">108</td>
      <td align="right">53064</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>484</strong></ins></td>
      <td align="right">42660</td>
      <td align="right">35925</td>
      <td align="right">6735</td>
      <td><a href="./sbox/aes-sbox-fwd-a45-ad4-g151-gd12-xx106-26.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-163-gd5.circ.txt">Link</a></td>
      <td align="right">3</td>
    </tr>
    <tr>
      <td><a href="./full/aes256-dec-a9936-ad108-g51296-gd549-xx41360-6527.circ.txt">circ.txt</a></td>
      <td align="right">256</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>9936</strong></ins></td>
      <td align="right">108</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>51296</strong></ins></td>
      <td align="right">549</td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>41360</strong></ins></td>
      <td align="right"><sub><sup>🟩</sup></sub> <ins><strong>34833</strong></ins></td>
      <td align="right">6527</td>
      <td><a href="./sbox/aes-sbox-fwd-a36-ad4-g138-gd14-xx102-22.circ.txt">Link</a></td>
      <td><a href="./sbox/aes-sbox-inv-a36-ad4-g148-gd14-xx112-24.circ.txt">Link</a></td>
      <td><a href="./mixcols/aes-mixcols-inv-xor-142-gd8.circ.txt">Link</a></td>
      <td align="right">4</td>
    </tr>
  </tbody>
</table>
<small>|k| = key size; A = #AND; AD = AND depth; G = #Gates; GD = gate depth; XX = X + X'; X = #XOR; X' = #XNOR; Fwd = Forward; Inv = Inverse; TM = tuple metric: 1 = A-AD-G-GD; 2 = AD-GD-G-A; 3 = GD-G-AD-A; 4 = G-A-GD-AD.</small>


</details>

</details>


</details>
