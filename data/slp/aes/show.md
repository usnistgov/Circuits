# Show-selection criterion 

This file explains the criteria used to select which circuits are linked in the tables of results.  Circuits added to the repository usually follow this criteria at the time of upload, but many eventually become superseded.  The criteria may be adjusted.

Notation: A = #AND; AD = AND depth; G = #Gates; GD = Gate depth; XX = #XOR + #XNOR; L = XXN = XX + #NOT.


## XOR-based aes-[inv]mixcols

Note: The described criteria applies separately the sets invmixcols and mixcols.

1. Exclude GD/G dominated circuits, i.e., where another circuit has (≤GD & <G) or (<GD & ≤G).
2. Include a non-dominated GD/G iff (GD - GDmin) ≤ (G[GDmin] - G), where GDmin is the minimum depth in the set of circuits.


## \{AND,XOR,XNOR\}-based aes-[inv]sbox

Note: The described criteria applies separately to the sets aes-sbox and aes-[inv]sbox.

### For A≤36:
1. Exclude: 
   - a. A/AD/G/GD dominated circuits, i.e., when another circuit has ≤A & ≤AD & ≤G & ≤GD.
   - b. circuits with G≥185.
2. For each A/AD, include non-excluded results for:
   - a. Every GD≤15.
   - b. GD16 iff G[GD16]≤G[GD15]-5.
   - c. GD20. If GD20 is unavailable, use another GDd, with minimal |20-d|.
   - d. GD25 iff G[GD25]≤G[GD20]-5. If GD25 is unavailable, use another GDd, with minimal |25-d|.
   - e. Minimum Gmin (whichever GD) satisfying (GD[Gmin] - 20 ≤ G[GD20] - Gmin) or (Gmin is the global minimum for this function)
3. Exclude A'/AD'/G'/GD' (even if fitting the above criteria) for which a conceivable A/AD/G'/GD' would not fit in at least one A/AD category with (A<A' & AD≤AD') or (A≤A' & AD<AD').

### For A≥37: 
1. Same criteria as [inv]mixcols (i.e., A/AD does not matter)


## \{AND,XOR,XNOR\}-based aes-sbox-sel

The criteria for 9-bit to 8-bit aes-sbox-sel circuits that jointly implement aes-sbox and aes-invsbox is the same criteria as for sbox/invsbox, once:
  - increasing every GD by 5
  - increasing every A by 8
  - excluding G>=200 and GD>=31