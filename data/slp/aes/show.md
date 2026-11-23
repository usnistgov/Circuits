# Show-selection criterion 

This file explains the criteria used to select which circuits are linked in the tables of results.  Circuits added to the repository usually follow this criteria at the time of upload, but many eventually become superseded.  The criteria may be adjusted.

Notation: A = #AND; AD = AND depth; G = #Gates; GD = Gate depth; XX = #XOR + #XNOR; L = XXN = XX + #NOT.


## XOR-based aes-[inv]mixcols

1. Exclude GD/G dominated circuits, i.e., where another circuit has (≤GD & <G) or (<GD & ≤G).


## \{AND,XOR,XNOR\}-based aes-[inv]sbox

Note: The described criteria applies separately to the sets aes-sbox and aes-[inv]sbox.

### For A≤36:

1. Exclude:
  - a. A/AD/G/GD dominated circuits, i.e., where another circuit has ≤A & ≤AD & ≤G & ≤GD.
  - b. circuits with G≥160.
  - c. circuits with A in {29,30,31} (for low A we'll choose A28 or A≥32), unless their AD/G/GD would still be dominant over some otherwise selected A32 or A33 circuits.
2. For each A/AD, include non-excluded results for:
  - a. Every GD≤15.
  - b. GD16 iff G[GD16]≤G[GD15]-5. Upload these circuits and omit them from the tables.
  - c. GD20. If GD20 is unavailable, use another GDd, with minimal |20-d|.
  - d. GD25 iff G[GD25]≤G[GD20]-5. If GD25 is unavailable, use another GDd, with minimal |25-d|.
  - e. Minimum Gmin (whichever GD) satisfying (GD[Gmin] - 20 < G[GD20] - Gmin) or (Gmin is the global minimum for this function)
3. Exclude A'/AD'/G'/GD' (even if fitting the above criteria) for which a conceivable A/AD/G'/GD' would not fit in at least one A/AD category with (A<A' & AD≤AD') or (A≤A' & AD<AD').

### For GD≤12 with A≥37:

1. Same criteria as [inv]mixcols (i.e., ignoring A/AD)
2. For AD4/GD11, upload the eligible Pareto points in each direction. Show the lowest-G point and the lowest-A point. Omit the lowest-A point when G[lowest A] > G[lowest G] - (A_max - A_min).


## \{AND,XOR,XNOR\}-based aes-sbox-sel

The criteria for 9-bit to 8-bit aes-sbox-sel circuits that jointly implement aes-sbox and aes-invsbox is the same criteria as for sbox/invsbox, once:
  - increasing every GD by 5
  - increasing every A by 8
  - excluding G≥200
This was used for an initial compilation of sbox-sel circuits, to illustrate the functionality. Better circuits can "easily" stem from [inv]sbox improvement, which are updated more frequently.