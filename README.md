<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>

# intent

Given an input structure as .xyz file, the moderator script
`xyz2mol_b.py` provides a representation in the .sdf format. The user is
expected to describe the overall charge by flag `-c` (or `--charge`) if
the molecule is not neutral. Based on Jan Jensen's work on
`xyz2mol`[^1], it depends on a recent release of RDKit (e.g.,
2023.03.2). It is suggested to resolve the dependencies in an activated
virtual environment for Python and a call of

``` shell
pip -r install requirements.txt
```

and eventually use the script in a pattern like

``` shell
python ./xyz2mol_b.py ./acetate.xyz --charge -1
```

For an easier systemwide deployment, the [release
page](https://github.com/nbehrnd/xyz2mol_b/releases) provides a platform
independent Python wheel which in turn resolves dependencies like RDKit
and numpy from the PyPI.

# background

Based on an algorithm by Kim and Kim[^2], Jan Jensen implemented
xyz2mol[^3] relying on RDKit to convert .xyz files into .sdf including
assignment of bond orders different than one. As presented in 2020
([recording](https://www.youtube.com/watch?v=HD6IpXMVKeo), [pdf
slides](https://github.com/rdkit/UGM_2020/blob/master/Presentations/JanJensen.pdf)),
the approach works best on small molecules which do not carry
organometallic bonds. With release of RDKit 2022.09.3, this
functionality became available in RDKit,[^4] and was presented in Greg
Landrum's blog,[^5] too.

This moderator script was written to provide this functionality until
RDKit as packaged by DebiChem (for Debian, Ubuntu, etc) would catch up
(cf. notes in repology[^6]) however can be "handy" for a rapid
conversion if one forgot the required syntax in RDKit.

# user notes

By default, the submitted structure in the input .xyz file is presumed
to be balanced and overall neutral (e.g., `ethane.xyz` in subfolder
`tests`). The implementation equally is capable to "recover" a .sdf file
of charged structures like the acetate anion (file `acetate.xyz`), or
tetrabutylammonium (file `tbab_cation.xyz`) for which the user is
required to explicitly assign the overall charge. As one can check with
file `C9H11.xyz`, one input file can lead to both a .sdf file of a
cation (`--charge 1`), and anion (`-c
  -1`).[^7]

[^1]: <https://github.com/jensengroup/xyz2mol>

[^2]: Kim, Y and Kim, W. Y. Universal Structure Conversion Method for
    Organic Molecules: From Atomic Connectivity to Three-Dimensional
    Geometry. *Bull. Korean Chem. Soc.* **2015**, *36*, 1769-1777, [doi
    10.1002/bkcs.10334](https://doi.org/10.1002/bkcs.10334).

[^3]: <https://github.com/jensengroup/xyz2mol>

[^4]: <https://github.com/jensengroup/xyz2mol/issues/40>

[^5]: <https://greglandrum.github.io/rdkit-blog/posts/2022-12-18-introducing-rdDetermineBonds.html>

[^6]: <https://repology.org/project/rdkit/packages>

[^7]: The same Hill formula equally applies to the neutral
    2-phenyl-2-propyl radical, PubChem [CID
    140141](https://pubchem.ncbi.nlm.nih.gov/compound/140141).
