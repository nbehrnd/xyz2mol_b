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

# background

Based on an algorithm by Kim and Kim[^2], Jan Jensen implemented
xyz2mol[^3] relying on RDKit. As presented in 2020
([recording](https://www.youtube.com/watch?v=HD6IpXMVKeo), [pdf
slides](https://github.com/rdkit/UGM_2020/blob/master/Presentations/JanJensen.pdf)),
the approach works best on small molecules which do not carry
organometallic bonds.

With release 2202.09.3, the functionality is available in RDKit.[^4]
This allows a shorter (moderator) script to convert xyz into .sdf
provided one actually has access to a version of RDKit modern enough:
contrasting to Greg Landrum's blog,[^5] RDKit 2022.09.3 as so far
packaged for Linux Debian does not yet include this functionality –
according to repology.org,[^6] this equally applies to relatives via
DebiChem (e.g. Ubuntu ecosystem). As checked, creating a virtual
environment for Python and a pip based installation of RDKit 2023.03.2
(and its dependencies such as numpy) offers one successful bypass. An
alternative environment not tested may be anaconda.

The script was written in an instance of Linux Debian 13/trixie (branch
testing) with Python (version 3.11.4) and pip installed RDKit
(2202.09.3). As inspected with openbabel and Jmol, Jan Jensen's four
test structures pass the conversion.

| file name              | charge |
|------------------------|-------:|
| acetate.xyz            |     -1 |
| ethane.xyz             |      0 |
| propylbenzene.xyz      |      1 |
| stereogenic_center.xyz |      0 |

[^1]: <https://github.com/jensengroup/xyz2mol>

[^2]: Kim, Y and Kim, W. Y. Universal Structure Conversion Method for
    Organic Molecules: From Atomic Connectivity to Three-Dimensional
    Geometry. *Bull. Korean Chem. Soc.* **2015**, *36*, 1769-1777, [doi
    10.1002/bkcs.10334](https://doi.org/10.1002/bkcs.10334).

[^3]: <https://github.com/jensengroup/xyz2mol>

[^4]: <https://github.com/jensengroup/xyz2mol/issues/40>

[^5]: <https://greglandrum.github.io/rdkit-blog/posts/2022-12-18-introducing-rdDetermineBonds.html>

[^6]: <https://repology.org/project/rdkit/packages>
