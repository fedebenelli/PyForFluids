# GERG-2008.f90

Program made to calculate the GERG-2008 equation [^1]


## Folder structure

```
.
├── concentrations		-> Each compound concentration
├── gerg.f90			-> GERG-2008
├── modules.f90			-> Subroutines to read the stored parameters
│   				     and calculate thermodynamic properties
│   				     
├── parameters			
│   │   
│   ├── compounds		-> List with each compound index
│   ├── critical		-> Pure compound properties
│   ├── departure		-> Folder with departure function parameters
│   │   ├── betaij		    All files are structured like:
│   │   ├── dij                      compound1_index compound2_index [parameters list]
│   │   ├── epsij
│   │   ├── etaij
│   │   ├── Fij
│   │   ├── gammij
│   │   ├── Kpolij-Kexpij
│   │   ├── nij
│   │   └── tij
│   ├── ideal			-> Folder with ideal gas parameters
│   │   ├── n			    All files are structured like:
│   │   └── v                        compound_index [parameters list]
│   ├── reducing_params		-> Reducing function parameters
│   └── residual		-> Residual Helmholtz Energy parameters
│       ├── c			    All files are structured like:
│       ├── d			     compound_index [parameters list]
│       ├── Kpol-Kexp
│       ├── n
│       └── t
└── README.md			-> This readme file
```

Right now it doesn't work like this, all the parameters are defined in the 
`parameters.f90` file, for more controlled testing.

[^1]: ![Paper link](https://pubs.acs.org/doi/10.1021/je300655b)
