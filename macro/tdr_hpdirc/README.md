### Procedure to recreate performance plots of the hpDIRC


#### Software installation

Prerequisites:
```
root (>= 6)
geant4 (>= v11.3.2)
```

Install hpDIRC standalone simulation/reconstruction:
```
git clone https://github.com/rdom/eicdirc
cd eicdirc
mkdir build
cd build
cmake ..
make -j4
```


#### Simulation/Reconstruction

Run simulation for 45k pi/K @ 6 GeV/c for every 5 degree polar angle step in [25,155] degree range (it will take ~1.3 cpu-hours per angle; number of threads can be adjusted inside script):

```
./simulate
```

Create PDFs (probability density functions) using last 40k events:

```
./create_pdf
```

Run Time Imaging reconstruction using first 5k events:

```
./reconstruct
```

Per-angle plots will be stored at "sim_data/reco". Create final plots:

```
root plot.C

```

Reference plots for the photon yield and for the separation power of time imaging reconstruction as a function of polar angle are located in "plots" folder.


