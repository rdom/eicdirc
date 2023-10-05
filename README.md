## Installation
```
#install and initialize geant4 with -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_QT=ON

git clone https://github.com/rdom/prttools
git clone https://github.com/rdom/eicdirc
cd eicdirc
mkdir build
cd build
cmake ..
make -j4

#neural network
install tensorflow libs, cppflow and run:
cmake -DAI=1 ..


#test event display
./eicdirc -r 0 -theta 30 -x "pi+" -p 6.0 -w 0 -h 11 -c 2031 -l 3 -v 0 -gz 1 -g 1 -ev 0 -e 1
```


## Synopsis
```
eicdirc [OPTION] [ARGUMENT] ... [OPTION] [ARGUMENT]

example:
./eicdirc -r 0 -theta 30 -x "pi+" -p 6.0 -w 0 -h 11 -c 2031 -l 3 -v 0 -gz 1 -g 1 -ev 0 -e 1
```
## Options
```
-o    output file name
-i    input file name
-u    look-up file name
-pdf  look-up file name
-nn   path to the neural network

-r    run type
                0    simulation
                1    look-up table generation
                2    geometrical reconstruction
                3    likelihood calculation
		4    create pdf
		5    simulate pdf

-field    field type
                0    no field (default)
                1    CORE b-field
                2    ePIC MARCO 1.7T
                3    ePIC MARCO 2.0T
		4    solenoidal 3.0T

-g    geometry configuration
                0    ATHENA one barbox
                1    ePIC one barbox
		2    CORE one barbox
		10   ATHENA whole DIRC
		11   ePIC whole DIRC
		12   CORE whole DIRC
		
-ev   expansion volume type
	        0    prism with lenses (default)
                1    BaBar wedge with focusing prism
                3    prism with plate, lens between bars and plate
		4    tilted EV
		5    prism with plate, lens between plate and EV

-h    number of bars in one radiator box

-c    MCP layout
                0    4x6 standard MCPs (pixel size == mcp size)
                1    4x6 standard MCPs (6.4x6.4 pixels)
		3    one MCP cover all FD plain (custom pixel size)
		4    2x3 LAPD
		2031    4x6 MCPs with compact packing and 32x32 pixels (2 mm x 2 mm)
		
-l    focusing system
                0    no lens
                1    spherical lens
                3    3-layer spherical lens
                6    3-layer cylindrical lens
                10   ideal lens (thickness = 0, ideal focusing)

-theta    polar between particle beam and bar radiator [deg]
      if theta == 0 then thata = [30,160]

-phi  azimuth angle between particle beam and bar radiator [deg]

-e    number of simulated events

-x    particle type
              "pi+"
              "proton"
              "kaon+"
                 ...
              "opticalphoton"
	      "mix_pie"  1 pion 1 electron mix
	      "mix_pimu" 1 pion 1 muon mix
              "mix_pik"  1 pion 1 kaon mix
	      "mix_pip"  1 pion 1 proton mix
	      "mix_kp"   1 kaon 1 proton mix
	      
-p    particle momentum [GeV/c]

-w    physical list
                0    standard EM
                1    without multiple scattering and bremsstrahlung
		2    standard EM with HAD (elastic/inelastic/absorption)
		3    FTFP_BERT
                10   monochromatic Cherenkov light
                11   10 + 1 

-seed seed number for the random generator 

-b    batch mode
               1    run silent (without GUI)

-d    display option
               use /Prt/geom/drawHits 2 
               0    standard (default)
               1    display hit occupancy of current run
               2    display hit occupancy of occuhits.root (needs to be generated)

-timeres   time resolution [ns]
               0.2  (default)

-timecut   time cut constant [ns]
               0.5  (default)  

-trackingres   tracking resolution [rad], applied at tracking layer
               0.0005  (default)
	       100 use realistic, mapped values

-gx target's X offset

-gz target's Z offset 


```

## LUT generation

```
eicdirc -o ../data/lut.root -r 1 -g 1 -h 11 -c 2031 -l 3 -v 0 -ev 0 -x "opticalphoton" -p "3.18 eV"  -e 1000000 -b 1
```

Visualization of 100 events:
eicdirc -o ../data/lut.root -r 1 -g 1 -h 11 -c 2031 -l 3 -v 0 -ev 0 -x "opticalphoton" -p "3.18 eV"  -e 100

![alt text](https://github.com/rdom/eicdirc/raw/master/pic/eicdirc_lut_gen.png)


LUT averaging:
```
root -q -b loadlib.C lutmean.C'("../data/lut")'
```

## Simulation:
```
eicdirc -r 1 -o hits.root -r 0 -theta 30 -x "mix_pik" -p 6 -w 0 -g 1 -h 11 -c 2031 -l 3 -trackres 0.0005 -v 0 -gz 1 -ev 0 -b 1 -e 2000
```

## Reconstruction:
```
eicdirc -s 2 -i hits.root -u ../data/lut_avr.root -s 2 -tc 0.5 -v 3 -e 0 -t1 0.0005
```


## Example of script usage from macro folder

Hit pattern:

```
root loadlib.C drawHP.C'("../build/hits.root")'
```
Example of 1k of 6 GeV/c pions @ 30 degree:
![alt text](https://github.com/rdom/eicdirc/raw/master/pic/hp_pi_1k.png)

Angle scan:
```
ba_scan -j6 -r0 -s5 -e50 -t1 -v0
root da_scan.C'("r_spr39498736070.root","ttt1.root")'
```