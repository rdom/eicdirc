## Instalation
```
git clone https://github.com/rdom/prttools
git clone https://github.com/rdom/eicdirc
cd eicdirc
mkdir build
cd build
cmake ..
make -j4

#test event display
./eicdirc -s 0 -a 30 -x "pi+" -p 6.0 -w 0 -h 11 -c 4 -l 3 -v 0 -zpos 1 -g 1 -ev 0 -e 1
```


## Synopsis
```
eicdirc [OPTION] [ARGUMENT] ... [OPTION] [ARGUMENT]

example:
./eicdirc -a 40 -l 0 -x "pi+" -p 1 -w 0 -g 0 -e 1
```
## Options
```
-o    output file name
-i    input file name
-u    look-up file name

-s    run type
                0    simulation
                1    look-up table generation
                2    reconstruction
                3    likelihood calculation for geometrical reconstruction
                4    likelihood calculation for time imaging
                5    calibration
                10   simulation of pions and kaons for 3,4


-g    geometry configuration
                0    whole DIRC
                1    one barbox (default)


-ev   expansion volume type
	        0    prism with lenses (default)
                1    BaBar wedge with focusing prism
                3    prism with plate
		4    tilted EV	

-h    number of bars in one radiator box

-c    MCP layout
                0    4x6 standard MCPs (pixel size == mcp size)
                1    4x6 standard MCPs (6.4x6.4 pixels)
		3    one MCP cover all FD plain (custom pixel size)
		4    4x6 MCPs with compact paking and 32x32 pixels (2 mm x 2 mm)
-l    focusing system
                0    no lens
                1    spherical lens
                3    3-layer spherical lens
                6    3-layer cylindrical lens
                10   ideal lens (thickness = 0, ideal focusing)

-a    angle between particle beam and bar radiator
      if a == 0 then a = [30,160]

-e    number of simulated events

-x    particle type
              "pi+"
              "proton"
              "kaon"
                 ...
              "opticalphoton"
              "mix_pik"  1 pion 1 kaon mix
	      "mix_pie"  1 pion 1 electron mix
	      "mix_pimu" 1 pion 1 muon mix

-p    particle momentum [GeV/c]

-w    physical list
                0    standard
                1    without multiple scattering, ionisation and bremsstrahlung
                10   monochromatic Cherenkov light
                11   10 + 1 

-r    seed number for the random generator 

-b    batch mode
               1    run silent (without GUI)

-d    display option
               0    standard (default)
               1    display hit occupancy of current run
               2    display hit occupancy of occuhits.root (needs to be generated)

-zpos target's Z offset 

-z    track smearing at vertex [rad]

```

## Example of script usage from macro folder
```
./ba_scan -j6 -r0 -s5 -e50 -t1 -v0
root da_scan.C'("r_spr39498736070.root","ttt1.root")'
```