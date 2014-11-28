eicdirc [OPTION] [ARGUMENT] ... [OPTION] [ARGUMENT]

example:
./eicdirc -a 40 -l 0 -x "pi+" -p 1 -w 0 -g 0 -e 1

Options

-o    output file name
-i    input file name
-u    look-up file name

-s    run type
                0    simulation
                1    look-up table generation
                2    reconstruction
                5    calibration

-g    geometry configuration
                1    in vacuum
                2    in air

-h    radiator type
                1    bar (default)
                2    plate

#-c   MCP layout
#                0    MCP covers all FD plain

#-l    focusing system
#                0    no lens
#                1    spherical lens
#                10   ideal lens (thickness = 0, ideal focusing)

-a    angle between particle beam and bar radiator

-e    number of simulated events

-x    particle type
              "pi+" 
              "proton"
              "kaon"
                 ...
              "opticalphoton"   

-p    particle momentum [GeV/c]

-w    physical list
                0    standard
                1    without multiple scattering
                10   monochromatic Cherenkov light
                11   10 + 1 

-r    seed number for the random generator 

-b    batch mode
               1    run silent (without GUI)

