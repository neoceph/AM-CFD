*KEYWORD_ID
Trial_SingleTrack
*DIRECTORY
Trial_SingleTrack
*OUTPUT_FINAL_GEOMETRY
finalGeom_SingleTrack.txt
*PARAMETERS
accelGravity 9.8
radBoltzmann 5.67e-8
ambientTemp  295
ambientPres  101.0e3
gasConstant  8.314
*SCALAR_OUT
temp
uvel
vvel
wvel
solfrac
csfrac
rho
diff
*CONTROL_TERMINATION
$ total_time  max_iter  nonlin_tolerance
  0.00067     100       5.0e-4
*CONTROL_TIMESTEP
$ sim_timestep  output_timestep
  2.0e-5        2.0e-5
*SINGLE_TRACK_TOOL
$ lasertime  xstart     ystart      vel
  1.000      0.500e-3   1.000e-3    1.5
*POWDER_BED_MODE
$ numlayers  ncvlayer  layerheight  density  cpa     cpb     thcona     thconb
  1          8         0.04e-3      4330     0.2508  357.7   0          0.995
*GAUSS_LASER_PARAMETERS
$ power  absorb   radius   factor
  100.0   0.20    40e-6    2
*MATERIAL_PROPERTIES
$ dens  denl   viscos  tsolid  tliquid  hsmelt  hlfriz
  8440  7640   0.007   1563    1623     861e3   1151e3
$ acpa    acpb    acpl
  0.2441  338.59  709.25
$ thconsa  thconsb   thconsc  thconla thconlb
  0.0      0.0155    5.0435   0.0     30.078
*UNDER_RELAXATION
$ urfu  urfv  urfw  urfp  urfh
  0.8   0.8   0.8   0.8   0.8
*BOUNDARY_TEMPERATURES
$ west  east  north  bottom  preheat
  353   353   353    353     353
*BOUNDARY_ENERGY_EQUATION
$ htci  htcj  htckl  htckn  emiss
  10     10     10     10   0.4
*BOUNDARY_MOMENTUM_EQUATION
$ beta  dgdtp  darcyKo  a_oxygen
  5e-5  -3.8e-4  1.6e4    -1
*MESH_FILE
../../MeshFiles/multiTrackmesh_01x01.k
