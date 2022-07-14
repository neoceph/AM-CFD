*KEYWORD_ID
Calibration
*DIRECTORY
Calibration
*OUTPUT_FINAL_GEOMETRY
finalGeom_MultiLayer.txt
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
mpPhase
mpLiqPhase
*CONTROL_TERMINATION
$ total_time  max_iter  nonlin_tolerance
 0.11e-3    100       1.0e-6
*CONTROL_TIMESTEP
$ sim_timestep  output_timestep
  1.0e-5	      1e-5
*POWDER_BED_MODE
$ numlayers  ncvlayer  layerheight  density  cpa     cpb     thcona     thconb
  1          10         0.05e-3      7734    1e-6    435     0          11.4
$*VOLUMETRIC_HEAT_SOURCE_TABLE
$./heatSourceParametersTable.k
*VOLUMETRIC_HEAT_SOURCE
$ power   apowseta   heatthick  heatrb   factor
  285.0   0.55       129.14e-6  130.92e-6     2
*MATERIAL_PROPERTIES
$ dens  denl   viscos  tsolid  tliquid  hsmelt  hlfriz
  7734  7578   5.30e-3   1533    1609     2.86e6  1.49e6
$ acpa    acpb    acpl
  1e-6  435     755
$ thconsa  thconsb   thconsc  thconla thconlb
  0.0000   0.0000    11.4     0.0     31.3
*UNDER_RELAXATION
$ urfu  urfv  urfw  urfp  urfh
  0.8   0.8   0.8   0.8   0.5
*BOUNDARY_TEMPERATURES
$ west  east  north  bottom  preheat
  353   353   353    353     353
*BOUNDARY_ENERGY_EQUATION
$ htci  htcj  htckl  htckn  emiss
  10     10     10     10   0.3
*BOUNDARY_MOMENTUM_EQUATION
$ beta  dgdtp  darcyKo  a_oxygen
  1.3e-5  -3.7e-4  1.6e4    -1
$*ACTIVATE_MOVING_MESH
$  2 2
*MESH_FILE
../Examples/MultiTrack/multiTrackmesh_01x01.k
*TOOLPATH_FILE
../Examples/MultiTrack/MultiLayerTest.crs
