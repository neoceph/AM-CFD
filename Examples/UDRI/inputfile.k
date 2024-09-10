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
mpPhase
mpLiqPhase
*CONTROL_TERMINATION
$ total_time  max_iter  nonlin_tolerance
 2.0e-3    100       1.0e-6
*CONTROL_TIMESTEP
$ sim_timestep  output_timestep
  1.0e-5	      1e-5
*POWDER_BED_MODE
$ numlayers  ncvlayer  layerheight  density  cpa     cpb     thcona     thconb
  1          10         0.04e-3      4500    0.1913  443.37  0.0138     11.42
$*VOLUMETRIC_HEAT_SOURCE_TABLE
$./heatSourceParametersTable.k
*VOLUMETRIC_HEAT_SOURCE
$ power   apowseta   heatDepth  heatrb   factor
  250.00  0.5e4      60.0e-06   40.0e-06      2
*MATERIAL_PROPERTIES
$ dens  denl   viscos  tsolid  tliquid  hsmelt  hlfriz
  7734  7578   5.30e-3   1533    1609     2.86e6  1.49e6
$ acpa    acpb    acpl
  0.1913  443.37  755
$ thconsa  thconsb   thconsc  thconla thconlb
  0.0000   0.0138    11.42     0.0     31.3
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
./multiTrackmesh_01x01.k
*TOOLPATH_FILE
./MultiLayerTest.crs
