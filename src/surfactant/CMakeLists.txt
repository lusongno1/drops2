set(HOME surfactant)

libs(ifacetransp)

target_link_libraries(surfactant-ifacetransp levelset-levelsetmapper)

exec_ser(surfactant geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer num-gradient_recovery levelset-levelsetmapper)

exec_ser(surfpatternfm geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer num-gradient_recovery levelset-levelsetmapper poisson-ale geom-boundary geom-builder out-ensightOut geom-simplex geom-multigrid num-unknowns geom-topo poisson-poisson misc-problem misc-progressaccu misc-utils out-output num-fe num-discretize num-interfacePatch geom-geomselect misc-params out-vtkOut misc-funcmap geom-deformation misc-scopetimer poisson-poissonParam misc-dynamicload)

exec_ser(surfacestokes geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer)

exec(heatequation poisson-ale geom-boundary geom-builder out-ensightOut geom-simplex geom-multigrid num-unknowns geom-topo poisson-poisson misc-problem misc-progressaccu misc-utils out-output num-fe num-discretize num-interfacePatch geom-geomselect misc-params out-vtkOut misc-funcmap geom-deformation misc-scopetimer poisson-poissonParam misc-dynamicload)


if(NOT MPI)
add_dependencies(surfactant misc-scalarFunctions misc-vectorFunctions)
add_dependencies(surfpatternfm misc-scalarFunctions misc-vectorFunctions poisson-poissonCoeff)
add_dependencies(surfacestokes misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff)
add_dependencies(heatequation misc-scalarFunctions misc-vectorFunctions poisson-poissonCoeff)
endif(NOT MPI)

add_my_custom_targets(surfactant surfacestokes)
