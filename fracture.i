[Mesh]
  type = FileMesh
  file = 'gold/geo.msh'
[]


[Variables]
  [d]
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [Gc]
    order = CONSTANT
    family = MONOMIAL
  []
[]

# [AuxKernels]
#   [psie_active]
#     type = ADMaterialRealAux
#     variable = psie_active
#     property = psie_active
#   []
#   [Gc]
#     type = ADMaterialRealAux
#     variable = Gc
#     property = Gc
#   []
# []

[Bounds]
  [irreversibility]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  []
  [upper]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1
  []
[]

[Kernels]
  [diff]
    type = ADPFFDiffusion
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0
  []
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
[]

[Materials]
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'Gc l  E ft'                   # for wu
    prop_values = '${Gc} ${l} ${E} ${ft}'
    # prop_names = 'Gc l E '                     # for power
    # prop_values = '${Gc} ${l} ${E}'
  []
  # [crack_geometric]
  #   type = CrackGeometricFunction
  #   f_name = alpha
  #   function = 'd^2'
  #   phase_field = d
  # []
  # [degradation]
  #   type = PowerDegradationFunction
  #   f_name = g
  #   function = (1-d)^p*(1-eta)+eta
  #   phase_field = d
  #   parameter_names = 'p eta '
  #   parameter_values = '2 1e-8'
  # []
  [crack_geometric]
    type = CrackGeometricFunction
    f_name = alpha
    function = '2*d-d^2'
    phase_field = d
  []
  [degradation]
    type = WUDegradationFunction
    f_name = g
    # function = (1-d)^2/((1-d)^2+a1*d*(1-0.5*d))
    function = (1-d)^2/((1-d)^2+(4*Gc*E/c0/l/ft^2)*d*(1-0.5*d))*(1-eta)+eta
    phase_field = d
    # material_property_names = 'a1 l '
    material_property_names = 'Gc E ft c0 l '
    parameter_names = 'eta '
    parameter_values = '0 '
  []
  # [degradation]
  #   type = RationalDegradationFunction
  #   f_name = g
  #   function = (1-d)^2/((1-d)^2+(Gc/psic*1/c0/l)*d*(1-0.5*d))*(1-eta)+eta
  #   phase_field = d
  #   material_property_names = 'Gc psic  c0 l '
  #   parameter_names = ' eta '
  #   parameter_values = '  1e-6'
  # []
  [a1]
    type = ADDerivativeParsedMaterial
    f_name = 'a1'
    function = '4*lch/c0/l'
    material_property_names = 'lch c0 l'
  []
  [lch]
    type = ADDerivativeParsedMaterial
    f_name = 'lch'
    function = 'E*Gc/ft^2'
    material_property_names = 'E Gc ft'
  []
  [psi]
    type = ADDerivativeParsedMaterial
    f_name = psi
    function = 'alpha*Gc/c0/l+g*psie_active'
    args = 'd psie_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l'
    derivative_order = 1
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'
  automatic_scaling = true

  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]
