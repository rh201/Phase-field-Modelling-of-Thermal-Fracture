# use quenching test material parameters in zhuang's paper
E = 3.4e5         # MPa
nu = 0.22
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
Gc = 42.47e-3     # 10^3 J/m^2
l = 0.01         # 10^-3 m
rho = 2450e-12    # 10^12 kg/m^3
c = 0.775e6       # 10^-6 J/kgK
k0 = 300          # W/mk
alpha_bulk = 8.0e-6
ft = 180        # MPa, NO ft for PowerDegradationFunction, # delete ft in material

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc=${Gc};l=${l};E=${E};ft=${ft}'  # for wu
    # cli_args = 'Gc=${Gc};l=${l};E=${E}'             # for power
    execute_on = 'TIMESTEP_END'
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    multi_app = fracture
    direction = from_multiapp
    variable = d
    source_variable = d
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    multi_app = fracture
    direction = to_multiapp
    variable = 'psie_active'
    source_variable = 'psie_active'
  []
[]

[Mesh]
  type = FileMesh
  file = 'gold/geo.msh'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [T]
    initial_condition = 300.0
  []
[]

[AuxVariables]
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [fy]
  []
  [d]
  []
  [./j_x]
    order = FIRST
    family = MONOMIAL
  [../]
  [./j_y]
    order = FIRST
    family = MONOMIAL
  [../]
  [./dT_dx]
    order = FIRST
    family = MONOMIAL
  [../]
  [./dT_dy]
    order = FIRST
    family = MONOMIAL
  [../]
  [./k_J]
    order = FIRST
    family = MONOMIAL
  [../]
  [./vonmises]
    order = FIRST
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [psie_active]
    type = ADMaterialRealAux
    variable = psie_active
    property = psie_active
  []
  [./j_x]
    type = ParsedAux
    function = '- k_J * dT_dx'
    variable = j_x
    args = 'k_J dT_dx'
  [../]
  [./j_y]
    type = ParsedAux
    function = '- k_J * dT_dy'
    variable = j_y
    args = 'k_J  dT_dy'
  [../]
  [./k_J]
    type = ADMaterialRealAux
    property = k
    variable = k_J
  [../]
  [./grad_T_x]
    type = VariableGradientComponent
    gradient_variable = T
    variable = dT_dx
    component = x
  [../]
  [./grad_T_y]
    type = VariableGradientComponent
    gradient_variable = T
    variable = dT_dy
    component = y
  [../]
  [./vonmises]
    type = ADRankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = VonMisesStress
  [../]
[]

# [Modules/TensorMechanics/Master]
#   [all]
#     # This block adds all of the proper Kernels, strain calculators, and Variables
#     # for TensorMechanics in the correct coordinate system (autodetected)
#     add_variables = true
#     strain = SMALL
#     eigenstrain_names = eigenstrain
#     use_automatic_differentiation = true
#     generate_output = 'vonmises_stress' # elastic_strain_xx elastic_strain_yy strain_xx strain_yy'
#   []
# []

[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    save_in = fy
  []
  [heat_conduction]
    type = ADHeatConduction
    variable = T
    thermal_conductivity = k
    outputs = exodus
  []
  [heat_conduction_time_derivative]
    type = ADHeatConductionTimeDerivative
    variable = T
    density_name = rho
    specific_heat = c
  []
  [T_coupled_dot_d]
    type = ADParsedCoupledTimeDerivative
    variable = T
    v = d
    args = 'd psie_active T'
    coef_name = 'coef_1'
  []
  [T_coupled_Grad_dot_d]
    type = ADParsedCoupledGradTimeDerivative
    variable  = T
    v = d
    args = 'd T'
    coef_name = 'coef_2'
  []
  [T_coupled_stress]
    type = ElasticHeatEnergy
    variable = T
  []
[]

[BCs]
  [ydisp]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = 't'
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  # [xfix]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = bottom
  #   value = 0
  # []
  # [top_T]
  #   type = FunctionDirichletBC
  #   variable = T
  #   boundary = top
  #   function = top_T_load
  # []
  # [bottom_T]
  #   type = DirichletBC
  #   variable = T
  #   boundary = bottom
  #   value = 300
  # []
[]

[Functions]
  [top_T_load]
    type = ParsedFunction
    value =  'if((t>0) & (t< 250e-9), t*(2e8)+300.0, 350)'
  []
[]

[Materials]
  [coef_1]
    # type = ADDerivativeParsedMaterial
    # f_name = 'coef_1'
    # function = 'T * df_1'
    # args = 'T'
    # material_property_names = 'df_1:=D[f_1(d,T),T,d]'
    # derivative_order = 2
    # outputs = exodus
    type = ADParsedMaterial
    f_name = 'coef_1'
    material_property_names = 'df_1:=D[f_1(d,T),T,d]'
    function = 'T * df_1'
    args = 'T'
    outputs = exodus
  []
  [f_1]
    type = ADDerivativeParsedMaterial
    f_name = 'f_1'
    function = 'alpha * Gc / c0 / l + g * psie_active'
    args = 'T d psie_active'
    material_property_names = 'alpha(d) g(d) Gc(T) c0 l'
    derivative_order = 2
    outputs = exodus
  []
  [coef_2]
    type = ADDerivativeParsedMaterial
    f_name = 'coef_2'
    function = 'T * df_2'
    args = 'T'
    material_property_names = 'df_2:=D[f_2(T),T]'
    derivative_order = 1
    outputs = exodus
  []
  [f_2]
    type = ADDerivativeParsedMaterial
    f_name = 'f_2'
    function = '2*Gc*l/c0'
    args = 'T'
    material_property_names = 'Gc(T) c0 l'
    derivative_order = 1
    outputs = exodus
  []
  [dStress_dT*StrainRate]
    type = ComputeElasticHeatEnergy
    T = 'T'
    K = K
    G = G
    g = g
    # outputs = exodus
  []
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
    # outputs = exodus
  []
  [bulk]
    type = ADGenericConstantMaterial
    # prop_names = 'K G E ft l psic rho c k0 thermal_expansion'         # for wu
    # prop_values = '${K} ${G} ${E} ${ft} ${l} ${psic} ${rho} ${c} ${k0} ${alpha_bulk}'
    prop_names = 'K G E l rho c k0 thermal_expansion ft'                   # for power
    prop_values = '${K} ${G} ${E} ${l} ${rho} ${c} ${k0} ${alpha_bulk} ${ft}'
  []
  [Gc]
    type = ADDerivativeParsedMaterial
    f_name = 'Gc'
    function = '42.47e-3 * (1 - 1.8*((T-300)/1000) + 1.1 * ((T-300)/1000)^2 )'
    args = 'T'
    # material_property_names = 'T0'
    derivative_order = 1
    outputs = exodus
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
  [strain]
    type = ADComputeSmallStrain
    eigenstrain_names = eigenstrain
    # outputs = exodus
  []
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = SPECTRAL
    output_properties = 'elastic_strain psie_active '
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
  []
  [thermal_conductivity]
    type = ADDerivativeParsedMaterial
    f_name = 'k'
    function = 'k0 * ((1-d)^2 + 1e-8)'
    args = 'd '
    material_property_names = 'k0'
    outputs = exodus
  []
  [thermal_strain]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = 300
    eigenstrain_name = eigenstrain
    temperature = T
    thermal_expansion_coeff = ${alpha_bulk}
    # outputs = exodus
  []
[]

[Postprocessors]
  [Fy]
    type = NodalSum
    variable = fy
    boundary = top
  []
  [vonmises]
    type = ElementAverageValue
    variable = 'vonmises'
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '
  automatic_scaling = true

  dt = 1e-5
  end_time = 1e-3
  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-08

  picard_max_its = 20
  picard_abs_tol = 1e-50
  picard_rel_tol = 1e-03
  accept_on_max_picard_iteration = true
[]

[Outputs]
  exodus = true
  # interval = 10
  print_linear_residuals = false
  [csv]
    type = CSV
    delimiter = ' '
    file_base = 'f'
  []
[]
