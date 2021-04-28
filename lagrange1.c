// This code conforms with the UFC specification version 2018.2.0.dev0
// and was automatically generated by FFCX version 0.1.0.
//
// This code was generated with the following parameters:
//
//  {'assume_aligned': -1,
//   'epsilon': 1e-14,
//   'output_directory': '.',
//   'padlen': 1,
//   'profile': False,
//   'scalar_type': 'double',
//   'table_atol': 1e-09,
//   'table_rtol': 1e-06,
//   'tabulate_tensor_void': False,
//   'ufl_file': ['lagrange.ufl'],
//   'verbosity': 30,
//   'visualise': False}

typedef double ufc_scalar_t;
#include <math.h>
#include <stdalign.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <ufc.h>

// Code for element element_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4

ufc_finite_element element_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4
    = {.signature = "FiniteElement('Lagrange', tetrahedron, 1)",
       .cell_shape = tetrahedron,
       .topological_dimension = 3,
       .geometric_dimension = 3,
       .space_dimension = 4,
       .value_rank = 0,
       .value_shape = NULL,
       .value_size = 1,
       .reference_value_rank = 0,
       .reference_value_shape = NULL,
       .reference_value_size = 1,
       .degree = 1,
       .family = "Lagrange",
       .block_size = 1,

       .needs_transformation_data = 0,
       .interpolation_is_identity = 1,

       .num_sub_elements = 0,
       .sub_elements = NULL};

// End of code for element element_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4

// Code for element element_f98876eb49550d506341b4c1d32d16b6eec2e3aa

int value_shape_element_f98876eb49550d506341b4c1d32d16b6eec2e3aa[1] = {3};
int reference_value_shape_element_f98876eb49550d506341b4c1d32d16b6eec2e3aa[1] = {3};
ufc_finite_element* sub_elements_element_f98876eb49550d506341b4c1d32d16b6eec2e3aa[3]
    = {&element_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
       &element_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
       &element_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4};

ufc_finite_element element_f98876eb49550d506341b4c1d32d16b6eec2e3aa = {
    .signature = "VectorElement(FiniteElement('Lagrange', tetrahedron, 1), dim=3)",
    .cell_shape = tetrahedron,
    .topological_dimension = 3,
    .geometric_dimension = 3,
    .space_dimension = 12,
    .value_rank = 1,
    .value_shape = value_shape_element_f98876eb49550d506341b4c1d32d16b6eec2e3aa,
    .value_size = 3,
    .reference_value_rank = 1,
    .reference_value_shape = reference_value_shape_element_f98876eb49550d506341b4c1d32d16b6eec2e3aa,
    .reference_value_size = 3,
    .degree = 1,
    .family = "Lagrange",
    .block_size = 3,

    .needs_transformation_data = 0,
    .interpolation_is_identity = 1,

    .num_sub_elements = 3,
    .sub_elements = sub_elements_element_f98876eb49550d506341b4c1d32d16b6eec2e3aa};

// End of code for element element_f98876eb49550d506341b4c1d32d16b6eec2e3aa

// Code for dofmap dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4

void tabulate_entity_dofs_dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4(int* restrict dofs, int d,
                                                                          int i)
{
  switch (d)
  {
  case 0:
    switch (i)
    {
    case 0:
      dofs[0] = 0;
      break;
    case 1:
      dofs[0] = 1;
      break;
    case 2:
      dofs[0] = 2;
      break;
    case 3:
      dofs[0] = 3;
      break;
    }
    break;
  }
}

ufc_dofmap dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4
    = {.signature = "FFCX dofmap for FiniteElement('Lagrange', tetrahedron, 1)",
       .block_size = 1,
       .num_global_support_dofs = 0,
       .num_element_support_dofs = 4,
       .num_entity_dofs[0] = 1,
       .num_entity_dofs[1] = 0,
       .num_entity_dofs[2] = 0,
       .num_entity_dofs[3] = 0,
       .tabulate_entity_dofs = tabulate_entity_dofs_dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
       .num_sub_dofmaps = 0,
       .sub_dofmaps = NULL};

// End of code for dofmap dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4

// Code for dofmap dofmap_f98876eb49550d506341b4c1d32d16b6eec2e3aa

ufc_dofmap* sub_dofmaps_dofmap_f98876eb49550d506341b4c1d32d16b6eec2e3aa[3]
    = {&dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
       &dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
       &dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4};

void tabulate_entity_dofs_dofmap_f98876eb49550d506341b4c1d32d16b6eec2e3aa(int* restrict dofs, int d,
                                                                          int i)
{
  switch (d)
  {
  case 0:
    switch (i)
    {
    case 0:
      dofs[0] = 0;
      break;
    case 1:
      dofs[0] = 1;
      break;
    case 2:
      dofs[0] = 2;
      break;
    case 3:
      dofs[0] = 3;
      break;
    }
    break;
  }
}

ufc_dofmap dofmap_f98876eb49550d506341b4c1d32d16b6eec2e3aa = {
    .signature = "FFCX dofmap for VectorElement(FiniteElement('Lagrange', tetrahedron, 1), dim=3)",
    .block_size = 3,
    .num_global_support_dofs = 0,
    .num_element_support_dofs = 4,
    .num_entity_dofs[0] = 1,
    .num_entity_dofs[1] = 0,
    .num_entity_dofs[2] = 0,
    .num_entity_dofs[3] = 0,
    .tabulate_entity_dofs = tabulate_entity_dofs_dofmap_f98876eb49550d506341b4c1d32d16b6eec2e3aa,
    .num_sub_dofmaps = 3,
    .sub_dofmaps = sub_dofmaps_dofmap_f98876eb49550d506341b4c1d32d16b6eec2e3aa};

// End of code for dofmap dofmap_f98876eb49550d506341b4c1d32d16b6eec2e3aa

// Code for integral integral_cell_otherwise_6169f1b8b02a174df5aa5b2135c811df067e58d5

void tabulate_tensor_integral_cell_otherwise_6169f1b8b02a174df5aa5b2135c811df067e58d5(
    ufc_scalar_t* restrict A, const ufc_scalar_t* restrict w, const ufc_scalar_t* restrict c,
    const double* restrict coordinate_dofs, const int* restrict entity_local_index,
    const uint8_t* restrict quadrature_permutation, const uint32_t cell_permutation)
{
  // Quadrature rules
  static const double weights_421[1] = {0.1666666666666667};
  // Precomputed values of basis functions and precomputations
  // FE* dimensions: [permutation][entities][points][dofs]
  static const double FE8_C0_D001_Q421[1][1][1][4] = {{{{-1.0, 0.0, 0.0, 1.0}}}};
  static const double FE8_C0_D010_Q421[1][1][1][4] = {{{{-1.0, 0.0, 1.0, 0.0}}}};
  static const double FE8_C0_D100_Q421[1][1][1][4] = {{{{-1.0, 1.0, 0.0, 0.0}}}};
  // Quadrature loop independent computations for quadrature rule 421
  const double J_c4 = coordinate_dofs[1] * FE8_C0_D010_Q421[0][0][0][0]
                      + coordinate_dofs[4] * FE8_C0_D010_Q421[0][0][0][1]
                      + coordinate_dofs[7] * FE8_C0_D010_Q421[0][0][0][2]
                      + coordinate_dofs[10] * FE8_C0_D010_Q421[0][0][0][3];
  const double J_c8 = coordinate_dofs[2] * FE8_C0_D001_Q421[0][0][0][0]
                      + coordinate_dofs[5] * FE8_C0_D001_Q421[0][0][0][1]
                      + coordinate_dofs[8] * FE8_C0_D001_Q421[0][0][0][2]
                      + coordinate_dofs[11] * FE8_C0_D001_Q421[0][0][0][3];
  const double J_c5 = coordinate_dofs[1] * FE8_C0_D001_Q421[0][0][0][0]
                      + coordinate_dofs[4] * FE8_C0_D001_Q421[0][0][0][1]
                      + coordinate_dofs[7] * FE8_C0_D001_Q421[0][0][0][2]
                      + coordinate_dofs[10] * FE8_C0_D001_Q421[0][0][0][3];
  const double J_c7 = coordinate_dofs[2] * FE8_C0_D010_Q421[0][0][0][0]
                      + coordinate_dofs[5] * FE8_C0_D010_Q421[0][0][0][1]
                      + coordinate_dofs[8] * FE8_C0_D010_Q421[0][0][0][2]
                      + coordinate_dofs[11] * FE8_C0_D010_Q421[0][0][0][3];
  const double J_c0 = coordinate_dofs[0] * FE8_C0_D100_Q421[0][0][0][0]
                      + coordinate_dofs[3] * FE8_C0_D100_Q421[0][0][0][1]
                      + coordinate_dofs[6] * FE8_C0_D100_Q421[0][0][0][2]
                      + coordinate_dofs[9] * FE8_C0_D100_Q421[0][0][0][3];
  const double J_c1 = coordinate_dofs[0] * FE8_C0_D010_Q421[0][0][0][0]
                      + coordinate_dofs[3] * FE8_C0_D010_Q421[0][0][0][1]
                      + coordinate_dofs[6] * FE8_C0_D010_Q421[0][0][0][2]
                      + coordinate_dofs[9] * FE8_C0_D010_Q421[0][0][0][3];
  const double J_c6 = coordinate_dofs[2] * FE8_C0_D100_Q421[0][0][0][0]
                      + coordinate_dofs[5] * FE8_C0_D100_Q421[0][0][0][1]
                      + coordinate_dofs[8] * FE8_C0_D100_Q421[0][0][0][2]
                      + coordinate_dofs[11] * FE8_C0_D100_Q421[0][0][0][3];
  const double J_c3 = coordinate_dofs[1] * FE8_C0_D100_Q421[0][0][0][0]
                      + coordinate_dofs[4] * FE8_C0_D100_Q421[0][0][0][1]
                      + coordinate_dofs[7] * FE8_C0_D100_Q421[0][0][0][2]
                      + coordinate_dofs[10] * FE8_C0_D100_Q421[0][0][0][3];
  const double J_c2 = coordinate_dofs[0] * FE8_C0_D001_Q421[0][0][0][0]
                      + coordinate_dofs[3] * FE8_C0_D001_Q421[0][0][0][1]
                      + coordinate_dofs[6] * FE8_C0_D001_Q421[0][0][0][2]
                      + coordinate_dofs[9] * FE8_C0_D001_Q421[0][0][0][3];
  ufc_scalar_t sp_421[80];
  sp_421[0] = J_c4 * J_c8;
  sp_421[1] = J_c5 * J_c7;
  sp_421[2] = sp_421[0] + -1 * sp_421[1];
  sp_421[3] = J_c0 * sp_421[2];
  sp_421[4] = J_c5 * J_c6;
  sp_421[5] = J_c3 * J_c8;
  sp_421[6] = sp_421[4] + -1 * sp_421[5];
  sp_421[7] = J_c1 * sp_421[6];
  sp_421[8] = sp_421[3] + sp_421[7];
  sp_421[9] = J_c3 * J_c7;
  sp_421[10] = J_c4 * J_c6;
  sp_421[11] = sp_421[9] + -1 * sp_421[10];
  sp_421[12] = J_c2 * sp_421[11];
  sp_421[13] = sp_421[8] + sp_421[12];
  sp_421[14] = sp_421[2] / sp_421[13];
  sp_421[15] = J_c3 * (-1 * J_c8);
  sp_421[16] = sp_421[4] + sp_421[15];
  sp_421[17] = sp_421[16] / sp_421[13];
  sp_421[18] = sp_421[11] / sp_421[13];
  sp_421[19] = sp_421[14] * sp_421[14];
  sp_421[20] = sp_421[14] * sp_421[17];
  sp_421[21] = sp_421[18] * sp_421[14];
  sp_421[22] = sp_421[17] * sp_421[17];
  sp_421[23] = sp_421[18] * sp_421[17];
  sp_421[24] = sp_421[18] * sp_421[18];
  sp_421[25] = J_c2 * J_c7;
  sp_421[26] = J_c8 * (-1 * J_c1);
  sp_421[27] = sp_421[25] + sp_421[26];
  sp_421[28] = sp_421[27] / sp_421[13];
  sp_421[29] = J_c0 * J_c8;
  sp_421[30] = J_c6 * (-1 * J_c2);
  sp_421[31] = sp_421[29] + sp_421[30];
  sp_421[32] = sp_421[31] / sp_421[13];
  sp_421[33] = J_c1 * J_c6;
  sp_421[34] = J_c0 * J_c7;
  sp_421[35] = sp_421[33] + -1 * sp_421[34];
  sp_421[36] = sp_421[35] / sp_421[13];
  sp_421[37] = sp_421[28] * sp_421[28];
  sp_421[38] = sp_421[28] * sp_421[32];
  sp_421[39] = sp_421[28] * sp_421[36];
  sp_421[40] = sp_421[32] * sp_421[32];
  sp_421[41] = sp_421[32] * sp_421[36];
  sp_421[42] = sp_421[36] * sp_421[36];
  sp_421[43] = sp_421[37] + sp_421[19];
  sp_421[44] = sp_421[38] + sp_421[20];
  sp_421[45] = sp_421[39] + sp_421[21];
  sp_421[46] = sp_421[40] + sp_421[22];
  sp_421[47] = sp_421[41] + sp_421[23];
  sp_421[48] = sp_421[24] + sp_421[42];
  sp_421[49] = J_c1 * J_c5;
  sp_421[50] = J_c2 * J_c4;
  sp_421[51] = sp_421[49] + -1 * sp_421[50];
  sp_421[52] = sp_421[51] / sp_421[13];
  sp_421[53] = J_c2 * J_c3;
  sp_421[54] = J_c0 * J_c5;
  sp_421[55] = sp_421[53] + -1 * sp_421[54];
  sp_421[56] = sp_421[55] / sp_421[13];
  sp_421[57] = J_c0 * J_c4;
  sp_421[58] = J_c1 * J_c3;
  sp_421[59] = sp_421[57] + -1 * sp_421[58];
  sp_421[60] = sp_421[59] / sp_421[13];
  sp_421[61] = sp_421[52] * sp_421[52];
  sp_421[62] = sp_421[52] * sp_421[56];
  sp_421[63] = sp_421[60] * sp_421[52];
  sp_421[64] = sp_421[56] * sp_421[56];
  sp_421[65] = sp_421[60] * sp_421[56];
  sp_421[66] = sp_421[60] * sp_421[60];
  sp_421[67] = sp_421[43] + sp_421[61];
  sp_421[68] = sp_421[44] + sp_421[62];
  sp_421[69] = sp_421[45] + sp_421[63];
  sp_421[70] = sp_421[46] + sp_421[64];
  sp_421[71] = sp_421[47] + sp_421[65];
  sp_421[72] = sp_421[48] + sp_421[66];
  sp_421[73] = fabs(sp_421[13]);
  sp_421[74] = sp_421[67] * sp_421[73];
  sp_421[75] = sp_421[68] * sp_421[73];
  sp_421[76] = sp_421[69] * sp_421[73];
  sp_421[77] = sp_421[70] * sp_421[73];
  sp_421[78] = sp_421[71] * sp_421[73];
  sp_421[79] = sp_421[72] * sp_421[73];
  for (int iq = 0; iq < 1; ++iq)
  {
    const ufc_scalar_t fw0 = sp_421[74] * weights_421[iq];
    const ufc_scalar_t fw1 = sp_421[75] * weights_421[iq];
    const ufc_scalar_t fw2 = sp_421[76] * weights_421[iq];
    const ufc_scalar_t fw3 = sp_421[77] * weights_421[iq];
    const ufc_scalar_t fw4 = sp_421[78] * weights_421[iq];
    const ufc_scalar_t fw5 = sp_421[79] * weights_421[iq];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
      {
        A[4 * i + j] += fw0 * FE8_C0_D100_Q421[0][0][0][i] * FE8_C0_D100_Q421[0][0][0][j];
        A[4 * i + j] += fw1 * FE8_C0_D100_Q421[0][0][0][i] * FE8_C0_D010_Q421[0][0][0][j];
        A[4 * i + j] += fw2 * FE8_C0_D100_Q421[0][0][0][i] * FE8_C0_D001_Q421[0][0][0][j];
        A[4 * i + j] += fw1 * FE8_C0_D010_Q421[0][0][0][i] * FE8_C0_D100_Q421[0][0][0][j];
        A[4 * i + j] += fw3 * FE8_C0_D010_Q421[0][0][0][i] * FE8_C0_D010_Q421[0][0][0][j];
        A[4 * i + j] += fw4 * FE8_C0_D010_Q421[0][0][0][i] * FE8_C0_D001_Q421[0][0][0][j];
        A[4 * i + j] += fw2 * FE8_C0_D001_Q421[0][0][0][i] * FE8_C0_D100_Q421[0][0][0][j];
        A[4 * i + j] += fw4 * FE8_C0_D001_Q421[0][0][0][i] * FE8_C0_D010_Q421[0][0][0][j];
        A[4 * i + j] += fw5 * FE8_C0_D001_Q421[0][0][0][i] * FE8_C0_D001_Q421[0][0][0][j];
      }
  }
}

ufc_integral integral_cell_otherwise_6169f1b8b02a174df5aa5b2135c811df067e58d5
    = {.enabled_coefficients = NULL,
       .tabulate_tensor
       = tabulate_tensor_integral_cell_otherwise_6169f1b8b02a174df5aa5b2135c811df067e58d5,
       .needs_transformation_data = 0};

// End of code for integral integral_cell_otherwise_6169f1b8b02a174df5aa5b2135c811df067e58d5

// Code for form form_6169f1b8b02a174df5aa5b2135c811df067e58d5

ufc_dofmap* dofmaps_form_6169f1b8b02a174df5aa5b2135c811df067e58d5[2]
    = {&dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
       &dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4};
ufc_finite_element* finite_elements_form_6169f1b8b02a174df5aa5b2135c811df067e58d5[2]
    = {&element_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
       &element_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4};

// Return a list of the coefficient names.
const char** coefficient_name_form_6169f1b8b02a174df5aa5b2135c811df067e58d5(void) { return NULL; }

// Return a list of the constant names.
const char** constant_name_form_6169f1b8b02a174df5aa5b2135c811df067e58d5(void) { return NULL; }

int* integral_ids_form_6169f1b8b02a174df5aa5b2135c811df067e58d5(ufc_integral_type integral_type)
{
  static int integral_ids_cell_form_6169f1b8b02a174df5aa5b2135c811df067e58d5[1] = {-1};
  switch (integral_type)
  {
  case cell:
    return integral_ids_cell_form_6169f1b8b02a174df5aa5b2135c811df067e58d5;
  default:
    return NULL;
  }
}

int num_integrals_form_6169f1b8b02a174df5aa5b2135c811df067e58d5(ufc_integral_type integral_type)
{
  switch (integral_type)
  {
  case cell:
    return 1;
  case interior_facet:
    return 0;
  case exterior_facet:
    return 0;
  default:
    return 0;
  }
}

ufc_integral**
integrals_form_6169f1b8b02a174df5aa5b2135c811df067e58d5(ufc_integral_type integral_type)
{
  static ufc_integral* integrals_cell_form_6169f1b8b02a174df5aa5b2135c811df067e58d5[1]
      = {&integral_cell_otherwise_6169f1b8b02a174df5aa5b2135c811df067e58d5};
  switch (integral_type)
  {
  case cell:
    return integrals_cell_form_6169f1b8b02a174df5aa5b2135c811df067e58d5;
  default:
    return NULL;
  }
}

ufc_form form_6169f1b8b02a174df5aa5b2135c811df067e58d5 = {

    .signature = "a4bfe56d456282b252665ff416d395c3487df4a19188830c511c169c360cd953886e4a17e40c4085b"
                 "0b3c619474987b1b91f13b4eca9859bee0dfe00da9d3b48",
    .rank = 2,
    .num_coefficients = 0,
    .num_constants = 0,
    .original_coefficient_position = NULL,

    .coefficient_name_map = coefficient_name_form_6169f1b8b02a174df5aa5b2135c811df067e58d5,
    .constant_name_map = constant_name_form_6169f1b8b02a174df5aa5b2135c811df067e58d5,

    .finite_elements = finite_elements_form_6169f1b8b02a174df5aa5b2135c811df067e58d5,
    .dofmaps = dofmaps_form_6169f1b8b02a174df5aa5b2135c811df067e58d5,

    .integral_ids = integral_ids_form_6169f1b8b02a174df5aa5b2135c811df067e58d5,
    .num_integrals = num_integrals_form_6169f1b8b02a174df5aa5b2135c811df067e58d5,

    .integrals = integrals_form_6169f1b8b02a174df5aa5b2135c811df067e58d5};

// Alias name
ufc_form* form_lagrange_a = &form_6169f1b8b02a174df5aa5b2135c811df067e58d5;

ufc_function_space* functionspace_form_lagrange_a(const char* function_name)
{
  static ufc_function_space functionspace_v
      = {.finite_element = &element_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
         .dofmap = &dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
         .geometry_family = "Lagrange",
         .geometry_degree = 1};
  static ufc_function_space functionspace_u
      = {.finite_element = &element_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
         .dofmap = &dofmap_0867dce6bf5c3352d4f192cc2fe32c1c51c671a4,
         .geometry_family = "Lagrange",
         .geometry_degree = 1};
  if (strcmp(function_name, "v") == 0)
    return &functionspace_v;
  else if (strcmp(function_name, "u") == 0)
    return &functionspace_u;
  return NULL;
}

// End of code for form form_6169f1b8b02a174df5aa5b2135c811df067e58d5
