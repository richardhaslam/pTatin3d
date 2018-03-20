
#include <ptatin3d.h>
#include <ptatin_init.h>
#include <quantity.h>

PetscErrorCode test_quantity_utils(void)
{
  PetscErrorCode ierr;
  PetscReal length_si,length_geo,length_model;

  PetscOptionsInsertString(NULL,"-ptatin_input_units geo");
  PetscOptionsInsertString(NULL,"-ptatin_q_length_mag 1000.0");
  ierr = ptatinQuantityCreate();CHKERRQ(ierr);

  /* take a quantity (length) and check it gets coverted correctly */
  length_si = 3000.0 * 1.0e3;
  length_model = pQConvert_SIUnits2ModelUnits(QLength,length_si);
  PetscPrintf(PETSC_COMM_WORLD,"length [nd]: %1.4e \n",length_model);

  length_model = 3.0;
  length_si = pQConvert_ModelUnits2SIUnits(QLength,length_model);
  PetscPrintf(PETSC_COMM_WORLD,"length si: %1.4e (m)\n",length_si);

  length_model = 3.0;
  length_geo = pQConvert_ModelUnits2GeoUnits(QLength,length_model);
  PetscPrintf(PETSC_COMM_WORLD,"length geo: %1.4e (km)\n",length_geo);

  length_geo = 3000.0;
  length_model = pQConvert_GeoUnits2ModelUnits(QLength,length_geo);
  PetscPrintf(PETSC_COMM_WORLD,"length [nd]: %1.4e (km)\n",length_model);

  ierr = ptatinQuantityDestroy();CHKERRQ(ierr);
  PetscOptionsClearValue(NULL,"-ptatin_input_units");

  PetscFunctionReturn(0);
}

PetscErrorCode test_quantity_utils_2(void)
{
  PetscErrorCode ierr;
  PetscReal length_si,length_geo,length_model;

  PetscOptionsInsertString(NULL,"-ptatin_input_units si");
  PetscOptionsInsertString(NULL,"-ptatin_q_length_mag 1000.0e3");
  ierr = ptatinQuantityCreate();CHKERRQ(ierr);

  /* take a quantity (length) and check it gets coverted correctly */
  length_si = 3000.0 * 1.0e3;
  length_model = pQConvert_SIUnits2ModelUnits(QLength,length_si);
  PetscPrintf(PETSC_COMM_WORLD,"length [nd]: %1.4e \n",length_model);

  length_model = 3.0;
  length_si = pQConvert_ModelUnits2SIUnits(QLength,length_model);
  PetscPrintf(PETSC_COMM_WORLD,"length si: %1.4e (m)\n",length_si);

  length_model = 3.0;
  length_geo = pQConvert_ModelUnits2GeoUnits(QLength,length_model);
  PetscPrintf(PETSC_COMM_WORLD,"length geo: %1.4e (km)\n",length_geo);

  length_geo = 3000.0;
  length_model = pQConvert_GeoUnits2ModelUnits(QLength,length_geo);
  PetscPrintf(PETSC_COMM_WORLD,"length [nd]: %1.4e (km)\n",length_model);

  ierr = ptatinQuantityDestroy();CHKERRQ(ierr);
  PetscOptionsClearValue(NULL,"-ptatin_input_units");

  PetscFunctionReturn(0);
}

PetscErrorCode test_quantity_utils_3(void)
{
  PetscErrorCode ierr;
  PetscReal visc_array[] = {1.0e18, 1.0e22, 1.0e24 };
  PetscReal visc_array2[3];

  PetscOptionsInsertString(NULL,"-ptatin_input_units si");
  PetscOptionsInsertString(NULL,"-ptatin_q_length_mag 1000.0e3");
  PetscOptionsInsertString(NULL,"-ptatin_q_velocity_mag 1.0e-10");
  PetscOptionsInsertString(NULL,"-ptatin_q_viscosity_mag 1.0e22");
  ierr = ptatinQuantityCreate();CHKERRQ(ierr);

  pQConvert_SIUnits2ModelUnits_Array(QViscosity,3,visc_array,visc_array2);
  printf("[0] visc %1.4e [-]\n",visc_array2[0]);
  printf("[1] visc %1.4e [-]\n",visc_array2[1]);
  printf("[2] visc %1.4e [-]\n",visc_array2[2]);

  pQConvert_SIUnits2ModelUnits_ArrayInPlace(QViscosity,3,visc_array);
  printf("[0] visc %1.4e [-]\n",visc_array[0]);
  printf("[1] visc %1.4e [-]\n",visc_array[1]);
  printf("[2] visc %1.4e [-]\n",visc_array[2]);

  pQConvert_ModelUnits2SIUnits_ArrayInPlace(QViscosity,3,visc_array);
  printf("[0] visc %1.4e [Pa s]\n",visc_array[0]);
  printf("[1] visc %1.4e [Pa s]\n",visc_array[1]);
  printf("[2] visc %1.4e [Pa s]\n",visc_array[2]);

  ierr = ptatinQuantityDestroy();CHKERRQ(ierr);
  PetscOptionsClearValue(NULL,"-ptatin_q_viscosity_mag");
  PetscOptionsClearValue(NULL,"-ptatin_q_velocity_mag");
  PetscOptionsClearValue(NULL,"-ptatin_q_length_mag");
  PetscOptionsClearValue(NULL,"-ptatin_input_units");

  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;

  ierr = pTatinInitialize(&argc,&argv,0,NULL);CHKERRQ(ierr);

  ierr = test_quantity_utils();CHKERRQ(ierr);
  ierr = test_quantity_utils_2();CHKERRQ(ierr);
  ierr = test_quantity_utils_3();CHKERRQ(ierr);

  ierr = pTatinFinalize();CHKERRQ(ierr);
  return 0;
}
