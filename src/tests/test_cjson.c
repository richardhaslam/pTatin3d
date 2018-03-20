/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   test_cjson.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

#include "petsc.h"
#include "ptatin3d.h"
#include "ptatin_init.h"
#include "geometry_object.h"
#include "geometry_object_parse.h"

PetscErrorCode test_GeometryObjectParse_cJSON(char filename[])
{
  PetscErrorCode  ierr;
  PetscInt        k,ngo;
  GeometryObject *golist;

  PetscPrintf(PETSC_COMM_WORLD,"** test_GeometryObjectParse_cJSON **\n");
  ierr = GeometryObjectLoadJSON(filename,&ngo,&golist);CHKERRQ(ierr);

  for (k=0; k<ngo; k++) {
    ierr = GeometryObjectView(golist[k]);CHKERRQ(ierr);
  }

  /* Gnuplot test viewer */
  /*
     {
     FILE *fp;
     int i,j,k,nn = 60;
     double dx,pos[3];

     fp = fopen("geom.gp","w");
     dx = 3.0/((double)nn-1.0);
     for (k=0; k<nn; k++) {
     for (j=0; j<nn; j++) {
     for (i=0; i<nn; i++) {
     int inside;

     pos[0] = i * dx;
     pos[1] = j * dx;
     pos[2] = k * dx;

     ierr = GeometryObjectPointInside(golist[0],pos,&inside);CHKERRQ(ierr);

     if (inside == 1) {
     fprintf(fp,"%1.4e %1.4e %1.4e\n",pos[0],pos[1],pos[2]);
     }
     }
     }
     }
     fclose(fp);
     }
     */

  for (k=0; k<ngo; k++) {
    ierr = GeometryObjectDestroy(&golist[k]);CHKERRQ(ierr);
  }
  free(golist);

  PetscFunctionReturn(0);
}

PetscErrorCode test_GeometryObjectParse2_cJSON(void)
{
  PetscErrorCode ierr;
  PetscInt       k,ngo;
  GeometryObject *golist;
  FILE           *fp;
  PetscLogDouble t0,t1;
  const char     filename[] = "fat_test_geom.json";

  PetscPrintf(PETSC_COMM_WORLD,"** test_GeometryObjectParse2_cJSON **\n");

  ngo = 100;
  PetscOptionsGetInt(NULL,NULL,"-ngo",&ngo,NULL);
  fp = fopen(filename,"w");
  fprintf(fp,"{\n");

  fprintf(fp,"\"GeometryObjectList\": [\n");
  for (k=0; k<ngo; k++) {

    fprintf(fp,"{\n");
    fprintf(fp,"    \"name\":     \"Ucrust-%d\",\n",k+1);
    fprintf(fp,"    \"type\":     \"GeomType_EllipticCylinder\",\n");
    fprintf(fp,"    \"centroid\": [1.5, 1.2, 1.3],\n");
    fprintf(fp,"    \"length\":   1.0,\n");
    fprintf(fp,"    \"radiusA\":  0.6,\n");
    fprintf(fp,"    \"radiusB\":  0.8,\n");
    fprintf(fp,"    \"axis\":     \"Y\"\n");
    if (k <= ngo-2) {
      fprintf(fp,"},\n");
    } else {
      fprintf(fp,"}\n");
    }
  }

  fprintf(fp,"]\n");

  fprintf(fp,"}\n");
  fclose(fp);

  PetscTime(&t0);
  ierr = GeometryObjectLoadJSON(filename,&ngo,&golist);CHKERRQ(ierr);
  PetscTime(&t1);
  PetscPrintf(PETSC_COMM_WORLD,"Time to parse %D geom objects: %1.4e (sec)\n",ngo,t1-t0);

  for (k=0; k<ngo; k++) {
    ierr = GeometryObjectDestroy(&golist[k]);CHKERRQ(ierr);
  }
  free(golist);

  PetscFunctionReturn(0);
}

PetscErrorCode render(void)
{
  char filename[] = "ptatin_mesh_state.json";
  FILE *file;
  cJSON *root,*parent,*item;
  char *j;

  root = cJSON_CreateObject();

  parent = cJSON_CreateObject();
    item = cJSON_CreateNumber(64);  cJSON_AddItemToObject(parent,"mx",item);
    item = cJSON_CreateNumber(32);  cJSON_AddItemToObject(parent,"my",item);
    item = cJSON_CreateNumber(64);  cJSON_AddItemToObject(parent,"mz",item);
  cJSON_AddItemToObject(root,"StokesDMDA",parent);

  item = cJSON_CreateString("rift3D");  cJSON_AddItemToObject(root,"ptatinModel",item);


  j = cJSON_Print(root);
  printf("[ptatin state]\n");
  printf("%s\n",j);
  file = fopen(filename,"w");
  fprintf(file,"%s",j);
  fclose(file);
  free(j);

  cJSON_Delete(root);

  PetscFunctionReturn(0);
}

PetscErrorCode render_units(void)
{
  char filename[] = "ptatin_units.json";
  FILE *file;
  cJSON *root,*parent,*array,*item;
  char *j;

  root = cJSON_CreateObject();

  {
    const char *units[] = { "length", "velocity", "viscosity" };
    item = cJSON_CreateStringArray(units,3);  cJSON_AddItemToObject(root,"ptatinReferenceScales",item);
  }

  array = cJSON_CreateArray();
  // length
  parent = cJSON_CreateObject();
  item = cJSON_CreateString("length");    cJSON_AddItemToObject(parent,"quantity",item);
  item = cJSON_CreateString("m");         cJSON_AddItemToObject(parent,"unit",item);
  item = cJSON_CreateNumber(1.0e5);       cJSON_AddItemToObject(parent,"scale",item);
  cJSON_AddItemToArray(array,parent);
  // velocity
  parent = cJSON_CreateObject();
  item = cJSON_CreateString("velocity");    cJSON_AddItemToObject(parent,"quantity",item);
  item = cJSON_CreateString("m/s");         cJSON_AddItemToObject(parent,"unit",item);
  item = cJSON_CreateNumber(1.0e-10);       cJSON_AddItemToObject(parent,"scale",item);
  cJSON_AddItemToArray(array,parent);

  cJSON_AddItemToObject(root,"ptatinScales",array);


  j = cJSON_Print(root);
  printf("[ptatin units]\n");
  printf("%s\n",j);
  file = fopen(filename,"w");
  fprintf(file,"%s",j);
  fclose(file);
  free(j);

  cJSON_Delete(root);

  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscInt       i,test_id,nloops;
  char           filename1[PETSC_MAX_PATH_LEN] = "src/tests/test_geom.json";

  ierr = pTatinInitialize(&argc,&argv,(char *)0,NULL);CHKERRQ(ierr);

  test_id = 1;
  PetscOptionsGetInt(NULL,NULL,"-test_id",&test_id,NULL);
  PetscOptionsGetString(NULL,NULL,"-filename1",filename1,PETSC_MAX_PATH_LEN,NULL);

  switch (test_id) {

    case 1:
      nloops = 1;
      PetscOptionsGetInt(NULL,NULL,"-nloops",&nloops,NULL);
      for (i=0; i<nloops; i++) {
        ierr = test_GeometryObjectParse_cJSON(filename1);CHKERRQ(ierr);
      }
      break;

    case 2:
      ierr = test_GeometryObjectParse2_cJSON();CHKERRQ(ierr);
      break;

    case 3:
      ierr = render();CHKERRQ(ierr);
      ierr = render_units();CHKERRQ(ierr);
      break;

    default:
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Undefined test id %d\n",test_id);
      break;
    }

  ierr = pTatinFinalize();CHKERRQ(ierr);
  return 0;
}
