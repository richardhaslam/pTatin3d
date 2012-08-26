
#ifndef _ptatin_units_h__
#define _ptatin_units_h__

typedef struct _p_Units *Units;

struct _p_Units {
	char   *quantity_name;
	char   *unit_name;
	double characteristic_scale;
};


void UnitsApplyScaling(Units u,double y,double *ystar);
void UnitsApplyScalingToArray(Units u,int n,double y[],double ystar[]);
void UnitsApplyInverseScaling(Units u,double y,double *ystar);
void UnitsApplyInverseScalingToArray(Units u,int n,double y[],double ystar[]);

void UnitsConvertSI2Kilo(double si,double *ksi);
void UnitsConvertSI2Mega(double si,double *ksi);
void UnitsConvertSI2Giga(double si,double *ksi);
void UnitsConvertSI2Milli(double si,double *ksi);
void UnitsConvertSI2Micro(double si,double *ksi);
void UnitsConvertSI2Pico(double si,double *ksi);

void UnitsCreate(Units *_u,const char quantity_name[],const char unit_name[],double scale);
void UnitsSetScale(Units u,double scale);
void UnitsView(Units u);

void UnitsCreateSIVelocity(Units *u);
void UnitsCreateSILength(Units *u);
void UnitsCreateSITime(Units *u);
void UnitsCreateSIStress(Units *u);
void UnitsCreateSIViscosity(Units *u);
void UnitsCreateSITemperature(Units *u);

#endif
